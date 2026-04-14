#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ----------------------------
# User-configurable paths
# ----------------------------
ARAB_FILE = Path("/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS_annotate/arabidopsis_TN_0_5_mean.annot_25000bp.tsv")
RICE_FILE = Path("/Users/nirwantandukar/Documents/Research/results/Rice_3001/GWAS_annotate/rice_N_mod_sub_rice_gwas_phenotype_TN_rice3000_gwas_qc.snp.assoc.annot_25000bp.tsv")
MAIZE_FILE = Path("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/annotation_maize_LMM_nitrogen_0_5.csv")
SORGHUM_FILE = Path("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/annotation_sorghum_LMM_nitrogen_0-5_sorghum.csv")
BARLEY_FILE = Path("/Users/nirwantandukar/Documents/Research/results/Barley/GWAS_annotate/barley_soilN_mod_sub_barley_soilN_gwas_phenotype_gbs_landrace_12129inds_beagle_imputed_isec558kSNP.chrfix.assoc.annot_25000bp.tsv")
MAIZE_RAW_FILE = Path("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_maize_LMM.txt")
SORGHUM_RAW_FILE = Path("/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/nitrogen/GWAS_results/nitrogen_0-5cm_sorghum_LMM.txt")

OG2GENES_FILE = Path("/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_OG2genes.tab.gz")
GENES_FILE = Path("/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_genes.tab.gz")

OUT_DIR = Path("/Users/nirwantandukar/Documents/Github/Landadapt/results")
P_THRESHOLD = 1e-3
TOP_FRACTION: float | None = None
PREFILTER_THRESHOLD: float | None = None
USE_SNP_RANK_CUTOFF: bool = False
ALLOWED_TAXA = {"3702_0", "4577_0", "4558_0", "4557_0", "39947_0", "4530_0", "4513_0", "112509_0"}


def p_label(p: float) -> str:
    s = f"{p:.0e}".replace("+0", "").replace("+", "")
    return s.lower()


def top_label(f: float) -> str:
    pct = f * 100.0
    if abs(pct - round(pct)) < 1e-9:
        return f"{int(round(pct))}pct"
    return f"{pct:g}pct"


def normalize_gene_id(gene_id: str) -> str:
    g = str(gene_id).strip().upper()
    g = re.sub(r"^(GENE:|TRANSCRIPT:)", "", g)
    g = re.sub(r"_T\d+$", "", g)
    g = re.sub(r"_P\d+$", "", g)
    g = re.sub(r"\.[0-9]+$", "", g)
    g = re.sub(r"-[A-Z]{1,3}$", "", g)
    return g


def build_aliases(species: str, gene_id: str) -> set[str]:
    aliases = {str(gene_id).strip().upper()}
    norm = normalize_gene_id(gene_id)
    aliases.add(norm)

    if species == "Rice":
        m = re.match(r"^LOC_OS(\d{2})G(\d+)$", norm)
        if m:
            chrom, num = m.groups()
            aliases.add(f"OS{chrom}G{num}")
            aliases.add(f"OS{chrom}G{num}0")
            aliases.add(f"OS{chrom}G{num}00")
            aliases.add(f"OS{chrom}G{num.zfill(7)}")

    return aliases


def compute_species_snp_rank_cutoffs(top_fraction: float) -> tuple[dict[str, float], pd.DataFrame]:
    snp_cfg = [
        ("Arabidopsis", ARAB_FILE, "\t", "CHR", "BP", "P"),
        ("Rice", RICE_FILE, "\t", "chr", "ps", "p_wald"),
        ("Barley", BARLEY_FILE, "\t", "chr", "ps", "p_wald"),
        ("Maize", MAIZE_RAW_FILE, "\t", "chr", "ps", "p_wald"),
        ("Sorghum", SORGHUM_RAW_FILE, "\t", "chr", "ps", "p_wald"),
    ]

    cutoffs: dict[str, float] = {}
    rows = []

    for species, file, sep, chr_col, pos_col, p_col in snp_cfg:
        df = pd.read_csv(file, sep=sep, usecols=[chr_col, pos_col, p_col], low_memory=False)
        df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
        df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
        out = df.dropna(subset=[chr_col, pos_col, p_col]).copy()
        out = out[out[p_col] > 0]
        out[chr_col] = out[chr_col].astype(str).str.strip()
        out = out[out[chr_col] != ""]

        out = out.groupby([chr_col, pos_col], as_index=False)[p_col].min()

        n_total = len(out)
        n_rank = max(1, int(n_total * top_fraction)) if n_total > 0 else 0
        if n_rank > 0:
            p_cutoff = float(out.nsmallest(n_rank, p_col).iloc[-1][p_col])
            n_selected = int((out[p_col] <= p_cutoff).sum())
        else:
            p_cutoff = float("nan")
            n_selected = 0

        cutoffs[species] = p_cutoff
        rows.append(
            {
                "species": species,
                "snp_rows_total": int(n_total),
                "snp_rank_target": int(n_rank),
                "snp_selected_by_cutoff": int(n_selected),
                "p_cutoff": p_cutoff,
            }
        )

    return cutoffs, pd.DataFrame(rows)


def load_significant_genes(
    species_p_cutoff: dict[str, float] | None = None,
) -> tuple[dict[str, set[str]], dict[str, dict[str, float]], dict[str, dict[str, int]]]:
    species_genes: dict[str, set[str]] = {}
    species_gene_maxp: dict[str, dict[str, float]] = {}
    species_row_stats: dict[str, dict[str, int]] = {}

    def pick_sig(df: pd.DataFrame, gene_col: str, p_col: str, species_name: str) -> pd.DataFrame:
        out = df[[gene_col, p_col]].copy()
        out[p_col] = pd.to_numeric(out[p_col], errors="coerce")
        out = out.dropna(subset=[gene_col, p_col])
        out[gene_col] = out[gene_col].astype(str).str.strip()
        out = out[out[gene_col] != ""]
        n_total = len(out)

        if species_p_cutoff is not None:
            cutoff = species_p_cutoff.get(species_name)
            if cutoff is None or pd.isna(cutoff):
                out = out.iloc[0:0]
            else:
                out = out[out[p_col] <= cutoff]
            n_pref = len(out)
            n_keep = len(out)
        elif TOP_FRACTION is None:
            out = out[out[p_col] <= P_THRESHOLD]
            n_pref = len(out)
            n_keep = len(out)
        else:
            if PREFILTER_THRESHOLD is not None:
                out = out[out[p_col] <= PREFILTER_THRESHOLD]
            n_pref = len(out)
            n_keep = max(1, int(n_pref * TOP_FRACTION)) if n_pref > 0 else 0
            out = out.nsmallest(n_keep, p_col, keep="all") if n_keep > 0 else out.iloc[0:0]

        species_row_stats[species_name] = {
            "rows_total": int(n_total),
            "rows_after_prefilter": int(n_pref),
            "rows_selected": int(len(out)),
        }
        return out

    arab = pd.read_csv(ARAB_FILE, sep="\t", usecols=["closest_gene", "P"], low_memory=False)
    arab_sig_df = pick_sig(arab, "closest_gene", "P", "Arabidopsis")
    species_genes["Arabidopsis"] = set(arab_sig_df["closest_gene"].unique())
    species_gene_maxp["Arabidopsis"] = arab_sig_df.groupby("closest_gene", sort=False)["P"].max().to_dict()

    rice = pd.read_csv(RICE_FILE, sep="\t", usecols=["closest_gene", "p_wald"], low_memory=False)
    rice_sig_df = pick_sig(rice, "closest_gene", "p_wald", "Rice")
    species_genes["Rice"] = set(rice_sig_df["closest_gene"].unique())
    species_gene_maxp["Rice"] = rice_sig_df.groupby("closest_gene", sort=False)["p_wald"].max().to_dict()

    maize = pd.read_csv(MAIZE_FILE, usecols=["GeneID", "PValue"], low_memory=False)
    maize_sig_df = pick_sig(maize, "GeneID", "PValue", "Maize")
    species_genes["Maize"] = set(maize_sig_df["GeneID"].unique())
    species_gene_maxp["Maize"] = maize_sig_df.groupby("GeneID", sort=False)["PValue"].max().to_dict()

    sorghum = pd.read_csv(SORGHUM_FILE, usecols=["GeneID", "PValue"], low_memory=False)
    sorghum_sig_df = pick_sig(sorghum, "GeneID", "PValue", "Sorghum")
    species_genes["Sorghum"] = set(sorghum_sig_df["GeneID"].unique())
    species_gene_maxp["Sorghum"] = sorghum_sig_df.groupby("GeneID", sort=False)["PValue"].max().to_dict()

    barley = pd.read_csv(BARLEY_FILE, sep="\t", usecols=["closest_gene", "p_wald"], low_memory=False)
    barley_sig_df = pick_sig(barley, "closest_gene", "p_wald", "Barley")
    species_genes["Barley"] = set(barley_sig_df["closest_gene"].unique())
    species_gene_maxp["Barley"] = barley_sig_df.groupby("closest_gene", sort=False)["p_wald"].max().to_dict()

    return species_genes, species_gene_maxp, species_row_stats


def map_to_orthodb_gene_ids(species_genes: dict[str, set[str]]) -> tuple[dict[str, dict[str, set[str]]], dict[str, int]]:
    alias_to_species_genes = defaultdict(list)

    for species, genes in species_genes.items():
        for g in genes:
            for alias in build_aliases(species, g):
                alias_to_species_genes[alias].append((species, g))

    species_gene_to_odb: dict[str, dict[str, set[str]]] = {s: defaultdict(set) for s in species_genes.keys()}

    stats = {
        "genes_lines_scanned": 0,
        "genes_alias_hits": 0,
    }

    with gzip.open(GENES_FILE, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stats["genes_lines_scanned"] += 1
            if stats["genes_lines_scanned"] % 5000000 == 0:
                print(
                    f"[genes] scanned={stats['genes_lines_scanned']:,} alias_hits={stats['genes_alias_hits']:,}",
                    file=sys.stderr,
                    flush=True,
                )

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            tax = parts[1].strip()
            if tax not in ALLOWED_TAXA:
                continue

            odb_gene = parts[0].strip()
            candidates = [p.strip() for p in parts[2:7] if p.strip()]
            desc = parts[7]
            candidates.extend(re.findall(r"\b(?:LOC_Os\d{2}g\d+|Os\d{2}g\d{4,})\b", desc, flags=re.I))

            seen_this_line = set()
            for candidate in candidates:
                cand_u = candidate.strip().upper()
                if cand_u not in alias_to_species_genes:
                    continue
                for sp, original_gene in alias_to_species_genes[cand_u]:
                    key = (sp, original_gene)
                    if key in seen_this_line:
                        continue
                    species_gene_to_odb[sp][original_gene].add(odb_gene)
                    seen_this_line.add(key)
                    stats["genes_alias_hits"] += 1

    return species_gene_to_odb, stats


def map_odb_to_ogs(matched_odb_genes: set[str]) -> dict[str, set[str]]:
    odb_to_ogs: dict[str, set[str]] = defaultdict(set)
    scanned = 0

    with gzip.open(OG2GENES_FILE, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            scanned += 1
            if scanned % 5000000 == 0:
                print(
                    f"[og2genes] scanned={scanned:,} matched_odb_with_og={len(odb_to_ogs):,}",
                    file=sys.stderr,
                    flush=True,
                )
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            og, odb_gene = parts[0], parts[1]
            if odb_gene in matched_odb_genes:
                odb_to_ogs[odb_gene].add(og)

    return odb_to_ogs


def join_to_orthogroups(
    species_gene_to_odb: dict[str, dict[str, set[str]]],
    odb_to_ogs: dict[str, set[str]],
    species_gene_maxp: dict[str, dict[str, float]],
) -> pd.DataFrame:
    og_species_genes: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))

    for species, gene_map in species_gene_to_odb.items():
        for gene, odb_ids in gene_map.items():
            for odb_id in odb_ids:
                for og in odb_to_ogs.get(odb_id, set()):
                    og_species_genes[og][species].add(gene)

    rows = []
    ordered_species = ["Maize", "Sorghum", "Rice", "Barley", "Arabidopsis"]

    for og, sp_map in og_species_genes.items():
        row = {"Orthogroup": og}
        gene_count = 0
        for sp in ordered_species:
            genes = sorted(sp_map.get(sp, set()))
            if genes:
                row[f"{sp}_gene"] = ";".join(genes)
                p_entries = []
                for g in genes:
                    pmax = species_gene_maxp.get(sp, {}).get(g)
                    if pmax is None:
                        p_entries.append(f"{g}:NA")
                    else:
                        p_entries.append(f"{g}:{pmax:.6g}")
                row[f"{sp}_pvalue_max"] = ";".join(p_entries)
                gene_count += 1
            else:
                row[f"{sp}_gene"] = "NA"
                row[f"{sp}_pvalue_max"] = "NA"
        row["Gene_count"] = gene_count
        rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=[
                "Orthogroup",
                "Gene_count",
                "Maize_gene",
                "Maize_pvalue_max",
                "Sorghum_gene",
                "Sorghum_pvalue_max",
                "Rice_gene",
                "Rice_pvalue_max",
                "Barley_gene",
                "Barley_pvalue_max",
                "Arabidopsis_gene",
                "Arabidopsis_pvalue_max",
            ]
        )

    out = pd.DataFrame(rows)
    out = out[
        [
            "Orthogroup",
            "Gene_count",
            "Maize_gene",
            "Maize_pvalue_max",
            "Sorghum_gene",
            "Sorghum_pvalue_max",
            "Rice_gene",
            "Rice_pvalue_max",
            "Barley_gene",
            "Barley_pvalue_max",
            "Arabidopsis_gene",
            "Arabidopsis_pvalue_max",
        ]
    ]
    out = out.sort_values(["Gene_count", "Orthogroup"], ascending=[False, True])
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build cross-species GWAS orthogroup tables.")
    parser.add_argument("--p-threshold", type=float, default=1e-3, help="GWAS p-value threshold (default: 1e-3)")
    parser.add_argument(
        "--top-fraction",
        type=float,
        default=None,
        help="Select top fraction of lowest-p rows per species (e.g., 0.05 for top 5%%). Overrides --p-threshold.",
    )
    parser.add_argument(
        "--prefilter-threshold",
        type=float,
        default=None,
        help="Optional prefilter p-value before --top-fraction selection (e.g., 0.05).",
    )
    parser.add_argument(
        "--use-snp-rank-cutoff",
        action="store_true",
        help="When using --top-fraction, compute species-specific p-cutoffs from total SNP rank and apply those cutoffs to gene tables.",
    )
    args = parser.parse_args()

    global P_THRESHOLD, TOP_FRACTION, PREFILTER_THRESHOLD, USE_SNP_RANK_CUTOFF
    P_THRESHOLD = float(args.p_threshold)
    TOP_FRACTION = args.top_fraction
    PREFILTER_THRESHOLD = args.prefilter_threshold
    USE_SNP_RANK_CUTOFF = args.use_snp_rank_cutoff

    if TOP_FRACTION is not None and not (0 < TOP_FRACTION < 1):
        raise ValueError("--top-fraction must be between 0 and 1.")
    if PREFILTER_THRESHOLD is not None and TOP_FRACTION is None:
        raise ValueError("--prefilter-threshold is only used with --top-fraction.")
    if USE_SNP_RANK_CUTOFF and TOP_FRACTION is None:
        raise ValueError("--use-snp-rank-cutoff requires --top-fraction.")
    if USE_SNP_RANK_CUTOFF and PREFILTER_THRESHOLD is not None:
        raise ValueError("--use-snp-rank-cutoff cannot be combined with --prefilter-threshold.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("[stage] loading significant genes", file=sys.stderr, flush=True)

    species_p_cutoff = None
    snp_stats_df = pd.DataFrame()
    if USE_SNP_RANK_CUTOFF:
        print("[stage] computing SNP-rank cutoffs from total SNPs", file=sys.stderr, flush=True)
        species_p_cutoff, snp_stats_df = compute_species_snp_rank_cutoffs(TOP_FRACTION)

    species_genes, species_gene_maxp, species_row_stats = load_significant_genes(species_p_cutoff=species_p_cutoff)
    print("[stage] mapping gene IDs to OrthoDB gene IDs via genes.tab", file=sys.stderr, flush=True)

    species_gene_to_odb, xref_stats = map_to_orthodb_gene_ids(species_genes)

    matched_odb_genes = set()
    for _, gmap in species_gene_to_odb.items():
        for _, odbs in gmap.items():
            matched_odb_genes.update(odbs)

    print("[stage] mapping OrthoDB gene IDs to orthogroups via OG2genes", file=sys.stderr, flush=True)
    odb_to_ogs = map_odb_to_ogs(matched_odb_genes)

    print("[stage] building final orthogroup tables", file=sys.stderr, flush=True)
    result = join_to_orthogroups(species_gene_to_odb, odb_to_ogs, species_gene_maxp)
    n_species = 5
    commonN = result[result["Gene_count"] == n_species].copy()

    if TOP_FRACTION is None:
        tag = f"p{p_label(P_THRESHOLD)}"
    else:
        if USE_SNP_RANK_CUTOFF:
            tag = f"top{top_label(TOP_FRACTION)}_snpRankCutoff"
        elif PREFILTER_THRESHOLD is None:
            tag = f"top{top_label(TOP_FRACTION)}"
        else:
            tag = f"p{p_label(PREFILTER_THRESHOLD)}_top{top_label(TOP_FRACTION)}"

    out_all = OUT_DIR / f"gwas_orthogroups_{tag}_genesTab_riceFix_all.tsv"
    out_commonN = OUT_DIR / f"gwas_orthogroups_{tag}_genesTab_riceFix_common{n_species}.tsv"
    out_summary = OUT_DIR / f"gwas_orthogroups_{tag}_genesTab_riceFix_summary.txt"
    out_counts = OUT_DIR / f"gwas_orthogroups_{tag}_selection_counts.tsv"
    out_snp_counts = OUT_DIR / f"gwas_orthogroups_{tag}_snp_rank_cutoffs.tsv"

    result.to_csv(out_all, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    commonN.to_csv(out_commonN, sep="\t", index=False, quoting=csv.QUOTE_NONE)

    counts_df = pd.DataFrame(
        [
            {
                "species": sp,
                **species_row_stats.get(sp, {}),
                "unique_genes_selected": len(species_genes.get(sp, set())),
                "p_cutoff_from_snp_rank": (species_p_cutoff.get(sp) if species_p_cutoff else None),
            }
            for sp in ["Arabidopsis", "Rice", "Barley", "Maize", "Sorghum"]
        ]
    )
    counts_df.to_csv(out_counts, sep="\t", index=False)
    if USE_SNP_RANK_CUTOFF:
        snp_stats_df.to_csv(out_snp_counts, sep="\t", index=False)

    selection_mode = "p_threshold"
    if TOP_FRACTION is not None and PREFILTER_THRESHOLD is None:
        selection_mode = "top_fraction"
    if TOP_FRACTION is not None and PREFILTER_THRESHOLD is not None:
        selection_mode = "p_threshold_then_top_fraction"
    if USE_SNP_RANK_CUTOFF:
        selection_mode = "snp_rank_cutoff_then_gene_filter"

    summary_lines = [
        f"selection_mode={selection_mode}",
        f"p_threshold={P_THRESHOLD}",
        f"top_fraction={TOP_FRACTION}",
        f"prefilter_threshold={PREFILTER_THRESHOLD}",
        f"use_snp_rank_cutoff={USE_SNP_RANK_CUTOFF}",
        f"arabidopsis_unique_sig_genes={len(species_genes['Arabidopsis'])}",
        f"rice_unique_sig_genes={len(species_genes['Rice'])}",
        f"barley_unique_sig_genes={len(species_genes['Barley'])}",
        f"maize_unique_sig_genes={len(species_genes['Maize'])}",
        f"sorghum_unique_sig_genes={len(species_genes['Sorghum'])}",
        f"genes_lines_scanned={xref_stats['genes_lines_scanned']}",
        f"genes_alias_hits={xref_stats['genes_alias_hits']}",
        f"mapped_arabidopsis_genes={len(species_gene_to_odb['Arabidopsis'])}",
        f"mapped_rice_genes={len(species_gene_to_odb['Rice'])}",
        f"mapped_barley_genes={len(species_gene_to_odb['Barley'])}",
        f"mapped_maize_genes={len(species_gene_to_odb['Maize'])}",
        f"mapped_sorghum_genes={len(species_gene_to_odb['Sorghum'])}",
        f"matched_orthodb_genes={len(matched_odb_genes)}",
        f"orthogroups_total={len(result)}",
        f"orthogroups_common{n_species}={len(commonN)}",
        f"output_all={out_all}",
        f"output_common{n_species}={out_commonN}",
        f"output_counts={out_counts}",
        f"output_snp_rank_cutoffs={out_snp_counts if USE_SNP_RANK_CUTOFF else 'NA'}",
    ]
    out_summary.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"Wrote: {out_all}")
    print(f"Wrote: {out_commonN}")
    print(f"Wrote: {out_summary}")


if __name__ == "__main__":
    main()
