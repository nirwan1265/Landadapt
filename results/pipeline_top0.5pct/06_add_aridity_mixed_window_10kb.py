#!/usr/bin/env python3

import csv
import gzip
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path("/Users/nirwantandukar/Documents/Github/Landadapt/results")
OUT = ROOT / "mixed_window_10kb_csv"
ANNOT_DIR = OUT / "annotated_10kb_csv"
ORTHO_DIR = OUT / "orthogroups_csv"
TEX_FILE = OUT / "cross_species_top0p5_mixed_window_10kb_sections.tex"
SUMMARY_FILE = OUT / "orthogroups_summary_mixed_window_10kb.csv"

for d in [OUT, ANNOT_DIR, ORTHO_DIR]:
    d.mkdir(parents=True, exist_ok=True)

ANNOT_SCRIPT = Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/pipeline_top0.5pct/01_gene_annotate_generic.R")
GENES_FILE = Path("/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_genes.tab.gz")
OG2GENES_FILE = Path("/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_OG2genes.tab.gz")

TRAIT = "aridity_index"
LABEL = "aridity index"
TOP_FRAC = 0.005
SPECIES_ORDER = ["Maize", "Sorghum", "Rice", "Barley", "Arabidopsis"]
ALLOWED_TAXA = {"3702_0", "4577_0", "4558_0", "4557_0", "39947_0", "4530_0", "4513_0", "112509_0"}

GWAS = {
    "Arabidopsis": Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index/AT_aridity_index.txt"),
    "Rice": Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index/rice_aridity_index.txt"),
    "Maize": Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index/maize_aridity_index.txt"),
    "Sorghum": Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index/sorghum_aridity_index.txt"),
    "Barley": Path("/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index/barley_aridity_index.txt"),
}

GFF = {
    "Arabidopsis": Path("/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff"),
    "Rice": Path("/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3"),
    "Maize": Path("/Users/nirwantandukar/Documents/Research/data/maize_gene_annotation/ENSEMBLE_Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.chr.gff3"),
    "Sorghum": Path("/Users/nirwantandukar/Documents/Research/data/sorghum_annotation/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3"),
    "Barley": Path("/Users/nirwantandukar/Documents/Research/results/Barley/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3"),
}

WINDOW = {
    "Arabidopsis": 10000,
    "Rice": 10000,
    "Maize": 25000,
    "Sorghum": 25000,
    "Barley": 25000,
}


def resolve_col(columns, candidates):
    low = {c.lower(): c for c in columns}
    for cand in candidates:
        hit = low.get(cand.lower())
        if hit:
            return hit
    return None


def load_gwas(path: Path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    snp_col = resolve_col(df.columns, ["rs", "SNP", "snp"])
    chr_col = resolve_col(df.columns, ["chr", "Chr", "CHR", "chrom", "CHROM"])
    pos_col = resolve_col(df.columns, ["ps", "pos", "bp", "BP", "POS"])
    p_col = resolve_col(df.columns, ["p_wald", "P", "pvalue", "P.value"])
    if None in (snp_col, chr_col, pos_col, p_col):
        raise ValueError(f"Could not resolve GWAS columns in {path}")

    out = pd.DataFrame({
        "snp": df[snp_col].astype(str),
        "chr": df[chr_col].astype(str),
        "pos": pd.to_numeric(df[pos_col], errors="coerce"),
        "p": pd.to_numeric(df[p_col], errors="coerce"),
    })
    out = out.dropna(subset=["chr", "pos", "p"])
    out = out[(out["chr"] != "") & (out["p"] > 0)]
    out = out.groupby(["snp", "chr", "pos"], as_index=False)["p"].min()
    return out


def select_top(df):
    n_total = len(df)
    if n_total == 0:
        return df.copy(), 0, 0, float("nan"), float("nan")
    n_rank = max(1, math.floor(n_total * TOP_FRAC))
    cutoff = df["p"].sort_values().iloc[n_rank - 1]
    bonf_cutoff = 0.05 / n_total
    top = df[df["p"] <= cutoff].copy()
    n_bonf = int((df["p"] <= bonf_cutoff).sum())
    return top, n_total, n_bonf, float(cutoff), float(bonf_cutoff)


def annotate_species(species, top_df):
    tmp_in = ANNOT_DIR / f"{TRAIT}_{species}_selected_for_annot.tsv"
    tmp_out = ANNOT_DIR / f"{TRAIT}_{species}_annot_tmp.tsv"
    out_csv = ANNOT_DIR / f"{TRAIT}_{species}_annot_{'10kb' if WINDOW[species] == 10000 else '25kb'}.csv"

    w = pd.DataFrame({
        "chr": top_df["chr"].astype(str),
        "rs": top_df["snp"].astype(str),
        "ps": pd.to_numeric(top_df["pos"], errors="coerce"),
        "p_wald": pd.to_numeric(top_df["p"], errors="coerce"),
    }).dropna(subset=["ps", "p_wald"])
    w.to_csv(tmp_in, sep="\t", index=False)

    cmd = ["Rscript", str(ANNOT_SCRIPT), str(tmp_in), str(GFF[species]), str(tmp_out), str(WINDOW[species])]
    subprocess.run(cmd, check=True)

    ann = pd.read_csv(tmp_out, sep="\t", low_memory=False)
    ann.to_csv(out_csv, index=False)

    if tmp_in.exists():
        tmp_in.unlink()
    if tmp_out.exists():
        tmp_out.unlink()

    if "closest_gene" not in ann.columns or "p_wald" not in ann.columns:
        return {}

    ann["closest_gene"] = ann["closest_gene"].astype(str).str.strip()
    ann["p_wald"] = pd.to_numeric(ann["p_wald"], errors="coerce")
    ann = ann[(ann["closest_gene"] != "") & (~ann["closest_gene"].str.lower().isin(["na", "nan", "none"]))]
    ann = ann.dropna(subset=["p_wald"])

    gm = ann.groupby("closest_gene", as_index=False)["p_wald"].max()
    return {r["closest_gene"]: float(r["p_wald"]) for _, r in gm.iterrows()}


def normalize_gene_id(gene_id: str) -> str:
    g = str(gene_id).strip().upper()
    g = re.sub(r"^(GENE:|TRANSCRIPT:)", "", g)
    g = re.sub(r"_T\d+$", "", g)
    g = re.sub(r"_P\d+$", "", g)
    g = re.sub(r"\.[0-9]+$", "", g)
    g = re.sub(r"-[A-Z]{1,3}$", "", g)
    return g


def build_aliases(species: str, gene_id: str):
    aliases = {str(gene_id).strip().upper(), normalize_gene_id(gene_id)}
    if species == "Rice":
        m = re.match(r"^LOC_OS(\d{2})G(\d+)$", normalize_gene_id(gene_id))
        if m:
            chrom, num = m.groups()
            aliases.update({f"OS{chrom}G{num}", f"OS{chrom}G{num}0", f"OS{chrom}G{num}00", f"OS{chrom}G{num.zfill(7)}"})
    return aliases


def map_to_orthodb_gene_ids(species_genes):
    alias_to_species_genes = defaultdict(list)
    for species, genes in species_genes.items():
        for g in genes:
            for alias in build_aliases(species, g):
                alias_to_species_genes[alias].append((species, g))

    species_gene_to_odb = {s: defaultdict(set) for s in species_genes.keys()}
    with gzip.open(GENES_FILE, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
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

            seen = set()
            for candidate in candidates:
                cu = candidate.upper()
                if cu not in alias_to_species_genes:
                    continue
                for sp, original in alias_to_species_genes[cu]:
                    key = (sp, original)
                    if key in seen:
                        continue
                    species_gene_to_odb[sp][original].add(odb_gene)
                    seen.add(key)
    return species_gene_to_odb


def map_odb_to_ogs(matched_odb_genes):
    odb_to_ogs = defaultdict(set)
    with gzip.open(OG2GENES_FILE, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            og, odb_gene = parts[0], parts[1]
            if odb_gene in matched_odb_genes:
                odb_to_ogs[odb_gene].add(og)
    return odb_to_ogs


def join_to_orthogroups(species_gene_to_odb, odb_to_ogs, species_gene_pmax):
    og_species_genes = defaultdict(lambda: defaultdict(set))
    for species, gene_map in species_gene_to_odb.items():
        for gene, odb_ids in gene_map.items():
            for odb_id in odb_ids:
                for og in odb_to_ogs.get(odb_id, set()):
                    og_species_genes[og][species].add(gene)

    rows = []
    for og, sp_map in og_species_genes.items():
        row = {"Orthogroup": og}
        count = 0
        for sp in SPECIES_ORDER:
            genes = sorted(sp_map.get(sp, set()))
            if genes:
                row[f"{sp}_gene"] = ";".join(genes)
                p_entries = []
                for g in genes:
                    pmax = species_gene_pmax.get(sp, {}).get(g)
                    p_entries.append(f"{g}:{pmax:.6g}" if pmax is not None else f"{g}:NA")
                row[f"{sp}_pvalue_max"] = ";".join(p_entries)
                count += 1
            else:
                row[f"{sp}_gene"] = "NA"
                row[f"{sp}_pvalue_max"] = "NA"
        row["Gene_count"] = count
        rows.append(row)

    cols = ["Orthogroup", "Gene_count"] + [x for sp in SPECIES_ORDER for x in (f"{sp}_gene", f"{sp}_pvalue_max")]
    if not rows:
        return pd.DataFrame(columns=cols)
    out = pd.DataFrame(rows, columns=cols).sort_values(["Gene_count", "Orthogroup"], ascending=[False, True]).reset_index(drop=True)
    return out


def tex_escape(s):
    return str(s).replace("_", r"\_")


def pvals_to_logp_list(cell):
    t = str(cell).strip()
    if t == "" or t.upper() == "NA":
        return "NA"
    vals = []
    for part in t.split(";"):
        if ":" not in part:
            vals.append("NA")
            continue
        _, p = part.split(":", 1)
        try:
            pv = float(p)
            vals.append(f"{-math.log10(pv):.2f}" if pv > 0 else "NA")
        except Exception:
            vals.append("NA")
    return ";".join(vals) if vals else "NA"


def build_tex_subsection(info, common5_df):
    lines = []
    lines.append(rf"\subsection{{Cross-species genome-wide associations with {LABEL} identify shared orthologous genes across all species}}")
    lines.append("")
    lines.append("Aridity index represents long-term water availability relative to atmospheric evaporative demand, providing an integrated climatic dryness gradient across sampling locations.")
    lines.append("")
    lines.append(
        "We performed GWAS in Arabidopsis, barley, rice, maize, and sorghum using a mixed linear model and kinship correction in GEMMA. "
        f"Under Bonferroni correction, the number of significant SNPs was heterogeneous across species (Arabidopsis {info['bonf'].get('Arabidopsis',0)}, rice {info['bonf'].get('Rice',0)}, sorghum {info['bonf'].get('Sorghum',0)}, barley {info['bonf'].get('Barley',0)}, maize {info['bonf'].get('Maize',0)}). "
        "To evaluate cross-species convergence beyond the Bonferroni tail, we selected the top 0.5\\% of SNPs per species and mapped nearby genes to orthogroups, using \\pm10 kb annotation for Arabidopsis and rice and \\pm25 kb annotation for maize, sorghum, and barley. "
        f"This yielded {info['genes'].get('Arabidopsis',0):,} Arabidopsis genes, {info['genes'].get('Rice',0):,} rice genes, {info['genes'].get('Maize',0):,} maize genes, {info['genes'].get('Barley',0):,} barley genes, and {info['genes'].get('Sorghum',0):,} sorghum genes."
    )
    lines.append("")
    lines.append(
        f"Orthogroup integration identified {info['og_total']:,} total orthogroups, with {info['og_common5']:,} orthogroups represented in all five species. "
        f"Table~\\ref{{tab:shared_function_candidates_{TRAIT}_mixed10kb}} summarizes these shared orthogroups and reports, for each species, the full gene list within each orthogroup and corresponding $-\\log_{{10}}(P)$ values."
    )
    lines.append("")
    lines.append(r"\begin{table*}[t]")
    lines.append(rf"\caption{{Shared orthogroups across all five species from the top 0.5\% of SNPs for {LABEL} (Arabidopsis/Rice \pm10 kb; Maize/Sorghum/Barley \pm25 kb).\label{{tab:shared_function_candidates_{TRAIT}_mixed10kb}}}}")
    lines.append(r"\scriptsize")
    lines.append(r"\tabcolsep=3pt")
    lines.append(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllllllllll@{}}")
    lines.append(r"\toprule")
    lines.append(r"Orthogroup & Maize genes & $-\log_{10}(P)$ & Sorghum genes & $-\log_{10}(P)$ & Rice genes & $-\log_{10}(P)$ & Barley genes & $-\log_{10}(P)$ & Arabidopsis genes & $-\log_{10}(P)$ \\")
    lines.append(r"\midrule")

    if common5_df.empty:
        lines.append(r"NA & NA & NA & NA & NA & NA & NA & NA & NA & NA & NA \\")
    else:
        for _, r in common5_df.sort_values("Orthogroup").iterrows():
            row = [tex_escape(r["Orthogroup"])]
            for sp in SPECIES_ORDER:
                genes = str(r.get(f"{sp}_gene", "NA"))
                row.append(tex_escape(genes) if genes.strip() else "NA")
                row.append(pvals_to_logp_list(r.get(f"{sp}_pvalue_max", "NA")))
            lines.append(" & ".join(row) + r" \\")

    lines.append(r"\botrule")
    lines.append(r"\end{tabular*}")
    lines.append(r"\begin{tablenotes}%")
    lines.append(r"\item Note: All genes detected per species within each orthogroup are shown (semicolon-separated), alongside semicolon-separated $-\log_{10}(P)$ values in matching order.")
    lines.append(r"\end{tablenotes}")
    lines.append(r"\end{table*}")
    lines.append("")
    return "\n".join(lines)


def main():
    species_gene_pmax = {sp: {} for sp in SPECIES_ORDER}
    top_counts = {}
    bonf_counts = {}
    gene_counts = {}

    all_gene_rows = []
    count_rows = []

    for sp in SPECIES_ORDER:
        gdf = load_gwas(GWAS[sp])
        top_df, n_total, n_bonf, p_cutoff, bonf_cutoff = select_top(gdf)
        gene_map = annotate_species(sp, top_df)

        species_gene_pmax[sp] = gene_map
        top_counts[sp] = len(top_df)
        bonf_counts[sp] = n_bonf
        gene_counts[sp] = len(gene_map)

        for g, p in gene_map.items():
            all_gene_rows.append({"species": sp, "gene": g, "pvalue": p})

        count_rows.append({
            "trait": TRAIT,
            "species": sp,
            "top0p5_snps": len(top_df),
            "bonf_snps": n_bonf,
            "genes_mixed_window_10kb": len(gene_map),
            "annotation_window_bp": WINDOW[sp],
            "n_total_snps": n_total,
            "top0p5_p_cutoff": p_cutoff,
            "bonf_p_cutoff": bonf_cutoff,
        })

    gene_tbl = pd.DataFrame(all_gene_rows)
    if not gene_tbl.empty:
        gene_tbl = gene_tbl.sort_values(["species", "pvalue", "gene"])
    gene_tbl.to_csv(OUT / f"{TRAIT}_top0p5_gene_table_mixed_window_10kb.csv", index=False)
    pd.DataFrame(count_rows).to_csv(OUT / f"{TRAIT}_top0p5_counts_mixed_window_10kb.csv", index=False)

    species_genes = {sp: set(species_gene_pmax[sp].keys()) for sp in SPECIES_ORDER}
    species_gene_to_odb = map_to_orthodb_gene_ids(species_genes)
    matched_odb = set()
    for gm in species_gene_to_odb.values():
        for s in gm.values():
            matched_odb.update(s)
    odb_to_ogs = map_odb_to_ogs(matched_odb)

    all_df = join_to_orthogroups(species_gene_to_odb, odb_to_ogs, species_gene_pmax)
    ge2_df = all_df[all_df["Gene_count"] >= 2].copy() if not all_df.empty else all_df.copy()
    c5_df = all_df[all_df["Gene_count"] == 5].copy() if not all_df.empty else all_df.copy()

    all_df.to_csv(ORTHO_DIR / f"{TRAIT}_top0p5_orthogroups_all_mixed_window_10kb.csv", index=False)
    ge2_df.to_csv(ORTHO_DIR / f"{TRAIT}_top0p5_orthogroups_common_ge2_mixed_window_10kb.csv", index=False)
    c5_df.to_csv(ORTHO_DIR / f"{TRAIT}_top0p5_orthogroups_common5_mixed_window_10kb.csv", index=False)

    if SUMMARY_FILE.exists():
        s = pd.read_csv(SUMMARY_FILE)
    else:
        s = pd.DataFrame(columns=["trait", "orthogroups_total", "orthogroups_common_ge2", "orthogroups_common5"])
    s = s[s["trait"] != TRAIT]
    s = pd.concat(
        [s, pd.DataFrame([{
            "trait": TRAIT,
            "orthogroups_total": len(all_df),
            "orthogroups_common_ge2": len(ge2_df),
            "orthogroups_common5": len(c5_df),
        }])],
        ignore_index=True
    )
    s.to_csv(SUMMARY_FILE, index=False)

    info = {
        "bonf": bonf_counts,
        "genes": gene_counts,
        "og_total": len(all_df),
        "og_common5": len(c5_df),
    }
    new_sec = build_tex_subsection(info, c5_df)

    existing = TEX_FILE.read_text(encoding="utf-8") if TEX_FILE.exists() else ""
    anchor = r"\subsection{Cross-species genome-wide associations with aridity index identify shared orthologous genes across all species}"
    if anchor in existing:
        start = existing.find(anchor)
        nxt = existing.find(r"\subsection{", start + 1)
        if nxt == -1:
            updated = existing[:start].rstrip() + "\n\n" + new_sec + "\n"
        else:
            updated = existing[:start].rstrip() + "\n\n" + new_sec + "\n\n" + existing[nxt:].lstrip()
    else:
        updated = existing.rstrip() + "\n\n" + new_sec + "\n"
    TEX_FILE.write_text(updated, encoding="utf-8")

    print("DONE")
    print(f"Gene table: {OUT / (TRAIT + '_top0p5_gene_table_mixed_window_10kb.csv')}")
    print(f"Counts: {OUT / (TRAIT + '_top0p5_counts_mixed_window_10kb.csv')}")
    print(f"Orthogroups all: {ORTHO_DIR / (TRAIT + '_top0p5_orthogroups_all_mixed_window_10kb.csv')}")
    print(f"TeX updated: {TEX_FILE}")


if __name__ == "__main__":
    main()
