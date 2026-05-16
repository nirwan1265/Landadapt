from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ALL_DIR = BASE / "results" / "mixed_window_10kb_csv" / "orthogroups_csv" / "All"
PC_DEDUP = BASE / "results" / "tables" / "Supplementary" / "PC" / "SuppTable3A_PC_orthogroups_atleast_2_desc.csv"
NOPC_DEDUP = BASE / "results" / "tables" / "Supplementary" / "SuppTable3A_orthogroups_atleast_2_desc_noPC.csv"
OUTDIR = BASE / "results" / "tables" / "Null_model"
OUTDIR.mkdir(parents=True, exist_ok=True)

TRAITS = [
    "PC1",
    "PC2",
    "PC3",
    "pH",
    "soilN",
    "am_rel_abundance_colonization",
    "am_roots_colonized",
    "aridity_index",
]
ALL_FILE_TRAIT = {
    "PC1": "pc1_worldclim",
    "PC2": "pc2_worldclim",
    "PC3": "pc3_worldclim",
    "pH": "ph",
    "soilN": "soilN",
    "am_rel_abundance_colonization": "am_rel_abundance_colonization",
    "am_roots_colonized": "am_roots_colonized",
    "aridity_index": "aridity_index",
}
SPECIES_COLS = {
    "Arabidopsis": "Arabidopsis_gene",
    "Barley": "Barley_gene",
    "Rice": "Rice_gene",
    "Maize": "Maize_gene",
    "Sorghum": "Sorghum_gene",
}
PHENO_LABEL = {
    "PC1": "PC1",
    "PC2": "PC2",
    "PC3": "PC3",
    "pH": "pH",
    "soilN": "soilN",
    "am_rel_abundance_colonization": "AM_rel",
    "am_roots_colonized": "AM_roots",
    "aridity_index": "aridity",
}


def split_genes(cell) -> list[str]:
    if pd.isna(cell):
        return []
    genes = []
    for gene in str(cell).split(";"):
        gene = gene.strip()
        if gene and gene.lower() != "nan":
            genes.append(gene)
    return genes


def trait_occ_bin(n: int) -> int:
    return min(int(n), 4)


def og_mult_bin(n: int) -> int:
    if n <= 1:
        return 1
    if n == 2:
        return 2
    if n == 3:
        return 3
    return 4


def load_background_metadata() -> pd.DataFrame:
    trait_presence = defaultdict(set)
    og_presence = defaultdict(set)
    for trait_file in ALL_FILE_TRAIT.values():
        df = pd.read_csv(ALL_DIR / f"{trait_file}_top0p5_orthogroups_all_mixed_window_10kb.csv")
        for _, row in df.iterrows():
            og = row["Orthogroup"]
            for species, col in SPECIES_COLS.items():
                for gene in split_genes(row[col]):
                    key = (species, gene)
                    trait_presence[key].add(trait_file)
                    og_presence[key].add(og)
    records = []
    for (species, gene), traits in trait_presence.items():
        ogs = og_presence[(species, gene)]
        records.append(
            {
                "species": species,
                "gene": gene,
                "trait_occurrence": len(traits),
                "trait_bin": trait_occ_bin(len(traits)),
                "global_og_count": len(ogs),
                "og_bin": og_mult_bin(len(ogs)),
            }
        )
    return pd.DataFrame(records)


def load_observed_gene_sets() -> dict[str, dict[str, set[str]]]:
    observed = {}
    for trait, trait_file in ALL_FILE_TRAIT.items():
        df = pd.read_csv(ALL_DIR / f"{trait_file}_top0p5_orthogroups_all_mixed_window_10kb.csv")
        observed[trait] = {}
        for species, col in SPECIES_COLS.items():
            genes = set()
            for cell in df[col].dropna().astype(str):
                genes.update(split_genes(cell))
            observed[trait][species] = genes
    return observed


def load_dedup_rows() -> tuple[dict[str, pd.DataFrame], dict[str, dict[str, int]]]:
    pc = pd.read_csv(PC_DEDUP)
    no = pd.read_csv(NOPC_DEDUP)
    no = no[~no["Phenotype"].isin(["PC1", "PC2", "PC3"])].copy()
    all_rows = pd.concat([pc, no], ignore_index=True)

    rows = {}
    observed_counts = {}
    for trait in TRAITS:
        sub = all_rows[all_rows["Phenotype"] == trait].copy().reset_index(drop=True)
        rows[trait] = sub
        observed_counts[trait] = {
            "all5": int((sub["Gene_count"] == 5).sum()),
            "ge4": int((sub["Gene_count"] >= 4).sum()),
            "n_rows": len(sub),
        }
    return rows, observed_counts


def build_species_structures(meta: pd.DataFrame, observed_genes):
    species_data = {}
    for species, sub in meta.groupby("species"):
        sub = sub.copy().reset_index(drop=True)
        sub["gene_index"] = np.arange(len(sub), dtype=int)
        gene_to_idx = dict(zip(sub["gene"], sub["gene_index"]))

        pair_pools = {}
        trait_pools = {}
        for pair, pool in sub.groupby(["trait_bin", "og_bin"]):
            pair_pools[pair] = pool["gene_index"].to_numpy(dtype=int)
        for trait_bin, pool in sub.groupby("trait_bin"):
            trait_pools[int(trait_bin)] = pool["gene_index"].to_numpy(dtype=int)

        observed_by_trait = {}
        for trait in TRAITS:
            idx = np.array([gene_to_idx[g] for g in observed_genes[trait][species] if g in gene_to_idx], dtype=int)
            obs_sub = sub.iloc[idx]
            need_by_pair = Counter(zip(obs_sub["trait_bin"], obs_sub["og_bin"]))
            observed_by_trait[trait] = {
                "n": len(idx),
                "need_by_pair": need_by_pair,
            }

        species_data[species] = {
            "meta": sub,
            "gene_to_idx": gene_to_idx,
            "all_idx": sub["gene_index"].to_numpy(dtype=int),
            "pair_pools": pair_pools,
            "trait_pools": trait_pools,
            "observed_by_trait": observed_by_trait,
        }
    return species_data


def build_trait_row_indices(rows_by_trait: dict[str, pd.DataFrame], species_data):
    trait_row_idx = {}
    for trait, df in rows_by_trait.items():
        species_rows = {}
        for species, col in SPECIES_COLS.items():
            gene_to_idx = species_data[species]["gene_to_idx"]
            row_arrays = []
            for cell in df[col].tolist():
                idx = [gene_to_idx[g] for g in split_genes(cell) if g in gene_to_idx]
                row_arrays.append(np.array(idx, dtype=int))
            species_rows[species] = row_arrays
        trait_row_idx[trait] = species_rows
    return trait_row_idx


def sample_matched(species_info, trait: str, rng) -> np.ndarray:
    obs_info = species_info["observed_by_trait"][trait]
    selected_mask = np.zeros(len(species_info["all_idx"]), dtype=bool)
    remaining_by_trait = Counter()

    for pair, need in sorted(obs_info["need_by_pair"].items()):
        pool = species_info["pair_pools"].get(pair)
        if pool is None:
            remaining_by_trait[pair[0]] += need
            continue
        available = pool[~selected_mask[pool]]
        take = min(need, len(available))
        if take > 0:
            chosen = rng.choice(available, size=take, replace=False)
            selected_mask[chosen] = True
        if take < need:
            remaining_by_trait[pair[0]] += (need - take)

    for trait_bin, need in sorted(remaining_by_trait.items()):
        pool = species_info["trait_pools"].get(trait_bin)
        if pool is None:
            continue
        available = pool[~selected_mask[pool]]
        take = min(need, len(available))
        if take > 0:
            chosen = rng.choice(available, size=take, replace=False)
            selected_mask[chosen] = True

    still_needed = obs_info["n"] - int(selected_mask.sum())
    if still_needed > 0:
        pool = species_info["all_idx"]
        available = pool[~selected_mask[pool]]
        chosen = rng.choice(available, size=still_needed, replace=False)
        selected_mask[chosen] = True

    return selected_mask


def count_hits(trait: str, sampled_masks, trait_row_idx) -> tuple[int, int]:
    n_rows = len(next(iter(trait_row_idx[trait].values())))
    hits = np.zeros(n_rows, dtype=np.int16)
    for species in SPECIES_COLS:
        mask = sampled_masks[species]
        for i, idx in enumerate(trait_row_idx[trait][species]):
            if idx.size and mask[idx].any():
                hits[i] += 1
    ge4 = int((hits >= 4).sum())
    all5 = int((hits == 5).sum())
    return all5, ge4


def run_permutations(n_perm: int, seed: int):
    meta = load_background_metadata()
    observed_genes = load_observed_gene_sets()
    rows_by_trait, observed_counts = load_dedup_rows()
    species_data = build_species_structures(meta, observed_genes)
    trait_row_idx = build_trait_row_indices(rows_by_trait, species_data)
    rng = np.random.default_rng(seed)

    dist_rows = []
    summary_rows = []

    for trait in TRAITS:
        perm_all5 = np.zeros(n_perm, dtype=int)
        perm_ge4 = np.zeros(n_perm, dtype=int)
        for i in range(n_perm):
            sampled_masks = {}
            for species in SPECIES_COLS:
                sampled_masks[species] = sample_matched(species_data[species], trait, rng)
            a5, g4 = count_hits(trait, sampled_masks, trait_row_idx)
            perm_all5[i] = a5
            perm_ge4[i] = g4

        obs_a5 = observed_counts[trait]["all5"]
        obs_g4 = observed_counts[trait]["ge4"]
        p_a5 = (1 + int((perm_all5 >= obs_a5).sum())) / (n_perm + 1)
        p_g4 = (1 + int((perm_ge4 >= obs_g4).sum())) / (n_perm + 1)

        summary_rows.append(
            {
                "Phenotype": trait,
                "Trait_label": PHENO_LABEL[trait],
                "Observed_all5": obs_a5,
                "Null_mean_all5": perm_all5.mean(),
                "Null_sd_all5": perm_all5.std(ddof=1),
                "Null_q025_all5": np.quantile(perm_all5, 0.025),
                "Null_q975_all5": np.quantile(perm_all5, 0.975),
                "Empirical_p_all5": p_a5,
                "Observed_ge4": obs_g4,
                "Null_mean_ge4": perm_ge4.mean(),
                "Null_sd_ge4": perm_ge4.std(ddof=1),
                "Null_q025_ge4": np.quantile(perm_ge4, 0.025),
                "Null_q975_ge4": np.quantile(perm_ge4, 0.975),
                "Empirical_p_ge4": p_g4,
                "Observed_gene_universe_Arabidopsis": species_data["Arabidopsis"]["observed_by_trait"][trait]["n"],
                "Observed_gene_universe_Barley": species_data["Barley"]["observed_by_trait"][trait]["n"],
                "Observed_gene_universe_Rice": species_data["Rice"]["observed_by_trait"][trait]["n"],
                "Observed_gene_universe_Maize": species_data["Maize"]["observed_by_trait"][trait]["n"],
                "Observed_gene_universe_Sorghum": species_data["Sorghum"]["observed_by_trait"][trait]["n"],
                "Background_Arabidopsis": len(species_data["Arabidopsis"]["all_idx"]),
                "Background_Barley": len(species_data["Barley"]["all_idx"]),
                "Background_Rice": len(species_data["Rice"]["all_idx"]),
                "Background_Maize": len(species_data["Maize"]["all_idx"]),
                "Background_Sorghum": len(species_data["Sorghum"]["all_idx"]),
                "Permutations": n_perm,
                "Seed": seed,
            }
        )

        for j in range(n_perm):
            dist_rows.append(
                {
                    "Phenotype": trait,
                    "Trait_label": PHENO_LABEL[trait],
                    "perm_index": j + 1,
                    "all5_count": int(perm_all5[j]),
                    "ge4_count": int(perm_ge4[j]),
                }
            )

    return pd.DataFrame(summary_rows), pd.DataFrame(dist_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260513)
    args = parser.parse_args()

    summary, dist = run_permutations(args.n_perm, args.seed)
    summary_path = OUTDIR / "orthogroup_null_permutation_summary.csv"
    dist_path = OUTDIR / "orthogroup_null_permutation_distributions.csv"
    summary.to_csv(summary_path, index=False)
    dist.to_csv(dist_path, index=False)
    print("WROTE", summary_path)
    print("WROTE", dist_path)
    print(
        summary[
            [
                "Trait_label",
                "Observed_all5",
                "Null_mean_all5",
                "Empirical_p_all5",
                "Observed_ge4",
                "Null_mean_ge4",
                "Empirical_p_ge4",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
