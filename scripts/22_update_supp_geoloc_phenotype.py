from pathlib import Path

import pandas as pd


ROOT = Path(".")
WORLDCLIM = ROOT / "data/WorldClim/worldclim_all_crops_BIO_PCA.csv"
ARIDITY = ROOT / "data/Aridity/all_crops_aridity_values.csv"
PH = ROOT / "data/pH/all_crops_pH_values.csv"
SOILN = ROOT / "data/TN/all_crops_N_values.csv"
AM_REL = ROOT / "data/AM_rel_abundance_colonization/all_crops_AM_rel_abundance_colonization_values.csv"
AM_ROOTS = ROOT / "data/AM_roots_colonized/all_crops_AM_roots_colonized_values.csv"
OUT = ROOT / "results/tables/Supplementary/SuppTable1_geoloc_phenotype.csv"


species_map = {
    "arabidopsis": "Arabidopsis",
    "barley": "Barley",
    "maize": "Maize",
    "rice": "Rice",
    "sorghum": "Sorghum",
}


def main() -> None:
    base = pd.read_csv(WORLDCLIM)
    aridity = pd.read_csv(ARIDITY)[["Lines", "crop", "aridity_index"]]
    ph = pd.read_csv(PH)[["Lines", "crop", "ph_value"]]
    soiln = pd.read_csv(SOILN)[["Lines", "crop", "n_value"]]
    am_rel = pd.read_csv(AM_REL)[["Lines", "crop", "am_rel_abundance_colonization_value"]]
    am_roots = pd.read_csv(AM_ROOTS)[["Lines", "crop", "am_roots_colonized_value"]]

    merged = (
        base.merge(aridity, on=["Lines", "crop"], how="left")
        .merge(ph, on=["Lines", "crop"], how="left")
        .merge(soiln, on=["Lines", "crop"], how="left")
        .merge(am_rel, on=["Lines", "crop"], how="left")
        .merge(am_roots, on=["Lines", "crop"], how="left")
    )

    merged["Species"] = merged["crop"].map(species_map)

    ordered_cols = (
        ["Lines", "Long", "Lat", "Species"]
        + [f"BIO{i:02d}" for i in range(1, 20)]
        + [
            "PC1",
            "PC2",
            "PC3",
            "aridity_index",
            "ph_value",
            "n_value",
            "am_rel_abundance_colonization_value",
            "am_roots_colonized_value",
        ]
    )
    out = merged[ordered_cols].rename(
        columns={
            "Lines": "Taxa",
            "ph_value": "pH",
            "n_value": "soilN",
            "am_rel_abundance_colonization_value": "AM_rel_abundance_colonization",
            "am_roots_colonized_value": "AM_roots_colonized",
        }
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"Wrote: {OUT}")
    print(f"Rows: {len(out)}")


if __name__ == "__main__":
    main()
