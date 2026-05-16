import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary'
OUTDIR.mkdir(parents=True, exist_ok=True)

TRAITS = {
    'PC1': BASE / 'results' / 'pc1_worldclim' / 'pc1_worldclim_top0p5_selected_snps.csv',
    'PC2': BASE / 'results' / 'pc2_worldclim' / 'pc2_worldclim_top0p5_selected_snps.csv',
    'PC3': BASE / 'results' / 'pc3_worldclim' / 'pc3_worldclim_top0p5_selected_snps.csv',
    'pH': BASE / 'results' / 'ph' / 'ph_top0p5_selected_snps.csv',
    'soilN': BASE / 'results' / 'soilN' / 'soilN_top0p5_selected_snps.csv',
    'aridity_index': BASE / 'results' / 'aridity_index' / 'aridity_index_top0p5_selected_snps.csv',
    'AM fungal relative abundance colonization': BASE / 'results' / 'am_rel_abundance_colonization' / 'am_rel_abundance_colonization_top0p5_selected_snps.csv',
    'AM fungal roots colonized': BASE / 'results' / 'am_roots_colonized' / 'am_roots_colonized_top0p5_selected_snps.csv',
}


def load_with_phenotype(phenotype, path):
    df = pd.read_csv(path)
    df.insert(0, 'Phenotype', phenotype)
    return df


def main():
    pc = pd.concat([load_with_phenotype(k, TRAITS[k]) for k in ['PC1', 'PC2', 'PC3']], ignore_index=True)
    pc.to_csv(OUTDIR / 'SuppTable5_WorldClim_PC1_PC2_PC3_top0p5_GWAS_SNPs.csv', index=False)

    ph = load_with_phenotype('pH', TRAITS['pH'])
    ph.to_csv(OUTDIR / 'SuppTable6_pH_top0p5_GWAS_SNPs.csv', index=False)

    soiln = load_with_phenotype('soilN', TRAITS['soilN'])
    soiln.to_csv(OUTDIR / 'SuppTable7_soilN_top0p5_GWAS_SNPs.csv', index=False)

    aridity = load_with_phenotype('aridity_index', TRAITS['aridity_index'])
    aridity.to_csv(OUTDIR / 'SuppTable8_aridity_index_top0p5_GWAS_SNPs.csv', index=False)

    amf = pd.concat([
        load_with_phenotype('AM fungal relative abundance colonization', TRAITS['AM fungal relative abundance colonization']),
        load_with_phenotype('AM fungal roots colonized', TRAITS['AM fungal roots colonized'])
    ], ignore_index=True)
    amf.to_csv(OUTDIR / 'SuppTable9_AMF_top0p5_GWAS_SNPs.csv', index=False)

    summary = pd.DataFrame([
        {'Supplementary_table': 'SuppTable5', 'Phenotypes': 'PC1, PC2, PC3', 'Rows': len(pc)},
        {'Supplementary_table': 'SuppTable6', 'Phenotypes': 'pH', 'Rows': len(ph)},
        {'Supplementary_table': 'SuppTable7', 'Phenotypes': 'soilN', 'Rows': len(soiln)},
        {'Supplementary_table': 'SuppTable8', 'Phenotypes': 'aridity_index', 'Rows': len(aridity)},
        {'Supplementary_table': 'SuppTable9', 'Phenotypes': 'AM fungal relative abundance colonization; AM fungal roots colonized', 'Rows': len(amf)},
    ])
    summary.to_csv(OUTDIR / 'SuppTable5_9_top0p5_GWAS_SNPs_summary.csv', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
