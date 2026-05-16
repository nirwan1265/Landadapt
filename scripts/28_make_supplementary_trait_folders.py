import math
import shutil
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SUPP = BASE / 'results' / 'tables' / 'Supplementary'
ORTH = pd.read_csv(SUPP / 'SuppTable3A_orthogroups_atleast_2_desc_noPC.csv')
SUMM = pd.read_csv(SUPP / 'SuppTable3B_orthogroups_summary.csv')

TRAIT_GROUPS = {
    'PC': {
        'phenotypes': ['PC1', 'PC2', 'PC3'],
        'raw_files': {
            'PC1': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_PC1.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_PC1.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_PC1.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_PC1.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_PC1.txt'),
            },
            'PC2': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_PC2.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_PC2.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_PC2.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_PC2.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_PC2.txt'),
            },
            'PC3': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_PC3.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_PC3.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_PC3.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_PC3.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_PC3.txt'),
            },
        },
        'orth_name': 'SuppTable3A_PC_orthogroups_atleast_2_desc.csv',
        'sum_name': 'SuppTable3B_PC_orthogroups_summary.csv',
        'gwas_name': 'SuppTable5_WorldClim_PC1_PC2_PC3_top0p5_GWAS_SNPs.csv',
        'extra_copy': [
            'SuppTable1_WorldClim_BIO_variables.csv',
            'SuppTable10_WorldClim_PC1_PC2_PC3_top3_loci_per_species.csv',
            'SuppTable10_WorldClim_PC1_PC2_PC3_top3_loci_per_species.tex',
        ],
    },
    'pH': {
        'phenotypes': ['pH'],
        'raw_files': {
            'pH': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_ph.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_ph.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_ph.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_ph.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_ph.txt'),
            },
        },
        'orth_name': 'SuppTable3A_pH_orthogroups_atleast_2_desc.csv',
        'sum_name': 'SuppTable3B_pH_orthogroups_summary.csv',
        'gwas_name': 'SuppTable6_pH_top0p5_GWAS_SNPs.csv',
        'extra_copy': [
            'SuppTable6_pH_top3_loci_per_species.csv',
            'SuppTable6_pH_top3_loci_per_species.tex',
        ],
    },
    'soilN': {
        'phenotypes': ['soilN'],
        'raw_files': {
            'soilN': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_soilN.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_soilN.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_soilN.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_soilN.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_soilN.txt'),
            },
        },
        'orth_name': 'SuppTable3A_soilN_orthogroups_atleast_2_desc.csv',
        'sum_name': 'SuppTable3B_soilN_orthogroups_summary.csv',
        'gwas_name': 'SuppTable7_soilN_top0p5_GWAS_SNPs.csv',
        'extra_copy': [
            'SuppTable7_soilN_top3_loci_per_species.csv',
            'SuppTable7_soilN_top3_loci_per_species.tex',
        ],
    },
    'AMF': {
        'phenotypes': ['am_rel_abundance_colonization', 'am_roots_colonized'],
        'raw_files': {
            'AM fungal relative abundance colonization': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_am_rel_abundance_colonization.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_am_rel_abundance_colonization.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_am_rel_abundance_colonization.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_am_rel_abundance_colonization.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_am_rel_abundance_colonization.txt'),
            },
            'AM fungal roots colonized': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_am_roots_colonized.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_am_roots_colonized.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_am_roots_colonized.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_am_roots_colonized.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_am_roots_colonized.txt'),
            },
        },
        'orth_name': 'SuppTable3A_AMF_orthogroups_atleast_2_desc.csv',
        'sum_name': 'SuppTable3B_AMF_orthogroups_summary.csv',
        'gwas_name': 'SuppTable8_AMF_top0p5_GWAS_SNPs.csv',
        'extra_copy': [
            'SuppTable9_AMF_top3_loci_per_species.csv',
            'SuppTable9_AMF_top3_loci_per_species.tex',
        ],
    },
    'aridity': {
        'phenotypes': ['aridity_index'],
        'raw_files': {
            'aridity_index': {
                'Arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis/AT_aridity_index.txt'),
                'Barley': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/barley/barley_aridity_index.txt'),
                'Maize': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize/maize_aridity_index.txt'),
                'Rice': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice/rice_aridity_index.txt'),
                'Sorghum': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum/sorghum_aridity_index.txt'),
            },
        },
        'orth_name': 'SuppTable3A_aridity_orthogroups_atleast_2_desc.csv',
        'sum_name': 'SuppTable3B_aridity_orthogroups_summary.csv',
        'gwas_name': 'SuppTable9_aridity_index_top0p5_GWAS_SNPs.csv',
        'extra_copy': [
            'SuppTable8_aridity_top3_loci_per_species.csv',
            'SuppTable8_aridity_top3_loci_per_species.tex',
        ],
    },
}


def top05_from_raw(path):
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={'rs': 'snp', 'ps': 'pos', 'p_wald': 'p'})
    df['p'] = pd.to_numeric(df['p'], errors='coerce')
    df = df.dropna(subset=['p']).copy()
    n_keep = max(1, math.ceil(len(df) * 0.005))
    top = df.nsmallest(n_keep, 'p').copy()
    top['logp'] = -top['p'].map(lambda x: math.log10(x) if x > 0 else float('inf'))
    return top[['snp', 'chr', 'pos', 'p', 'logp']]


def build_group(group_name, cfg):
    outdir = SUPP / group_name
    outdir.mkdir(parents=True, exist_ok=True)

    orth = ORTH[ORTH['Phenotype'].isin(cfg['phenotypes'])].copy()
    orth.to_csv(outdir / cfg['orth_name'], index=False)

    summ = SUMM[SUMM['Phenotype'].isin(cfg['phenotypes'])].copy()
    summ.to_csv(outdir / cfg['sum_name'], index=False)

    frames = []
    for phenotype, species_files in cfg['raw_files'].items():
        for species, path in species_files.items():
            top = top05_from_raw(path)
            top.insert(0, 'species', species)
            top.insert(0, 'Phenotype', phenotype)
            frames.append(top)
    gwas = pd.concat(frames, ignore_index=True)
    gwas.to_csv(outdir / cfg['gwas_name'], index=False)

    for name in cfg['extra_copy']:
        src = SUPP / name
        if src.exists():
            shutil.copy2(src, outdir / name)

    return {
        'folder': outdir.name,
        'orth_rows': len(orth),
        'gwas_rows': len(gwas),
        'files': len(list(outdir.iterdir())),
    }

summary_rows = []
for group_name, cfg in TRAIT_GROUPS.items():
    summary_rows.append(build_group(group_name, cfg))

print(pd.DataFrame(summary_rows).to_string(index=False))
