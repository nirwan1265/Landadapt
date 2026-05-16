import math
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary' / 'PC'
OUTDIR.mkdir(parents=True, exist_ok=True)

SPECIES_FILES = {
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
}

rows = []
summary = []
for phenotype, files in SPECIES_FILES.items():
    for species, path in files.items():
        df = pd.read_csv(path, sep='\t')
        df = df.rename(columns={'rs': 'snp', 'ps': 'pos', 'p_wald': 'p', 'chr': 'chr'})
        df['p'] = pd.to_numeric(df['p'], errors='coerce')
        df = df.dropna(subset=['p']).copy()
        n_keep = max(1, math.ceil(len(df) * 0.005))
        top = df.nsmallest(n_keep, 'p').copy()
        top.insert(0, 'Phenotype', phenotype)
        top.insert(1, 'species', species)
        top['logp'] = -top['p'].map(lambda x: math.log10(x) if x > 0 else float('inf'))
        top = top[['Phenotype', 'species', 'snp', 'chr', 'pos', 'p', 'logp']]
        rows.append(top)
        summary.append({'Phenotype': phenotype, 'species': species, 'Total_SNPs': len(df), 'Top0.5pct_SNPs': len(top)})

out = pd.concat(rows, ignore_index=True)
out = out.sort_values(['Phenotype', 'species', 'p', 'snp']).reset_index(drop=True)
out.to_csv(OUTDIR / 'SuppTable5_WorldClim_PC1_PC2_PC3_top0p5_GWAS_SNPs_from_landadapt_raw.csv', index=False)
pd.DataFrame(summary).to_csv(OUTDIR / 'SuppTable5_WorldClim_PC1_PC2_PC3_top0p5_GWAS_SNPs_from_landadapt_raw_summary.csv', index=False)
print(out.head(10).to_string(index=False))
print('rows', len(out))
