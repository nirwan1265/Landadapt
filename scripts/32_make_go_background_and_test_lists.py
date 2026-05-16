from pathlib import Path
import pandas as pd
import numpy as np
import re

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary' / 'background_GO'
OUTDIR.mkdir(parents=True, exist_ok=True)

TRAITS = [
    ('PC1', 'pc1_worldclim'),
    ('PC2', 'pc2_worldclim'),
    ('PC3', 'pc3_worldclim'),
    ('pH', 'ph'),
    ('soilN', 'soilN'),
    ('AM_rel', 'am_rel_abundance_colonization'),
    ('AM_roots', 'am_roots_colonized'),
    ('aridity', 'aridity_index'),
]
RAW_SUFFIX = {
    'PC1': 'PC1',
    'PC2': 'PC2',
    'PC3': 'PC3',
    'pH': 'ph',
    'soilN': 'soilN',
    'AM_rel': 'am_rel_abundance_colonization',
    'AM_roots': 'am_roots_colonized',
    'aridity': 'aridity_index',
}
SPECIES = {
    'Arabidopsis': {
        'raw_prefix': 'AT',
        'gwas_dir': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/arabidopsis'),
        'gff': Path('/Users/nirwantandukar/Documents/Research/data/GFF3/Arabidopsis_thaliana.TAIR10.62.gff3'),
        'window': 10000,
        'chr_normalizer': lambda x: str(x).replace('Chr', '').replace('.0', ''),
    },
    'Rice': {
        'raw_prefix': 'rice',
        'gwas_dir': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/rice'),
        'gff': Path('/Users/nirwantandukar/Documents/Research/data/GFF3/Oryza_sativa.IRGSP-1.0.62.chr.gff3'),
        'window': 10000,
        'chr_normalizer': lambda x: str(x).replace('chr', '').replace('Chr', ''),
    },
    'Maize': {
        'raw_prefix': 'maize',
        'gwas_dir': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/maize'),
        'gff': Path('/Users/nirwantandukar/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3'),
        'window': 25000,
        'chr_normalizer': lambda x: str(x).replace('chr', '').replace('Chr', ''),
    },
    'Sorghum': {
        'raw_prefix': 'sorghum',
        'gwas_dir': Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt/sorghum'),
        'gff': Path('/Users/nirwantandukar/Documents/Research/data/sorghum_annotation/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3'),
        'window': 25000,
        'chr_normalizer': lambda x: str(x).replace('chr', '').replace('Chr', ''),
    },
}

# Barley fallback: use the empirical study-wide mapped universe since a local GFF3 was not found.
EMPIRICAL_ALL_TRAITS = [
    BASE / 'results' / 'mixed_window_10kb_csv' / 'pc1_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'pc2_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'pc3_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'ph_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'soilN_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'am_rel_abundance_colonization_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'am_roots_colonized_top0p5_gene_table_mixed_window_10kb.csv',
    BASE / 'results' / 'mixed_window_10kb_csv' / 'aridity_index_top0p5_gene_table_mixed_window_10kb.csv',
]


def parse_attr(attr, keys):
    for key in keys:
        m = re.search(rf'{key}=([^;]+)', attr)
        if m:
            return m.group(1)
    return None


def load_gff_genes(gff_path, normalize_chr):
    genes = {}
    with open(gff_path, 'r') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            seqid, _, feature, start, end, _, _, _, attrs = parts
            if feature != 'gene':
                continue
            gene_id = parse_attr(attrs, ['ID', 'gene_id', 'Name', 'gene'])
            if gene_id is None:
                continue
            chrom = normalize_chr(seqid)
            genes.setdefault(chrom, []).append((int(start), int(end), gene_id))
    for chrom in genes:
        genes[chrom].sort()
    return genes


def load_gwas_positions(gwas_path, normalize_chr):
    df = pd.read_csv(gwas_path, sep='\t', usecols=['chr', 'ps'])
    pos_by_chr = {}
    for chrom, sub in df.groupby('chr'):
        norm = normalize_chr(chrom)
        arr = np.sort(sub['ps'].astype(int).to_numpy())
        pos_by_chr[norm] = np.unique(arr)
    return pos_by_chr


def map_positions_to_genes(pos_by_chr, gene_index, window):
    bg = set()
    for chrom, genes in gene_index.items():
        positions = pos_by_chr.get(chrom)
        if positions is None or len(positions) == 0:
            continue
        for start, end, gene_id in genes:
            left = start - window
            right = end + window
            i = np.searchsorted(positions, left, side='left')
            if i < len(positions) and positions[i] <= right:
                bg.add(gene_id)
    return sorted(bg)


def write_list(path, genes):
    with open(path, 'w') as fh:
        for gene in genes:
            fh.write(f'{gene}\n')


def empirical_background(species):
    genes = set()
    for path in EMPIRICAL_ALL_TRAITS:
        df = pd.read_csv(path)
        genes.update(df.loc[df['species'] == species, 'gene'].dropna().astype(str).unique())
    return sorted(genes)


def main():
    summary = []

    # test lists from top-0.5 mapped genes
    for trait_label, trait_file in TRAITS:
        df = pd.read_csv(BASE / 'results' / 'mixed_window_10kb_csv' / f'{trait_file}_top0p5_gene_table_mixed_window_10kb.csv')
        for species in ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']:
            genes = sorted(df.loc[df['species'] == species, 'gene'].dropna().astype(str).unique())
            out = OUTDIR / f'test_{trait_label}_{species}_top0p5_genes.txt'
            write_list(out, genes)
            summary.append({'file_type': 'test', 'species': species, 'trait': trait_label, 'n_genes': len(genes), 'path': str(out)})

    # species-wide backgrounds from raw GWAS + GFF where available.
    # The tested SNP universe is shared across traits within a species, so one representative GWAS file is sufficient.
    for species, cfg in SPECIES.items():
        if species == 'Rice':
            continue
        gene_index = load_gff_genes(cfg['gff'], cfg['chr_normalizer'])
        representative = cfg['gwas_dir'] / f"{cfg['raw_prefix']}_{RAW_SUFFIX['PC1']}.txt"
        pos_by_chr = load_gwas_positions(representative, cfg['chr_normalizer'])
        bg = map_positions_to_genes(pos_by_chr, gene_index, cfg['window'])
        out = OUTDIR / f'background_{species}_all_GWAS_mappable_genes.txt'
        write_list(out, bg)
        summary.append({'file_type': 'background', 'species': species, 'trait': 'all_traits', 'n_genes': len(bg), 'path': str(out)})

    # empirical fallbacks where a matching local annotation resource is absent or ID namespaces differ
    for species in ['Barley', 'Rice']:
        bg = empirical_background(species)
        out = OUTDIR / f'background_{species}_all_GWAS_mappable_genes_fallback_from_empirical_mapping.txt'
        write_list(out, bg)
        summary.append({'file_type': 'background', 'species': species, 'trait': 'all_traits', 'n_genes': len(bg), 'path': str(out)})

    pd.DataFrame(summary).to_csv(OUTDIR / 'background_GO_file_summary.csv', index=False)
    print('WROTE', OUTDIR / 'background_GO_file_summary.csv')
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == '__main__':
    main()
