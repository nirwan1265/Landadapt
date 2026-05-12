#!/usr/bin/env python3
from pathlib import Path
import subprocess
import pandas as pd
import re
import math

ROOT = Path('/Users/nirwantandukar/Documents/Github/Landadapt')
OUT_BASE = ROOT / 'results' / 'tables' / 'GWAS_tables'
OUT_DIR = OUT_BASE / 'top0p5_significant'
TMP_DIR = OUT_DIR / 'tmp'
OUT_DIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)

ANNOT = ROOT / 'results' / 'pipeline_top0.5pct' / '01_gene_annotate_generic.R'

GFF = {
    'arabidopsis': Path('/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff'),
    'rice': Path('/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3'),
    'maize': Path('/Users/nirwantandukar/Documents/Research/data/maize_gene_annotation/ENSEMBLE_Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.chr.gff3'),
    'sorghum': Path('/Users/nirwantandukar/Documents/Research/data/sorghum_annotation/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3'),
    'barley': Path('/Users/nirwantandukar/Documents/Research/results/Barley/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3'),
}

WINDOW = {
    'arabidopsis': 10000,
    'rice': 10000,
    'maize': 25000,
    'sorghum': 25000,
    'barley': 25000,
}

TRAIT_MAP = {
    'PC1': 'pc1_worldclim',
    'PC2': 'pc2_worldclim',
    'PC3': 'pc3_worldclim',
    'soilN': 'soilN',
    'ph': 'ph',
    'cec': 'cec',
    'am_rel_abundance_colonization': 'am_rel_abundance_colonization',
    'am_roots_colonized': 'am_roots_colonized',
    'aridity_index': 'aridity_index',
}

TOP_FRAC = 0.005


def norm_chr(x):
    s = str(x).strip()
    s = re.sub(r'^9311_', '', s)
    s = re.sub(r'^(chromosome|chrom|chr)', '', s, flags=re.I)
    if re.fullmatch(r'[0-9]+', s):
        s = str(int(s))
    return s


def attr_get(attr, key):
    m = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', attr)
    return m.group(1) if m else None


def read_gene_ranges(gff_path):
    rows = []
    with open(gff_path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9 or p[2] != 'gene':
                continue
            try:
                start = int(p[3]); end = int(p[4])
            except Exception:
                continue
            attr = p[8]
            name = attr_get(attr, 'Name')
            gene_id = attr_get(attr, 'gene_id')
            gid = attr_get(attr, 'ID')
            if gid and gid.startswith('gene:'):
                gid = gid[5:]
            gene = name or gene_id or gid
            if not gene:
                continue
            rows.append((norm_chr(p[0]), gene, start, end))
    df = pd.DataFrame(rows, columns=['chr_norm', 'gene', 'start', 'end'])
    if df.empty:
        return {}
    gr = df.groupby(['chr_norm', 'gene'], as_index=False).agg(start=('start', 'min'), end=('end', 'max'))
    return {(r['chr_norm'], r['gene']): (int(r['start']), int(r['end'])) for _, r in gr.iterrows()}


def find_col(cols, candidates):
    low = {c.lower(): c for c in cols}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


def parse_trait_and_species(path):
    name = path.name.replace('.txt', '')
    l = name.lower()
    if l.startswith('at_'):
        sp = 'arabidopsis'
        tr = name[3:]
    elif l.startswith('rice_'):
        sp = 'rice'
        tr = name[5:]
    elif l.startswith('maize_'):
        sp = 'maize'
        tr = name[6:]
    elif l.startswith('sorghum_'):
        sp = 'sorghum'
        tr = name[8:]
    elif l.startswith('barley_'):
        sp = 'barley'
        tr = name[7:]
    else:
        return None, None
    trait = TRAIT_MAP.get(tr, tr)
    return sp, trait


def list_gwas_files():
    files = []
    base = Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt')
    for spdir in sorted(base.iterdir()):
        if not spdir.is_dir():
            continue
        for f in sorted(spdir.glob('*.txt')):
            sp, trait = parse_trait_and_species(f)
            if sp and trait:
                files.append((sp, trait, f))

    # override/add aridity from the explicit folder user gave
    ar = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index')
    for f in sorted(ar.glob('*_aridity_index.txt')):
        sp, trait = parse_trait_and_species(f)
        if sp and trait:
            files = [x for x in files if not (x[0] == sp and x[1] == trait)]
            files.append((sp, trait, f))

    # stable order
    files.sort(key=lambda x: (x[1], x[0]))
    return files


def top0p5_snps(gwas_file):
    d = pd.read_csv(gwas_file, sep='\t', low_memory=False)
    snp_col = find_col(d.columns, ['rs', 'SNP', 'snp'])
    chr_col = find_col(d.columns, ['chr', 'Chr', 'CHR', 'chrom', 'CHROM'])
    pos_col = find_col(d.columns, ['ps', 'pos', 'Pos', 'BP', 'bp'])
    p_col = find_col(d.columns, ['p_wald', 'P', 'pvalue', 'P.value'])
    if None in [snp_col, chr_col, pos_col, p_col]:
        return pd.DataFrame(columns=['chr', 'rs', 'ps', 'p_wald'])

    x = pd.DataFrame({
        'rs': d[snp_col].astype(str),
        'chr': d[chr_col].astype(str),
        'ps': pd.to_numeric(d[pos_col], errors='coerce'),
        'p_wald': pd.to_numeric(d[p_col], errors='coerce')
    })
    x = x.dropna(subset=['chr', 'ps', 'p_wald'])
    x = x[(x['chr'] != '') & (x['p_wald'] > 0)]
    x = x.groupby(['rs', 'chr', 'ps'], as_index=False)['p_wald'].min()

    n = len(x)
    if n == 0:
        return pd.DataFrame(columns=['chr', 'rs', 'ps', 'p_wald'])
    k = max(1, math.floor(n * TOP_FRAC))
    cutoff = x['p_wald'].sort_values().iloc[k - 1]
    sel = x[x['p_wald'] <= cutoff].copy()
    return sel.sort_values(['p_wald', 'chr', 'ps'])


def build_relation_rows(ann_df, species, trait, gene_ranges):
    cols = {c.lower(): c for c in ann_df.columns}
    snp_col = cols.get('rs')
    chr_col = cols.get('chr')
    pos_col = cols.get('ps')
    p_col = cols.get('p_wald')
    if None in [snp_col, chr_col, pos_col, p_col] or 'closest_gene' not in ann_df.columns:
        return pd.DataFrame(columns=['Species', 'Trait', 'Gene', 'SNP', 'Pvalue', 'Relation', 'Distance_bp', 'Chr', 'Pos'])

    sub = pd.DataFrame({
        'Species': species,
        'Trait': trait,
        'Gene': ann_df['closest_gene'].astype(str).str.strip(),
        'SNP': ann_df[snp_col].astype(str),
        'Pvalue': pd.to_numeric(ann_df[p_col], errors='coerce'),
        'Chr': ann_df[chr_col].astype(str),
        'Pos': pd.to_numeric(ann_df[pos_col], errors='coerce'),
        'Distance_bp': pd.to_numeric(ann_df.get('closest_gene_distance_bp', pd.NA), errors='coerce')
    })

    sub = sub[(sub['Gene'] != '') & (~sub['Gene'].str.lower().isin(['na', 'nan', 'none']))]
    sub = sub.dropna(subset=['Pvalue', 'Pos'])

    rel = []
    dist = []
    gm = gene_ranges[species]
    for _, r in sub.iterrows():
        key = (norm_chr(r['Chr']), r['Gene'])
        gr = gm.get(key)
        if gr is None:
            rel.append('unknown')
            dist.append(int(r['Distance_bp']) if pd.notna(r['Distance_bp']) else pd.NA)
            continue
        gs, ge = gr
        p = int(r['Pos'])
        if gs <= p <= ge:
            rel.append('in_gene')
            dist.append(0)
        elif p < gs:
            rel.append('upstream')
            dist.append(gs - p)
        else:
            rel.append('downstream')
            dist.append(p - ge)

    sub['Relation'] = rel
    sub['Distance_bp'] = dist
    sub['Species'] = sub['Species'].str.lower()
    return sub[['Species', 'Trait', 'Gene', 'SNP', 'Pvalue', 'Relation', 'Distance_bp', 'Chr', 'Pos']].sort_values('Pvalue')


def main():
    files = list_gwas_files()
    gm = {sp: read_gene_ranges(gff) for sp, gff in GFF.items()}

    master = []
    manifest = []

    for sp, trait, gwas in files:
        sel = top0p5_snps(gwas)
        if sel.empty:
            continue
        tmp_in = TMP_DIR / f'{sp}_{trait}_top0p5.tsv'
        tmp_out = TMP_DIR / f'{sp}_{trait}_top0p5.annot.tsv'
        sel.to_csv(tmp_in, sep='\t', index=False)

        cmd = ['Rscript', str(ANNOT), str(tmp_in), str(GFF[sp]), str(tmp_out), str(WINDOW[sp])]
        subprocess.run(cmd, check=True)

        ann = pd.read_csv(tmp_out, sep='\t', low_memory=False)
        out = build_relation_rows(ann, sp, trait, gm)

        out_file = OUT_DIR / f'{sp}_{trait}_GWAS_top0p5_annotation_table.csv'
        out.to_csv(out_file, index=False)

        master.append(out)
        manifest.append({'species': sp, 'trait': trait, 'rows': len(out), 'gwas_file': str(gwas), 'output': str(out_file)})

    if master:
        all_df = pd.concat(master, ignore_index=True)
        all_df.to_csv(OUT_DIR / 'all_species_all_traits_GWAS_top0p5_annotation_table.csv', index=False)
        summary = all_df.groupby(['Trait', 'Species'], as_index=False).agg(n_rows=('SNP', 'count'), n_genes=('Gene', 'nunique'), min_p=('Pvalue', 'min'))
        summary.to_csv(OUT_DIR / 'all_species_all_traits_GWAS_top0p5_annotation_summary.csv', index=False)

    pd.DataFrame(manifest).to_csv(OUT_DIR / 'manifest_top0p5_annotation_tables.csv', index=False)


if __name__ == '__main__':
    main()
