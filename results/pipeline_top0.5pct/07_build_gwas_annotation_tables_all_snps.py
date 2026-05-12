#!/usr/bin/env python3
from pathlib import Path
import subprocess
import pandas as pd
import re

ROOT = Path('/Users/nirwantandukar/Documents/Github/Landadapt')
OUTDIR = ROOT / 'results' / 'tables' / 'GWAS_tables'
TMPDIR = OUTDIR / 'annotated_tmp'
OUTDIR.mkdir(parents=True, exist_ok=True)
TMPDIR.mkdir(parents=True, exist_ok=True)

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

TRAIT_LABELS = {
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


def normalize_chr(x):
    s = str(x).strip()
    s = re.sub(r'^9311_', '', s)
    s = re.sub(r'^(chromosome|chrom|chr)', '', s, flags=re.I)
    if re.fullmatch(r'[0-9]+', s):
        s = str(int(s))
    return s


def extract_attr(attr, key):
    m = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', attr)
    return m.group(1) if m else None


def read_gene_ranges(gff_path):
    records = []
    with open(gff_path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'gene':
                continue
            chr_raw = parts[0]
            try:
                start = int(parts[3]); end = int(parts[4])
            except Exception:
                continue
            attr = parts[8]
            name = extract_attr(attr, 'Name')
            gene_id = extract_attr(attr, 'gene_id')
            gid = extract_attr(attr, 'ID')
            if gid and gid.startswith('gene:'):
                gid = gid[len('gene:'):]
            gene = name or gene_id or gid
            if not gene:
                continue
            records.append((normalize_chr(chr_raw), gene, start, end))
    df = pd.DataFrame(records, columns=['chr_norm','gene','start','end'])
    if df.empty:
        return {}
    # collapse duplicates per chr+gene to broad span
    g = df.groupby(['chr_norm','gene'], as_index=False).agg(start=('start','min'), end=('end','max'))
    return {(r['chr_norm'], r['gene']):(int(r['start']), int(r['end'])) for _,r in g.iterrows()}


def list_gwas_files():
    files = []
    base = Path('/Users/nirwantandukar/Documents/Research/results/GWAS/landadapt')
    for spdir in sorted(base.iterdir()):
        if not spdir.is_dir():
            continue
        sp = spdir.name.lower()
        for f in sorted(spdir.glob('*.txt')):
            files.append((sp, f))
    ar_base = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/aridity_index')
    for f in sorted(ar_base.glob('*_aridity_index.txt')):
        nm = f.name.lower()
        if nm.startswith('at_'):
            sp = 'arabidopsis'
        elif nm.startswith('rice_'):
            sp = 'rice'
        elif nm.startswith('maize_'):
            sp = 'maize'
        elif nm.startswith('sorghum_'):
            sp = 'sorghum'
        elif nm.startswith('barley_'):
            sp = 'barley'
        else:
            continue
        files.append((sp, f))
    # dedupe by species+basename priority to aridity folder for aridity trait
    seen = {}
    for sp,f in files:
        key = (sp, f.name)
        seen[key] = f
    out = [(sp, path) for (sp,_), path in sorted(seen.items(), key=lambda x:(x[0][0], x[0][1]))]
    return out


def parse_trait(fname):
    b = fname.replace('.txt','')
    if b.startswith('AT_'):
        trait = b[len('AT_'):]
    elif b.startswith('rice_'):
        trait = b[len('rice_'):]
    elif b.startswith('maize_'):
        trait = b[len('maize_'):]
    elif b.startswith('sorghum_'):
        trait = b[len('sorghum_'):]
    elif b.startswith('barley_'):
        trait = b[len('barley_'):]
    else:
        trait = b
    return TRAIT_LABELS.get(trait, trait)


def main():
    gene_maps = {sp: read_gene_ranges(gff) for sp,gff in GFF.items()}
    all_rows = []

    files = list_gwas_files()
    manifest = []

    for sp, gwas in files:
        trait = parse_trait(gwas.name)
        window = WINDOW[sp]
        out_annot = TMPDIR / f'{sp}_{trait}.annot.tsv'

        cmd = ['Rscript', str(ANNOT), str(gwas), str(GFF[sp]), str(out_annot), str(window)]
        subprocess.run(cmd, check=True)

        df = pd.read_csv(out_annot, sep='\t', low_memory=False)

        cols = {c.lower(): c for c in df.columns}
        snp_col = cols.get('rs') or cols.get('snp')
        chr_col = cols.get('chr')
        pos_col = cols.get('ps') or cols.get('bp') or cols.get('pos')
        p_col = cols.get('p_wald') or cols.get('p')

        if not all([snp_col, chr_col, pos_col, p_col, 'closest_gene' in df.columns]):
            continue

        sub = pd.DataFrame({
            'Species': sp,
            'Trait': trait,
            'Gene': df['closest_gene'].astype(str),
            'SNP': df[snp_col].astype(str),
            'Chr': df[chr_col].astype(str),
            'Pos': pd.to_numeric(df[pos_col], errors='coerce'),
            'Pvalue': pd.to_numeric(df[p_col], errors='coerce'),
            'Distance_bp': pd.to_numeric(df.get('closest_gene_distance_bp', pd.NA), errors='coerce')
        })

        sub['Gene'] = sub['Gene'].str.strip()
        sub = sub[(sub['Gene'] != '') & (~sub['Gene'].str.lower().isin(['na','nan','none']))]
        sub = sub.dropna(subset=['Pos','Pvalue'])

        # relation by coordinate (not strand-aware)
        rels = []
        for _, r in sub.iterrows():
            key = (normalize_chr(r['Chr']), r['Gene'])
            gr = gene_maps[sp].get(key)
            if gr is None:
                rels.append(('unknown', r['Distance_bp']))
                continue
            gs, ge = gr
            pos = int(r['Pos'])
            if gs <= pos <= ge:
                rels.append(('in_gene', 0))
            elif pos < gs:
                rels.append(('upstream', gs - pos))
            else:
                rels.append(('downstream', pos - ge))

        sub['Relation'] = [x[0] for x in rels]
        sub['Distance_bp'] = [int(x[1]) if pd.notna(x[1]) else pd.NA for x in rels]

        out_csv = OUTDIR / f'{sp}_{trait}_GWAS_annotation_table.csv'
        sub = sub[['Species','Trait','Gene','SNP','Pvalue','Relation','Distance_bp','Chr','Pos']]
        sub.to_csv(out_csv, index=False)

        all_rows.append(sub)
        manifest.append({'species': sp, 'trait': trait, 'gwas_file': str(gwas), 'rows': len(sub), 'output': str(out_csv)})

    if all_rows:
        master = pd.concat(all_rows, ignore_index=True)
        master.to_csv(OUTDIR / 'all_species_all_traits_GWAS_annotation_table.csv', index=False)

        summary = master.groupby(['Trait','Species'], as_index=False).agg(
            n_rows=('SNP','count'),
            n_genes=('Gene','nunique'),
            min_p=('Pvalue','min')
        )
        summary.to_csv(OUTDIR / 'all_species_all_traits_GWAS_annotation_summary.csv', index=False)

    pd.DataFrame(manifest).to_csv(OUTDIR / 'manifest_GWAS_annotation_tables.csv', index=False)


if __name__ == '__main__':
    main()
