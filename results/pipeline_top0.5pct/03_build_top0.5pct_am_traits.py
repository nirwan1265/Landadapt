#!/usr/bin/env python3

import csv
import gzip
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT_OUT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/top0.5pct_traits_am')
ANNOT_DIR = ROOT_OUT / 'annotated_25000bp'
GO_DIR = ROOT_OUT / 'go_inputs'
ORTHO_DIR = ROOT_OUT / 'orthogroups'
for d in [ANNOT_DIR, GO_DIR, ORTHO_DIR]:
    d.mkdir(parents=True, exist_ok=True)

ANNOT_SCRIPT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/pipeline_top0.5pct/01_gene_annotate_generic.R')

GENES_FILE = Path('/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_genes.tab.gz')
OG2GENES_FILE = Path('/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_OG2genes.tab.gz')
TOP_FRAC = 0.005
WINDOW_BP = 25000

BASE = '/Users/nirwantandukar/Documents/Github/Landadapt/results/AT_am_rel_abundance_colonization'

TRAITS = {
    'am_rel_abundance_colonization': {
        'Arabidopsis': f'{BASE}/AT_am_rel_abundance_colonization.txt',
        'Rice': f'{BASE}/rice_am_rel_abundance_colonization.txt',
        'Maize': f'{BASE}/maize_am_rel_abundance_colonization.txt',
        'Sorghum': f'{BASE}/sorghum_am_rel_abundance_colonization.txt',
        'Barley': f'{BASE}/barley_am_rel_abundance_colonization.txt',
    },
    'am_roots_colonized': {
        'Arabidopsis': f'{BASE}/AT_am_roots_colonized.txt',
        'Rice': f'{BASE}/rice_am_roots_colonized.txt',
        'Maize': f'{BASE}/maize_am_roots_colonized.txt',
        'Sorghum': f'{BASE}/sorghum_am_roots_colonized.txt',
        'Barley': f'{BASE}/barley_am_roots_colonized.txt',
    },
}

GFF = {
    'Arabidopsis': '/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff',
    'Rice': '/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3',
    'Maize': '/Users/nirwantandukar/Documents/Research/data/maize_gene_annotation/ENSEMBLE_Zea_mays.Zm-B73-REFERENCE-NAM-5.0.60.chr.gff3',
    'Sorghum': '/Users/nirwantandukar/Documents/Research/data/sorghum_annotation/gene_annotation/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.54.gff3',
    'Barley': '/Users/nirwantandukar/Documents/Research/results/Barley/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3',
}

SPECIES_ORDER = ['Arabidopsis', 'Rice', 'Maize', 'Sorghum', 'Barley']


def run_annotation(gwas_file: str, gff_file: str, out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = ['Rscript', str(ANNOT_SCRIPT), gwas_file, gff_file, str(out_file), str(WINDOW_BP)]
    subprocess.run(cmd, check=True)


def normalize_gene_id(gene_id: str) -> str:
    g = str(gene_id).strip().upper()
    g = re.sub(r'_T\d+$', '', g)
    g = re.sub(r'_P\d+$', '', g)
    g = re.sub(r'\.[0-9]+$', '', g)
    g = re.sub(r'-[A-Z]{1,3}$', '', g)
    return g


def build_aliases(species: str, gene_id: str) -> set[str]:
    aliases = {str(gene_id).strip().upper()}
    norm = normalize_gene_id(gene_id)
    aliases.add(norm)
    if species == 'Rice':
        m = re.match(r'^LOC_OS(\d{2})G(\d+)$', norm)
        if m:
            ch, num = m.groups()
            aliases.add(f'OS{ch}G{num}')
            aliases.add(f'OS{ch}G{num}0')
            aliases.add(f'OS{ch}G{num}00')
            aliases.add(f'OS{ch}G{num.zfill(7)}')
    return aliases


def detect_col(cands, cols):
    for c in cands:
        if c in cols:
            return c
    return None


def top_genes_from_annotated(annot_file: Path):
    df = pd.read_csv(annot_file, sep='\t', low_memory=False)
    snp_col = detect_col(['rs', 'SNP', 'snp'], df.columns)
    p_col = detect_col(['p_wald', 'P.value', 'PValue', 'P'], df.columns)
    if snp_col is None or p_col is None or 'closest_gene' not in df.columns:
        raise ValueError(f'Cannot detect required columns in {annot_file}')

    d = df[[snp_col, p_col, 'closest_gene']].copy().dropna(subset=[snp_col, p_col])
    snp_best = d.groupby(snp_col, as_index=False)[p_col].min().sort_values(p_col, ascending=True)
    n_unique = len(snp_best)
    n_top = max(1, math.ceil(TOP_FRAC * n_unique))
    top_snps = set(snp_best.head(n_top)[snp_col].astype(str))

    sub = d[d[snp_col].astype(str).isin(top_snps)].copy()
    sub['closest_gene'] = sub['closest_gene'].astype(str).str.strip()
    sub = sub[(sub['closest_gene'] != '') & (sub['closest_gene'].str.upper() != 'NA')]

    gene_maxp = sub.groupby('closest_gene', as_index=False)[p_col].max()
    genes = sorted(gene_maxp['closest_gene'].unique().tolist())
    g2p = dict(zip(gene_maxp['closest_gene'], gene_maxp[p_col]))
    return {'genes': genes, 'gene_maxp': g2p, 'n_unique_snps': n_unique, 'n_top_snps': n_top}


def build_global_mapping(all_trait_species_genes):
    alias_to_species_genes = defaultdict(list)
    for trait in all_trait_species_genes:
        for species, genes in all_trait_species_genes[trait].items():
            for g in genes:
                for a in build_aliases(species, g):
                    alias_to_species_genes[a].append((species, g))

    species_gene_to_odb = {sp: defaultdict(set) for sp in SPECIES_ORDER}

    with gzip.open(GENES_FILE, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            odb_gene = parts[0].strip()
            candidates = [p.strip() for p in parts[2:7] if p.strip()]
            desc = parts[7]
            candidates.extend(re.findall(r'\b(?:LOC_Os\d{2}g\d+|Os\d{2}g\d{4,})\b', desc, flags=re.I))

            seen = set()
            for cand in candidates:
                cu = cand.upper()
                if cu not in alias_to_species_genes:
                    continue
                for sp, original_gene in alias_to_species_genes[cu]:
                    key = (sp, original_gene)
                    if key in seen:
                        continue
                    species_gene_to_odb[sp][original_gene].add(odb_gene)
                    seen.add(key)

    matched_odb = set()
    for sp in species_gene_to_odb:
        for g in species_gene_to_odb[sp]:
            matched_odb.update(species_gene_to_odb[sp][g])

    odb_to_ogs = defaultdict(set)
    with gzip.open(OG2GENES_FILE, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            og, odb_gene = parts[0], parts[1]
            if odb_gene in matched_odb:
                odb_to_ogs[odb_gene].add(og)

    return species_gene_to_odb, odb_to_ogs


def build_trait_orthogroups(trait_gene_sets, trait_gene_pmax, species_gene_to_odb, odb_to_ogs):
    og_species_genes = defaultdict(lambda: defaultdict(set))

    for sp in SPECIES_ORDER:
        for g in trait_gene_sets[sp]:
            for odb in species_gene_to_odb.get(sp, {}).get(g, set()):
                for og in odb_to_ogs.get(odb, set()):
                    og_species_genes[og][sp].add(g)

    rows = []
    for og, sp_map in og_species_genes.items():
        row = {'Orthogroup': og}
        gc = 0
        for sp in SPECIES_ORDER:
            genes = sorted(sp_map.get(sp, set()))
            if genes:
                row[f'{sp}_gene'] = ';'.join(genes)
                row[f'{sp}_pvalue_max'] = ';'.join([f"{g}:{trait_gene_pmax[sp].get(g, 'NA')}" for g in genes])
                gc += 1
            else:
                row[f'{sp}_gene'] = 'NA'
                row[f'{sp}_pvalue_max'] = 'NA'
        row['Gene_count'] = gc
        rows.append(row)

    if not rows:
        cols = ['Orthogroup', 'Gene_count']
        for sp in SPECIES_ORDER:
            cols += [f'{sp}_gene', f'{sp}_pvalue_max']
        return pd.DataFrame(columns=cols)

    df = pd.DataFrame(rows)
    ordered_cols = ['Orthogroup', 'Gene_count']
    for sp in SPECIES_ORDER:
        ordered_cols += [f'{sp}_gene', f'{sp}_pvalue_max']
    return df[ordered_cols].sort_values(['Gene_count', 'Orthogroup'], ascending=[False, True])


def main():
    print('[A] Annotating GWAS files...')
    annotated_paths = defaultdict(dict)
    for trait, sp_map in TRAITS.items():
        for sp, gwas in sp_map.items():
            out_f = ANNOT_DIR / trait / f'{sp}_{trait}.annot_25000bp.tsv'
            run_annotation(gwas, GFF[sp], out_f)
            annotated_paths[trait][sp] = out_f

    print('[B] Extract top 0.5% SNP genes per trait/species...')
    trait_species_genes = defaultdict(dict)
    trait_species_gene_pmax = defaultdict(dict)

    for trait in TRAITS:
        summary = [f'trait={trait}', f'top_fraction={TOP_FRAC}']
        species_lists = {}
        for sp in SPECIES_ORDER:
            info = top_genes_from_annotated(annotated_paths[trait][sp])
            trait_species_genes[trait][sp] = info['genes']
            trait_species_gene_pmax[trait][sp] = info['gene_maxp']
            species_lists[sp] = info['genes']
            summary.append(f"{sp}: unique_snps={info['n_unique_snps']}, top0.5pct_snps={info['n_top_snps']}, genes={len(info['genes'])}")

        max_len = max(len(v) for v in species_lists.values())
        pad = lambda v: v + [''] * (max_len - len(v))
        go_df = pd.DataFrame({
            'Arabidopsis_gene': pad(species_lists['Arabidopsis']),
            'Rice_gene': pad(species_lists['Rice']),
            'Maize_gene': pad(species_lists['Maize']),
            'Sorghum_gene': pad(species_lists['Sorghum']),
            'Barley_gene': pad(species_lists['Barley']),
        })
        go_df.to_csv(GO_DIR / f'go_input_top0.5pct_{trait}_by_species.csv', index=False)
        (GO_DIR / f'go_input_top0.5pct_{trait}_by_species.summary.txt').write_text('\n'.join(summary) + '\n', encoding='utf-8')

    print('[C] Build global OrthoDB mapping once...')
    species_gene_to_odb, odb_to_ogs = build_global_mapping(trait_species_genes)

    print('[D] Build orthogroup/common tables per trait...')
    for trait in TRAITS:
        all_df = build_trait_orthogroups(
            trait_species_genes[trait],
            trait_species_gene_pmax[trait],
            species_gene_to_odb,
            odb_to_ogs,
        )

        all_path = ORTHO_DIR / f'orthogroups_top0.5pct_{trait}_all.tsv'
        ge2_path = ORTHO_DIR / f'orthogroups_top0.5pct_{trait}_common_ge2.tsv'
        all5_path = ORTHO_DIR / f'orthogroups_top0.5pct_{trait}_common_all5.tsv'
        sum_path = ORTHO_DIR / f'orthogroups_top0.5pct_{trait}.summary.txt'

        all_df.to_csv(all_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)
        ge2 = all_df[all_df['Gene_count'] >= 2].copy()
        ge2.to_csv(ge2_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)
        all5 = all_df[all_df['Gene_count'] == 5].copy()
        all5.to_csv(all5_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)

        summary = [
            f'trait={trait}',
            f'top_fraction={TOP_FRAC}',
            f'orthogroups_total={len(all_df)}',
            f'orthogroups_common_ge2={len(ge2)}',
            f'orthogroups_common_all5={len(all5)}',
        ]
        for sp in SPECIES_ORDER:
            summary.append(f'{sp}_top_genes={len(trait_species_genes[trait][sp])}')
        sum_path.write_text('\n'.join(summary) + '\n', encoding='utf-8')

    print('DONE')


if __name__ == '__main__':
    main()
