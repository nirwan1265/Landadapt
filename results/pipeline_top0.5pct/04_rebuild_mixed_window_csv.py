#!/usr/bin/env python3

import csv
import gzip
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results')
WITHPC = ROOT / 'withPC'
OUT = ROOT / 'mixed_window_csv'
ANNOT_DIR = OUT / 'annotated_no_window_csv'
ORTHO_DIR = OUT / 'orthogroups_csv'
TEX_PATH = OUT / 'cross_species_top0p5_mixed_window_sections.tex'

for d in [OUT, ANNOT_DIR, ORTHO_DIR]:
    d.mkdir(parents=True, exist_ok=True)

ANNOT_SCRIPT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/pipeline_top0.5pct/01_gene_annotate_generic.R')
GENES_FILE = Path('/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_genes.tab.gz')
OG2GENES_FILE = Path('/Users/nirwantandukar/Documents/Research/data/orthologous_genes/odb12v2_OG2genes.tab.gz')

ARAB_RICE_WINDOW = 0
OTHERS_WINDOW = 25000
TOP_FRACTION = 0.005

SPECIES_ORDER = ['Maize', 'Sorghum', 'Rice', 'Barley', 'Arabidopsis']
KEEP_25KB_SPECIES = {'Maize', 'Sorghum', 'Barley'}
REANNOTATE_SPECIES = {'Arabidopsis', 'Rice'}

TRAITS = {
    'soilN': {'label': 'soil nitrogen (0--5 cm)'},
    'pc1_worldclim': {'label': 'WorldClim PC1'},
    'pc2_worldclim': {'label': 'WorldClim PC2'},
    'pc3_worldclim': {'label': 'WorldClim PC3'},
    'ph': {'label': 'soil pH'},
    'am_rel_abundance_colonization': {'label': 'AM fungal relative abundance colonization'},
    'am_roots_colonized': {'label': 'AM fungal roots colonized'},
}

GFF = {
    'Arabidopsis': '/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff',
    'Rice': '/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3',
}

# OrthoDB taxa among our species + close Poaceae taxa used previously
ALLOWED_TAXA = {
    '3702_0',  # Arabidopsis thaliana
    '39947_0', # Oryza sativa
    '4577_0',  # Zea mays
    '4558_0',  # Sorghum bicolor
    '4513_0',  # Hordeum vulgare
    '4557_0', '4530_0', '112509_0',
}


def tex_escape(s: str) -> str:
    return str(s).replace('gene:', '').replace('_', r'\_')


def normalize_gene_id(gene_id: str) -> str:
    g = str(gene_id).strip().upper()
    g = re.sub(r'^(GENE:|TRANSCRIPT:)', '', g)
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


def parse_best_gene_p(cell: str):
    txt = str(cell).strip()
    if txt == '' or txt.upper() == 'NA':
        return ('NA', None)
    best_gene, best_p = 'NA', None
    for part in txt.split(';'):
        part = part.strip()
        if ':' not in part:
            continue
        g, p = part.split(':', 1)
        g = g.strip()
        try:
            pv = float(p)
        except Exception:
            continue
        if best_p is None or pv < best_p:
            best_gene, best_p = g, pv
    return (best_gene, best_p)


def fmt_logp(p):
    if p is None or p <= 0:
        return 'NA'
    import math
    return f"{-math.log10(p):.2f}"


def annotate_selected_snps_no_window(trait: str, species: str, selected_df: pd.DataFrame):
    """Annotate selected SNPs with 0 bp window for Arabidopsis/Rice.
    Returns dict gene->max_p among selected rows.
    """
    sub = selected_df[selected_df['species'] == species].copy()
    if sub.empty:
        return {}

    # Prepare small GWAS-like TSV for annotation script
    tmp_in_tsv = ANNOT_DIR / f'{trait}_{species}_selected_input.tsv'
    tmp_out_tsv = ANNOT_DIR / f'{trait}_{species}_annot_0bp.tsv'
    out_csv = ANNOT_DIR / f'{trait}_{species}_annot_0bp.csv'

    w = pd.DataFrame({
        'chr': sub['chr'].astype(str),
        'rs': sub['snp'].astype(str),
        'ps': pd.to_numeric(sub['pos'], errors='coerce'),
        'p_wald': pd.to_numeric(sub['p'], errors='coerce'),
    }).dropna(subset=['ps', 'p_wald'])

    w.to_csv(tmp_in_tsv, sep='\t', index=False)

    cmd = [
        'Rscript', str(ANNOT_SCRIPT),
        str(tmp_in_tsv),
        GFF[species],
        str(tmp_out_tsv),
        str(ARAB_RICE_WINDOW),
    ]
    subprocess.run(cmd, check=True)

    ann = pd.read_csv(tmp_out_tsv, sep='\t', low_memory=False)
    ann.to_csv(out_csv, index=False)

    # remove TSVs per user request (CSV only)
    if tmp_in_tsv.exists():
        tmp_in_tsv.unlink()
    if tmp_out_tsv.exists():
        tmp_out_tsv.unlink()

    if 'closest_gene' not in ann.columns or 'p_wald' not in ann.columns:
        return {}

    ann['closest_gene'] = ann['closest_gene'].astype(str).str.strip()
    ann = ann[(ann['closest_gene'] != '') & (~ann['closest_gene'].str.lower().isin(['na', 'nan', 'none']))]
    ann['p_wald'] = pd.to_numeric(ann['p_wald'], errors='coerce')
    ann = ann.dropna(subset=['p_wald'])

    gm = ann.groupby('closest_gene', as_index=False)['p_wald'].max()
    return {r['closest_gene']: float(r['p_wald']) for _, r in gm.iterrows()}


def build_global_mapping(all_trait_gene_maps):
    alias_to_species_genes = defaultdict(list)
    for trait in all_trait_gene_maps:
        for sp, g2p in all_trait_gene_maps[trait].items():
            for g in g2p.keys():
                for a in build_aliases(sp, g):
                    alias_to_species_genes[a].append((sp, g))

    species_gene_to_odb = {sp: defaultdict(set) for sp in SPECIES_ORDER}

    with gzip.open(GENES_FILE, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            tax = parts[1].strip()
            if tax not in ALLOWED_TAXA:
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
                for sp, g in alias_to_species_genes[cu]:
                    key = (sp, g)
                    if key in seen:
                        continue
                    species_gene_to_odb[sp][g].add(odb_gene)
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


def build_orthogroups_for_trait(trait, trait_species_gene_maxp, species_gene_to_odb, odb_to_ogs):
    og_species_genes = defaultdict(lambda: defaultdict(set))

    for sp in SPECIES_ORDER:
        gmap = trait_species_gene_maxp.get(sp, {})
        for g in gmap.keys():
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
                row[f'{sp}_pvalue_max'] = ';'.join([f"{g}:{trait_species_gene_maxp[sp].get(g, 'NA'):.6g}" for g in genes])
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
    cols = ['Orthogroup', 'Gene_count']
    for sp in SPECIES_ORDER:
        cols += [f'{sp}_gene', f'{sp}_pvalue_max']
    return df[cols].sort_values(['Gene_count', 'Orthogroup'], ascending=[False, True])


def main():
    print('[A] Build mixed-window gene tables (Arabidopsis/Rice 0 bp; others 25 kb)...')

    trait_gene_maps = {}  # trait -> species -> {gene:max_p}
    trait_counts = {}

    for trait, meta in TRAITS.items():
        tdir = WITHPC / trait
        selected_file = tdir / f'{trait}_top0p5_selected_snps.csv'
        gene_table_file = tdir / f'{trait}_top0p5_gene_table.csv'
        count_file = tdir / f'{trait}_top0p5_and_bonf_selection_counts.csv'

        if not selected_file.exists() or not gene_table_file.exists() or not count_file.exists():
            print(f'  [skip] missing files for {trait}')
            continue

        selected = pd.read_csv(selected_file)
        old_gene_table = pd.read_csv(gene_table_file)
        counts = pd.read_csv(count_file)

        species_map = {sp: {} for sp in SPECIES_ORDER}

        # keep existing 25kb gene set for maize/sorghum/barley
        for sp in KEEP_25KB_SPECIES:
            sub = old_gene_table[old_gene_table['species'] == sp].copy()
            if sub.empty:
                continue
            sub['pvalue'] = pd.to_numeric(sub['pvalue'], errors='coerce')
            sub = sub.dropna(subset=['gene', 'pvalue'])
            gm = sub.groupby('gene', as_index=False)['pvalue'].max()
            species_map[sp] = {r['gene']: float(r['pvalue']) for _, r in gm.iterrows()}

        # re-annotate Arabidopsis + Rice with 0bp
        for sp in REANNOTATE_SPECIES:
            species_map[sp] = annotate_selected_snps_no_window(trait, sp, selected)

        # write trait gene table CSV
        rows = []
        for sp in SPECIES_ORDER:
            for g, p in species_map[sp].items():
                rows.append({'species': sp, 'gene': g, 'pvalue': p})
        gdf = pd.DataFrame(rows).sort_values(['species', 'pvalue', 'gene']) if rows else pd.DataFrame(columns=['species','gene','pvalue'])
        gdf.to_csv(OUT / f'{trait}_top0p5_gene_table_mixed_window.csv', index=False)

        # save summary counts
        top_snps = {r['species']: int(r['n_top0p5_snps']) for _, r in counts.iterrows()}
        bonf_snps = {r['species']: int(r['n_bonf_snps']) for _, r in counts.iterrows()}
        trait_counts[trait] = {
            'top_snps': top_snps,
            'bonf_snps': bonf_snps,
            'gene_counts': {sp: len(species_map[sp]) for sp in SPECIES_ORDER},
            'label': meta['label'],
        }

        # write per-trait counts CSV
        crows = []
        for sp in SPECIES_ORDER:
            crows.append({
                'trait': trait,
                'species': sp,
                'top0p5_snps': top_snps.get(sp, 0),
                'bonf_snps': bonf_snps.get(sp, 0),
                'genes_mixed_window': len(species_map[sp]),
                'annotation_window_bp': ARAB_RICE_WINDOW if sp in REANNOTATE_SPECIES else OTHERS_WINDOW,
            })
        pd.DataFrame(crows).to_csv(OUT / f'{trait}_top0p5_counts_mixed_window.csv', index=False)

        trait_gene_maps[trait] = species_map
        print(f'  [ok] {trait}: ' + ', '.join([f"{sp}={len(species_map[sp])}" for sp in SPECIES_ORDER]))

    print('[B] Build global OrthoDB mapping...')
    species_gene_to_odb, odb_to_ogs = build_global_mapping(trait_gene_maps)

    print('[C] Build orthogroup tables (CSV only)...')
    ortho_summary_rows = []

    for trait in trait_gene_maps:
        all_df = build_orthogroups_for_trait(trait, trait_gene_maps[trait], species_gene_to_odb, odb_to_ogs)
        ge2 = all_df[all_df['Gene_count'] >= 2].copy()
        all5 = all_df[all_df['Gene_count'] == 5].copy()

        all_path = ORTHO_DIR / f'{trait}_top0p5_orthogroups_all_mixed_window.csv'
        ge2_path = ORTHO_DIR / f'{trait}_top0p5_orthogroups_common_ge2_mixed_window.csv'
        all5_path = ORTHO_DIR / f'{trait}_top0p5_orthogroups_common5_mixed_window.csv'

        all_df.to_csv(all_path, index=False)
        ge2.to_csv(ge2_path, index=False)
        all5.to_csv(all5_path, index=False)

        ortho_summary_rows.append({
            'trait': trait,
            'orthogroups_total': len(all_df),
            'orthogroups_common_ge2': len(ge2),
            'orthogroups_common5': len(all5),
        })

        trait_counts[trait]['orthogroups_total'] = len(all_df)
        trait_counts[trait]['orthogroups_common5'] = len(all5)

        print(f'  [ok] {trait}: total={len(all_df)}, common5={len(all5)}')

    pd.DataFrame(ortho_summary_rows).to_csv(OUT / 'orthogroups_summary_mixed_window.csv', index=False)

    print('[D] Build LaTeX sections...')
    lines = []
    lines.append('% Auto-generated mixed-window report (Arabidopsis/Rice 0 bp; others 25 kb)')
    lines.append('')

    # Trait order for manuscript text
    trait_order = ['soilN', 'pc1_worldclim', 'pc2_worldclim', 'pc3_worldclim', 'ph', 'am_rel_abundance_colonization', 'am_roots_colonized']

    for trait in trait_order:
        if trait not in trait_counts:
            continue
        info = trait_counts[trait]
        label = info['label']

        lines.append(f"\\subsection{{Cross-species genome-wide associations with {label} identify shared orthologous genes across all species}}")
        lines.append('')
        lines.append(
            'We performed GWAS in Arabidopsis, barley, rice, maize, and sorghum using a mixed linear model and kinship correction in GEMMA. '
            f"Under Bonferroni correction, the number of significant SNPs was heterogeneous across species (Arabidopsis {info['bonf_snps'].get('Arabidopsis',0)}, rice {info['bonf_snps'].get('Rice',0)}, sorghum {info['bonf_snps'].get('Sorghum',0)}, barley {info['bonf_snps'].get('Barley',0)}, maize {info['bonf_snps'].get('Maize',0)}). "
            'To evaluate cross-species convergence beyond the Bonferroni tail, we selected the top 0.5\\% of SNPs per species and mapped nearby genes to orthogroups, using gene-body-only annotation (0 bp) for Arabidopsis and rice and \\pm25 kb annotation for maize, sorghum, and barley. '
            f"This yielded {info['gene_counts'].get('Arabidopsis',0):,} Arabidopsis genes, {info['gene_counts'].get('Rice',0):,} rice genes, {info['gene_counts'].get('Maize',0):,} maize genes, {info['gene_counts'].get('Barley',0):,} barley genes, and {info['gene_counts'].get('Sorghum',0):,} sorghum genes."
        )
        lines.append('')
        lines.append(
            f"Orthogroup integration identified {info.get('orthogroups_total',0):,} total orthogroups, with {info.get('orthogroups_common5',0):,} orthogroups represented in all five species. "
            f"Table~\\ref{{tab:shared_function_candidates_{trait}_mixed}} summarizes these shared orthogroups and reports, for each species, the strongest gene-level signal within each orthogroup as $-\\log_{{10}}(P)$."
        )
        lines.append('')

        common5_file = ORTHO_DIR / f'{trait}_top0p5_orthogroups_common5_mixed_window.csv'
        cdf = pd.read_csv(common5_file) if common5_file.exists() else pd.DataFrame()

        lines.append('\\begin{table*}[t]')
        lines.append(f"\\caption{{Shared orthogroups across all five species from the top 0.5\\% of SNPs for {label} (mixed-window annotation).\\label{{tab:shared_function_candidates_{trait}_mixed}}}}")
        lines.append('\\scriptsize')
        lines.append('\\tabcolsep=3pt')
        lines.append('\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lllllllllll@{}}')
        lines.append('\\toprule')
        lines.append('Orthogroup & Maize gene & $-\\log_{10}(P)$ & Sorghum gene & $-\\log_{10}(P)$ & Rice gene & $-\\log_{10}(P)$ & Barley gene & $-\\log_{10}(P)$ & Arabidopsis gene & $-\\log_{10}(P)$ \\\\')
        lines.append('\\midrule')

        if not cdf.empty:
            for _, r in cdf.sort_values('Orthogroup').iterrows():
                row = [tex_escape(r['Orthogroup'])]
                for sp in SPECIES_ORDER:
                    g, p = parse_best_gene_p(r.get(f'{sp}_pvalue_max', 'NA'))
                    if g == 'NA':
                        rawg = str(r.get(f'{sp}_gene', 'NA'))
                        g = rawg.split(';')[0] if rawg and rawg.upper() != 'NA' else 'NA'
                    row.append(tex_escape(g))
                    row.append(fmt_logp(p))
                lines.append(' & '.join(row) + ' \\\\')
        else:
            lines.append('NA & NA & NA & NA & NA & NA & NA & NA & NA & NA & NA \\\\')

        lines.append('\\botrule')
        lines.append('\\end{tabular*}')
        lines.append('\\begin{tablenotes}%')
        lines.append('\\item Note: One representative gene per species is shown for each orthogroup (the gene with the smallest $P$ within that species/orthogroup).')
        lines.append('\\end{tablenotes}')
        lines.append('\\end{table*}')
        lines.append('')

    TEX_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print(f'[E] Wrote LaTeX: {TEX_PATH}')
    print(f'[E] Output directory: {OUT}')


if __name__ == '__main__':
    main()
