#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import pandas as pd
import math
import subprocess

ROOT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results')
WITHPC = ROOT / 'withPC'
OUT = ROOT / 'mixed_window_10kb_csv'
ANNOT_DIR = OUT / 'annotated_10kb_csv'
ORTHO_DIR = OUT / 'orthogroups_csv'
TEX_FILE = OUT / 'cross_species_top0p5_mixed_window_10kb_sections.tex'

for d in [OUT, ANNOT_DIR, ORTHO_DIR]:
    d.mkdir(parents=True, exist_ok=True)

ANNOT_SCRIPT = Path('/Users/nirwantandukar/Documents/Github/Landadapt/results/pipeline_top0.5pct/01_gene_annotate_generic.R')
WINDOW_BP = 10000

SPECIES_ORDER = ['Maize', 'Sorghum', 'Rice', 'Barley', 'Arabidopsis']
REANNOTATE = {'Arabidopsis', 'Rice'}
KEEP_25KB = {'Maize', 'Sorghum', 'Barley'}

GFF = {
    'Arabidopsis': '/Users/nirwantandukar/Documents/Research/results/Arabidopsis/gff3_ref/TAIR10_GFF3_genes.gff',
    'Rice': '/Users/nirwantandukar/Documents/Research/results/Rice_3001/gff3_ref/osa1_r7.all_models.gff3',
}

TRAITS = {
    'soilN': 'soil nitrogen (0--5 cm)',
    'pc1_worldclim': 'WorldClim PC1',
    'pc2_worldclim': 'WorldClim PC2',
    'pc3_worldclim': 'WorldClim PC3',
    'ph': 'soil pH',
    'am_rel_abundance_colonization': 'AM fungal relative abundance colonization',
    'am_roots_colonized': 'AM fungal roots colonized',
}


def split_genes(x):
    t = str(x).strip()
    if t == '' or t.upper() == 'NA' or t.lower() == 'nan':
        return []
    out = []
    for g in t.split(';'):
        g = g.strip().replace('gene:', '')
        if not g or g.upper() in {'NA', 'NAN', 'NONE'}:
            continue
        out.append(g)
    return out


def tex_escape(s):
    return str(s).replace('gene:', '').replace('_', r'\_')


def parse_best(cell):
    t = str(cell).strip()
    if t == '' or t.upper() == 'NA':
        return ('NA', None)
    bg, bp = 'NA', None
    for part in t.split(';'):
        if ':' not in part:
            continue
        g, p = part.split(':', 1)
        try:
            pv = float(p)
        except Exception:
            continue
        if bp is None or pv < bp:
            bg, bp = g, pv
    return (bg, bp)


def fmt_logp(p):
    if p is None or p <= 0:
        return 'NA'
    return f"{-math.log10(p):.2f}"


def annotate_selected_for_species(trait, species, selected_df):
    sub = selected_df[selected_df['species'] == species].copy()
    if sub.empty:
        return {}

    tmp_in = ANNOT_DIR / f'{trait}_{species}_selected.tsv'
    tmp_out = ANNOT_DIR / f'{trait}_{species}_annot_10kb.tsv'
    out_csv = ANNOT_DIR / f'{trait}_{species}_annot_10kb.csv'

    w = pd.DataFrame({
        'chr': sub['chr'].astype(str),
        'rs': sub['snp'].astype(str),
        'ps': pd.to_numeric(sub['pos'], errors='coerce'),
        'p_wald': pd.to_numeric(sub['p'], errors='coerce'),
    }).dropna(subset=['ps', 'p_wald'])
    w.to_csv(tmp_in, sep='\t', index=False)

    cmd = [
        'Rscript', str(ANNOT_SCRIPT), str(tmp_in), GFF[species], str(tmp_out), str(WINDOW_BP)
    ]
    subprocess.run(cmd, check=True)

    ann = pd.read_csv(tmp_out, sep='\t', low_memory=False)
    ann.to_csv(out_csv, index=False)

    # keep csv-only
    if tmp_in.exists():
        tmp_in.unlink()
    if tmp_out.exists():
        tmp_out.unlink()

    if 'closest_gene' not in ann.columns or 'p_wald' not in ann.columns:
        return {}

    ann['closest_gene'] = ann['closest_gene'].astype(str).str.strip()
    ann['p_wald'] = pd.to_numeric(ann['p_wald'], errors='coerce')
    ann = ann[(ann['closest_gene'] != '') & (~ann['closest_gene'].str.lower().isin(['na', 'nan', 'none']))]
    ann = ann.dropna(subset=['p_wald'])

    gm = ann.groupby('closest_gene', as_index=False)['p_wald'].max()
    return {r['closest_gene']: float(r['p_wald']) for _, r in gm.iterrows()}


def build_trait_orthogroups_from_old_projection(trait, species_gene_pmax):
    old_all = WITHPC / trait / f'{trait}_top0p5_orthogroups_all.tsv'
    odf = pd.read_csv(old_all, sep='\t', low_memory=False)

    gene_to_og = {sp: defaultdict(set) for sp in SPECIES_ORDER}
    for _, r in odf.iterrows():
        og = str(r['Orthogroup'])
        for sp in SPECIES_ORDER:
            for g in split_genes(r.get(f'{sp}_gene', 'NA')):
                gene_to_og[sp][g].add(og)

    og_species_genes = defaultdict(lambda: defaultdict(set))
    for sp in SPECIES_ORDER:
        for g in species_gene_pmax[sp].keys():
            for og in gene_to_og[sp].get(g, set()):
                og_species_genes[og][sp].add(g)

    rows = []
    for og, spm in og_species_genes.items():
        row = {'Orthogroup': og}
        gc = 0
        for sp in SPECIES_ORDER:
            genes = sorted(spm.get(sp, set()))
            if genes:
                row[f'{sp}_gene'] = ';'.join(genes)
                row[f'{sp}_pvalue_max'] = ';'.join([f"{g}:{species_gene_pmax[sp][g]:.6g}" for g in genes])
                gc += 1
            else:
                row[f'{sp}_gene'] = 'NA'
                row[f'{sp}_pvalue_max'] = 'NA'
        row['Gene_count'] = gc
        rows.append(row)

    cols = ['Orthogroup', 'Gene_count'] + [x for sp in SPECIES_ORDER for x in (f'{sp}_gene', f'{sp}_pvalue_max')]
    all_df = pd.DataFrame(rows, columns=cols)
    if not all_df.empty:
        all_df = all_df.sort_values(['Gene_count', 'Orthogroup'], ascending=[False, True])
    return all_df


def main():
    trait_counts = {}

    print('[1/4] Rebuild mixed-window (Arabidopsis/Rice 10kb) gene tables...')
    for trait, label in TRAITS.items():
        selected_file = WITHPC / trait / f'{trait}_top0p5_selected_snps.csv'
        old_gene_file = WITHPC / trait / f'{trait}_top0p5_gene_table.csv'
        count_file = WITHPC / trait / f'{trait}_top0p5_and_bonf_selection_counts.csv'

        if not selected_file.exists() or not old_gene_file.exists() or not count_file.exists():
            print(f'  [skip] {trait} missing input files')
            continue

        selected = pd.read_csv(selected_file)
        old_gene = pd.read_csv(old_gene_file)
        counts = pd.read_csv(count_file)

        sp_map = {sp: {} for sp in SPECIES_ORDER}

        # keep 25kb species from existing top0p5 gene table
        for sp in KEEP_25KB:
            s = old_gene[old_gene['species'] == sp].copy()
            if s.empty:
                continue
            s['pvalue'] = pd.to_numeric(s['pvalue'], errors='coerce')
            s = s.dropna(subset=['gene', 'pvalue'])
            gm = s.groupby('gene', as_index=False)['pvalue'].max()
            sp_map[sp] = {r['gene']: float(r['pvalue']) for _, r in gm.iterrows()}

        # re-annotate arab + rice at 10kb
        for sp in REANNOTATE:
            sp_map[sp] = annotate_selected_for_species(trait, sp, selected)

        # write gene table csv
        rows = []
        for sp in SPECIES_ORDER:
            for g, p in sp_map[sp].items():
                rows.append({'species': sp, 'gene': g, 'pvalue': p})
        gdf = pd.DataFrame(rows).sort_values(['species', 'pvalue', 'gene']) if rows else pd.DataFrame(columns=['species', 'gene', 'pvalue'])
        gdf.to_csv(OUT / f'{trait}_top0p5_gene_table_mixed_window_10kb.csv', index=False)

        top = {r['species']: int(r['n_top0p5_snps']) for _, r in counts.iterrows()}
        bonf = {r['species']: int(r['n_bonf_snps']) for _, r in counts.iterrows()}

        count_rows = []
        for sp in SPECIES_ORDER:
            count_rows.append({
                'trait': trait,
                'species': sp,
                'top0p5_snps': top.get(sp, 0),
                'bonf_snps': bonf.get(sp, 0),
                'genes_mixed_window_10kb': len(sp_map[sp]),
                'annotation_window_bp': WINDOW_BP if sp in REANNOTATE else 25000,
            })
        pd.DataFrame(count_rows).to_csv(OUT / f'{trait}_top0p5_counts_mixed_window_10kb.csv', index=False)

        trait_counts[trait] = {
            'label': label,
            'top_snps': top,
            'bonf_snps': bonf,
            'gene_counts': {sp: len(sp_map[sp]) for sp in SPECIES_ORDER},
            'species_gene_pmax': sp_map,
        }
        print(f"  [ok] {trait}: " + ', '.join([f"{sp}={len(sp_map[sp])}" for sp in SPECIES_ORDER]))

    print('[2/4] Build orthogroup tables (CSV)...')
    ortho_summary = []
    for trait in trait_counts:
        all_df = build_trait_orthogroups_from_old_projection(trait, trait_counts[trait]['species_gene_pmax'])
        ge2 = all_df[all_df['Gene_count'] >= 2].copy() if not all_df.empty else all_df.copy()
        all5 = all_df[all_df['Gene_count'] == 5].copy() if not all_df.empty else all_df.copy()

        all_df.to_csv(ORTHO_DIR / f'{trait}_top0p5_orthogroups_all_mixed_window_10kb.csv', index=False)
        ge2.to_csv(ORTHO_DIR / f'{trait}_top0p5_orthogroups_common_ge2_mixed_window_10kb.csv', index=False)
        all5.to_csv(ORTHO_DIR / f'{trait}_top0p5_orthogroups_common5_mixed_window_10kb.csv', index=False)

        trait_counts[trait]['orthogroups_total'] = len(all_df)
        trait_counts[trait]['orthogroups_common5'] = len(all5)

        ortho_summary.append({
            'trait': trait,
            'orthogroups_total': len(all_df),
            'orthogroups_common_ge2': len(ge2),
            'orthogroups_common5': len(all5),
        })
        print(f"  [ok] {trait}: total={len(all_df)}, common5={len(all5)}")

    pd.DataFrame(ortho_summary).to_csv(OUT / 'orthogroups_summary_mixed_window_10kb.csv', index=False)

    print('[3/4] Write LaTeX sections...')
    lines = ['% Auto-generated mixed-window report (Arabidopsis/Rice 10 kb; others 25 kb)', '']
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
            'To evaluate cross-species convergence beyond the Bonferroni tail, we selected the top 0.5\\% of SNPs per species and mapped nearby genes to orthogroups, using \\pm10 kb annotation for Arabidopsis and rice and \\pm25 kb annotation for maize, sorghum, and barley. '
            f"This yielded {info['gene_counts'].get('Arabidopsis',0):,} Arabidopsis genes, {info['gene_counts'].get('Rice',0):,} rice genes, {info['gene_counts'].get('Maize',0):,} maize genes, {info['gene_counts'].get('Barley',0):,} barley genes, and {info['gene_counts'].get('Sorghum',0):,} sorghum genes."
        )
        lines.append('')
        lines.append(
            f"Orthogroup integration identified {info.get('orthogroups_total',0):,} total orthogroups, with {info.get('orthogroups_common5',0):,} orthogroups represented in all five species. "
            f"Table~\\ref{{tab:shared_function_candidates_{trait}_mixed10kb}} summarizes these shared orthogroups and reports, for each species, the strongest gene-level signal within each orthogroup as $-\\log_{{10}}(P)$."
        )
        lines.append('')

        cfile = ORTHO_DIR / f'{trait}_top0p5_orthogroups_common5_mixed_window_10kb.csv'
        cdf = pd.read_csv(cfile) if cfile.exists() and cfile.stat().st_size > 0 else pd.DataFrame()

        lines.append('\\begin{table*}[t]')
        lines.append(f"\\caption{{Shared orthogroups across all five species from the top 0.5\\% of SNPs for {label} (Arabidopsis/Rice \\pm10 kb; Maize/Sorghum/Barley \\pm25 kb).\\label{{tab:shared_function_candidates_{trait}_mixed10kb}}}}")
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
                    g, p = parse_best(r.get(f'{sp}_pvalue_max', 'NA'))
                    if g == 'NA':
                        raw = str(r.get(f'{sp}_gene', 'NA'))
                        g = raw.split(';')[0] if raw and raw.upper() != 'NA' else 'NA'
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

    TEX_FILE.write_text('\n'.join(lines) + '\n', encoding='utf-8')

    print('[4/4] Done')
    print(f'Output: {OUT}')
    print(f'TeX: {TEX_FILE}')


if __name__ == '__main__':
    main()
