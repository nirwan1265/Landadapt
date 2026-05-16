import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary'
OUTDIR.mkdir(parents=True, exist_ok=True)

WINDOW_BP = 500_000
SPECIES_ORDER = ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']

TRAITS = {
    'ph': {
        'selected_snps': BASE / 'results' / 'ph' / 'ph_top0p5_selected_snps.csv',
        'label': 'soil pH',
        'slug': 'pH',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'ph_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'ph_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'ph_Barley_annot_25kb.csv', ',', 25),
            'Maize':       (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'ph_Maize_annot_25kb.csv', ',', 25),
            'Sorghum':     (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'ph_Sorghum_annot_25kb.csv', ',', 25),
        },
    },
    'soilN': {
        'selected_snps': BASE / 'results' / 'soilN' / 'soilN_top0p5_selected_snps.csv',
        'label': 'soil nitrogen (0--5 cm)',
        'slug': 'soilN',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'soilN_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'soilN_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'soilN' / 'Barley_soilN.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'soilN' / 'Maize_soilN.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'soilN' / 'Sorghum_soilN.annot_25000bp.tsv', '\t', 25),
        },
    },
    'aridity_index': {
        'selected_snps': BASE / 'results' / 'aridity_index' / 'aridity_index_top0p5_selected_snps.csv',
        'label': 'aridity index',
        'slug': 'aridity_index',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'aridity_index_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'aridity_index_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'aridity_index_Barley_annot_25kb.csv', ',', 25),
            'Maize':       (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'aridity_index_Maize_annot_25kb.csv', ',', 25),
            'Sorghum':     (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'aridity_index_Sorghum_annot_25kb.csv', ',', 25),
        },
    },
    'am_rel_abundance_colonization': {
        'selected_snps': BASE / 'results' / 'am_rel_abundance_colonization' / 'am_rel_abundance_colonization_top0p5_selected_snps.csv',
        'label': 'AM fungal relative abundance colonization',
        'slug': 'am_rel_abundance_colonization',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'am_rel_abundance_colonization_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'am_rel_abundance_colonization_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_rel_abundance_colonization' / 'Barley_am_rel_abundance_colonization.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_rel_abundance_colonization' / 'Maize_am_rel_abundance_colonization.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_rel_abundance_colonization' / 'Sorghum_am_rel_abundance_colonization.annot_25000bp.tsv', '\t', 25),
        },
    },
    'am_roots_colonized': {
        'selected_snps': BASE / 'results' / 'am_roots_colonized' / 'am_roots_colonized_top0p5_selected_snps.csv',
        'label': 'AM fungal roots colonized',
        'slug': 'am_roots_colonized',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'am_roots_colonized_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'am_roots_colonized_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_roots_colonized' / 'Barley_am_roots_colonized.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_roots_colonized' / 'Maize_am_roots_colonized.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits_am' / 'annotated_25000bp' / 'am_roots_colonized' / 'Sorghum_am_roots_colonized.annot_25000bp.tsv', '\t', 25),
        },
    },
}


def latex_escape(text):
    if pd.isna(text):
        return 'NA'
    text = str(text)
    repl = {
        '\\': r'\textbackslash{}',
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
    }
    for k, v in repl.items():
        text = text.replace(k, v)
    return text


def load_annotations(config):
    anns = {}
    for sp, (path, sep, window_kb) in config['sources'].items():
        df = pd.read_csv(path, sep=sep)
        chr_col = 'chr' if 'chr' in df.columns else 'Chr'
        pos_col = 'ps' if 'ps' in df.columns else 'Pos'
        genes_col = 'genes_in_25kb'
        cg_col = 'closest_gene_id'
        dist_col = 'closest_gene_distance_bp'
        df['species'] = sp
        df['chr_str'] = df[chr_col].astype(str).str.replace(r'^chr', '', regex=True)
        df['pos_int'] = pd.to_numeric(df[pos_col], errors='coerce').astype('Int64')
        anns[sp] = df[['species', 'chr_str', 'pos_int', cg_col, genes_col, dist_col]].copy()
        anns[sp].columns = ['species', 'chr_str', 'pos_int', 'closest_gene_id', 'genes_in_25kb', 'closest_gene_distance_bp']
        anns[sp]['window_kb'] = window_kb
    return anns


def top3_gene_loci(config):
    sel = pd.read_csv(config['selected_snps'])
    anns = load_annotations(config)
    rows = []
    for sp in SPECIES_ORDER:
        g = sel[sel['species'] == sp].copy()
        g['chr_str'] = g['chr'].astype(str).str.replace(r'^chr', '', regex=True)
        g['pos_int'] = pd.to_numeric(g['pos'], errors='coerce')
        g = g.sort_values(['logp', 'p'], ascending=[False, True])
        chosen = []
        for _, r in g.iterrows():
            if any(r['chr_str'] == c['chr_str'] and abs(r['pos_int'] - c['pos_int']) <= WINDOW_BP for c in chosen):
                continue
            ann_match = anns[sp][(anns[sp]['chr_str'] == r['chr_str']) & (anns[sp]['pos_int'] == r['pos_int'])]
            has_gene = False
            if not ann_match.empty:
                cg = ann_match.iloc[0]['closest_gene_id']
                has_gene = pd.notna(cg) and str(cg) != 'NA' and str(cg).strip() != ''
            if not has_gene:
                continue
            chosen.append(r)
            if len(chosen) == 3:
                break
        top = pd.DataFrame(chosen)
        top = top.merge(anns[sp], on=['species', 'chr_str', 'pos_int'], how='left')
        top['locus_rank'] = range(1, len(top) + 1)
        top['window_bp'] = top['window_kb'].fillna(25).astype(int) * 1000
        top['lead_gene_or_status'] = top['closest_gene_id'].fillna('Intergenic (no gene within annotation window)')
        top['candidate_genes'] = top['genes_in_25kb'].fillna('NA')
        top['Trait'] = config['label']
        rows.append(top)
    out = pd.concat(rows, ignore_index=True)
    out = out[['Trait', 'species', 'locus_rank', 'chr_str', 'pos_int', 'snp', 'logp', 'lead_gene_or_status', 'candidate_genes', 'window_bp']]
    out.columns = ['Trait', 'Species', 'Locus_rank', 'Chr', 'Lead_SNP_Position', 'Lead_SNP_ID', 'NegLog10P', 'Lead_gene_or_status', 'Candidate_genes_in_window', 'Annotation_window_bp']
    return out


def write_single_table(df, title, label, stem):
    csv_path = OUTDIR / f'{stem}.csv'
    tex_path = OUTDIR / f'{stem}.tex'
    df.to_csv(csv_path, index=False)
    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(rf'\caption{{{title}}}')
    lines.append(rf'\label{{{label}}}')
    lines.append(r'\scriptsize')
    lines.append(r'\tabcolsep=3pt')
    lines.append(r'\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllclp{4.2cm}p{5.4cm}@{}}')
    lines.append(r'\toprule')
    lines.append(r'Species & Rank & Chr & Lead SNP position & $-\log_{10}(P)$ & Lead gene & Candidate genes in annotation window \\')
    lines.append(r'\midrule')
    for _, r in df.iterrows():
        lines.append(
            f"{latex_escape(r['Species'])} & {int(r['Locus_rank'])} & {latex_escape(r['Chr'])} & {int(r['Lead_SNP_Position'])} & {r['NegLog10P']:.2f} & {latex_escape(r['Lead_gene_or_status'])} & {latex_escape(r['Candidate_genes_in_window'])} \\\\"
        )
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular*}')
    lines.append(r'\end{table*}')
    tex_path.write_text('\n'.join(lines) + '\n')


def write_amf_combined_table(df, stem='SuppTable9_AMF_top3_loci_per_species'):
    csv_path = OUTDIR / f'{stem}.csv'
    tex_path = OUTDIR / f'{stem}.tex'
    df.to_csv(csv_path, index=False)
    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\caption{Top three non-overlapping lead loci per species for the two AM fungal phenotypes. Loci were defined by ranking top 0.5\% SNPs within each species and phenotype and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\pm$10 kb; maize, barley, and sorghum: $\pm$25 kb).}')
    lines.append(r'\label{tab:supptable9_amf_top3_loci}')
    lines.append(r'\scriptsize')
    lines.append(r'\tabcolsep=3pt')
    lines.append(r'\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllcclp{3.6cm}p{4.8cm}@{}}')
    lines.append(r'\toprule')
    lines.append(r'Phenotype & Species & Rank & Chr & Lead SNP position & $-\log_{10}(P)$ & Lead gene & Candidate genes in annotation window \\')
    lines.append(r'\midrule')
    for _, r in df.iterrows():
        lines.append(
            f"{latex_escape(r['Trait'])} & {latex_escape(r['Species'])} & {int(r['Locus_rank'])} & {latex_escape(r['Chr'])} & {int(r['Lead_SNP_Position'])} & {r['NegLog10P']:.2f} & {latex_escape(r['Lead_gene_or_status'])} & {latex_escape(r['Candidate_genes_in_window'])} \\\\"
        )
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular*}')
    lines.append(r'\end{table*}')
    tex_path.write_text('\n'.join(lines) + '\n')


def main():
    ph = top3_gene_loci(TRAITS['ph'])
    soiln = top3_gene_loci(TRAITS['soilN'])
    arid = top3_gene_loci(TRAITS['aridity_index'])
    amrel = top3_gene_loci(TRAITS['am_rel_abundance_colonization'])
    amroot = top3_gene_loci(TRAITS['am_roots_colonized'])

    write_single_table(
        ph,
        'Top three non-overlapping soil pH lead loci per species. Loci were defined by ranking top 0.5\\% SNPs within each species and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\\pm$10 kb; maize, barley, and sorghum: $\\pm$25 kb).',
        'tab:supptable6_ph_top3_loci',
        'SuppTable6_pH_top3_loci_per_species',
    )
    write_single_table(
        soiln,
        'Top three non-overlapping soil nitrogen (0--5 cm) lead loci per species. Loci were defined by ranking top 0.5\\% SNPs within each species and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\\pm$10 kb; maize, barley, and sorghum: $\\pm$25 kb).',
        'tab:supptable7_soiln_top3_loci',
        'SuppTable7_soilN_top3_loci_per_species',
    )
    write_single_table(
        arid,
        'Top three non-overlapping aridity index lead loci per species. Loci were defined by ranking top 0.5\\% SNPs within each species and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\\pm$10 kb; maize, barley, and sorghum: $\\pm$25 kb).',
        'tab:supptable8_aridity_top3_loci',
        'SuppTable8_aridity_top3_loci_per_species',
    )
    amf = pd.concat([amrel, amroot], ignore_index=True)
    write_amf_combined_table(amf)
    print('Wrote:', 'SuppTable6_pH_top3_loci_per_species', 'SuppTable7_soilN_top3_loci_per_species', 'SuppTable8_aridity_top3_loci_per_species', 'SuppTable9_AMF_top3_loci_per_species')


if __name__ == '__main__':
    main()
