import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary'
OUTDIR.mkdir(parents=True, exist_ok=True)

WINDOW_BP = 500_000
SPECIES_ORDER = ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']
TRAITS = {
    'PC1': {
        'selected_snps': BASE / 'results' / 'pc1_worldclim' / 'pc1_worldclim_top0p5_selected_snps.csv',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc1_worldclim_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc1_worldclim_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Barley_PC1.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Maize_PC1.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Sorghum_PC1.annot_25000bp.tsv', '\t', 25),
        },
    },
    'PC2': {
        'selected_snps': BASE / 'results' / 'pc2_worldclim' / 'pc2_worldclim_top0p5_selected_snps.csv',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc2_worldclim_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc2_worldclim_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC2' / 'Barley_PC2.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC2' / 'Maize_PC2.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC2' / 'Sorghum_PC2.annot_25000bp.tsv', '\t', 25),
        },
    },
    'PC3': {
        'selected_snps': BASE / 'results' / 'pc3_worldclim' / 'pc3_worldclim_top0p5_selected_snps.csv',
        'sources': {
            'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc3_worldclim_Arabidopsis_annot_10kb.csv', ',', 10),
            'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc3_worldclim_Rice_annot_10kb.csv', ',', 10),
            'Barley':      (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC3' / 'Barley_PC3.annot_25000bp.tsv', '\t', 25),
            'Maize':       (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC3' / 'Maize_PC3.annot_25000bp.tsv', '\t', 25),
            'Sorghum':     (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC3' / 'Sorghum_PC3.annot_25000bp.tsv', '\t', 25),
        },
    },
}


def latex_escape(text):
    if pd.isna(text):
        return 'NA'
    text = str(text)
    repl = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}'}
    for k, v in repl.items():
        text = text.replace(k, v)
    return text


def load_annotations(config):
    anns = {}
    for sp, (path, sep, window_kb) in config['sources'].items():
        df = pd.read_csv(path, sep=sep)
        chr_col = 'chr' if 'chr' in df.columns else 'Chr'
        pos_col = 'ps' if 'ps' in df.columns else 'Pos'
        df['species'] = sp
        df['chr_str'] = df[chr_col].astype(str).str.replace(r'^chr', '', regex=True)
        df['pos_int'] = pd.to_numeric(df[pos_col], errors='coerce').astype('Int64')
        keep = df[['species', 'chr_str', 'pos_int', 'closest_gene_id', 'genes_in_25kb']].copy()
        keep['window_kb'] = window_kb
        anns[sp] = keep
    return anns


def top3_gene_loci(phenotype, config):
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
        top['Phenotype'] = phenotype
        top['Locus_rank'] = range(1, len(top) + 1)
        top['Lead_gene'] = top['closest_gene_id']
        top['Candidate_genes_in_window'] = top['genes_in_25kb'].fillna('NA')
        rows.append(top[['Phenotype', 'species', 'Locus_rank', 'chr_str', 'pos_int', 'snp', 'logp', 'Lead_gene', 'Candidate_genes_in_window']])
    out = pd.concat(rows, ignore_index=True)
    out.columns = ['Phenotype', 'Species', 'Locus_rank', 'Chr', 'Lead_SNP_Position', 'Lead_SNP_ID', 'NegLog10P', 'Lead_gene', 'Candidate_genes_in_window']
    return out


def main():
    frames = [top3_gene_loci(ph, cfg) for ph, cfg in TRAITS.items()]
    out = pd.concat(frames, ignore_index=True)
    csv_path = OUTDIR / 'SuppTable10_WorldClim_PC1_PC2_PC3_top3_loci_per_species.csv'
    tex_path = OUTDIR / 'SuppTable10_WorldClim_PC1_PC2_PC3_top3_loci_per_species.tex'
    out.to_csv(csv_path, index=False)

    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\caption{Top three non-overlapping lead loci per species for WorldClim PC1, PC2, and PC3. Loci were defined by ranking top 0.5\% SNPs within each species and principal component and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\pm$10 kb; maize, barley, and sorghum: $\pm$25 kb).}')
    lines.append(r'\label{tab:supptable10_worldclim_pc123_top3_loci}')
    lines.append(r'\scriptsize')
    lines.append(r'\tabcolsep=3pt')
    lines.append(r'\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllcclp{3.4cm}p{4.8cm}@{}}')
    lines.append(r'\toprule')
    lines.append(r'Phenotype & Species & Rank & Chr & Lead SNP position & $-\log_{10}(P)$ & Lead gene & Candidate genes in annotation window \\')
    lines.append(r'\midrule')
    for _, r in out.iterrows():
        line = (
            f"{latex_escape(r['Phenotype'])} & {latex_escape(r['Species'])} & {int(r['Locus_rank'])} & "
            f"{latex_escape(r['Chr'])} & {int(r['Lead_SNP_Position'])} & {r['NegLog10P']:.2f} & "
            f"{latex_escape(r['Lead_gene'])} & {latex_escape(r['Candidate_genes_in_window'])} \\\\" 
        )
        lines.append(line)
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular*}')
    lines.append(r'\end{table*}')
    tex_path.write_text('\n'.join(lines) + '\n')
    print('WROTE', csv_path)
    print('WROTE', tex_path)
    print(out.head(15).to_string(index=False))


if __name__ == '__main__':
    main()
