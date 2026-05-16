import pandas as pd
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUTDIR = BASE / 'results' / 'tables' / 'Supplementary'
OUTDIR.mkdir(parents=True, exist_ok=True)

WINDOW_BP = 500_000
SPECIES_ORDER = ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']
ANNOTATION_SOURCES = {
    'Arabidopsis': (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc1_worldclim_Arabidopsis_annot_10kb.csv', ',', 10),
    'Rice':        (BASE / 'results' / 'mixed_window_10kb_csv' / 'annotated_10kb_csv' / 'pc1_worldclim_Rice_annot_10kb.csv', ',', 10),
    'Barley':      (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Barley_PC1.annot_25000bp.tsv', '\t', 25),
    'Maize':       (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Maize_PC1.annot_25000bp.tsv', '\t', 25),
    'Sorghum':     (BASE / 'results' / 'top0.5pct_traits' / 'annotated_25000bp' / 'PC1' / 'Sorghum_PC1.annot_25000bp.tsv', '\t', 25),
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


def build_table():
    sel = pd.read_csv(BASE / 'results' / 'pc1_worldclim' / 'pc1_worldclim_top0p5_selected_snps.csv')
    anns = {}
    for sp, (path, sep, window_kb) in ANNOTATION_SOURCES.items():
        df = pd.read_csv(path, sep=sep)
        df['species'] = sp
        df['chr_str'] = df['chr'].astype(str).str.replace(r'^chr', '', regex=True)
        df['pos_int'] = pd.to_numeric(df['ps'], errors='coerce').astype('Int64')
        anns[sp] = df[['species', 'chr_str', 'pos_int', 'closest_gene_id', 'genes_in_25kb', 'closest_gene_distance_bp']].copy()
        anns[sp]['window_kb'] = window_kb

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
        rows.append(top)

    out = pd.concat(rows, ignore_index=True)
    out = out[['species', 'locus_rank', 'chr_str', 'pos_int', 'snp', 'logp', 'lead_gene_or_status', 'candidate_genes', 'window_bp']]
    out.columns = ['Species', 'Locus_rank', 'Chr', 'Lead_SNP_Position', 'Lead_SNP_ID', 'NegLog10P', 'Lead_gene_or_status', 'Candidate_genes_in_window', 'Annotation_window_bp']
    return out


def write_outputs(out):
    out.to_csv(OUTDIR / 'SuppTable5_PC1_top3_loci_per_species.csv', index=False)

    lines = []
    lines.append(r'\begin{table*}[t]')
    lines.append(r'\caption{Top three non-overlapping WorldClim PC1 lead loci per species. Loci were defined by ranking top 0.5\% SNPs within each species and retaining the top three lead SNPs separated by at least 500 kb on the same chromosome. Candidate genes were assigned using the same species-specific annotation windows used in the manuscript analyses (Arabidopsis and rice: $\pm$10 kb; maize, barley, and sorghum: $\pm$25 kb).}')
    lines.append(r'\label{tab:supptable5_pc1_top3_loci}')
    lines.append(r'\scriptsize')
    lines.append(r'\tabcolsep=3pt')
    lines.append(r'\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lllclp{4.2cm}p{5.4cm}@{}}')
    lines.append(r'\toprule')
    lines.append(r'Species & Rank & Chr & Lead SNP position & $-\log_{10}(P)$ & Lead gene / status & Candidate genes in annotation window \\')
    lines.append(r'\midrule')
    for _, r in out.iterrows():
        line = (
            f"{latex_escape(r['Species'])} & {int(r['Locus_rank'])} & {latex_escape(r['Chr'])} & "
            f"{int(r['Lead_SNP_Position'])} & {r['NegLog10P']:.2f} & "
            f"{latex_escape(r['Lead_gene_or_status'])} & {latex_escape(r['Candidate_genes_in_window'])} \\\\" 
        )
        lines.append(line)
    lines.append(r'\bottomrule')
    lines.append(r'\end{tabular*}')
    lines.append(r'\end{table*}')
    (OUTDIR / 'SuppTable5_PC1_top3_loci_per_species.tex').write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    out = build_table()
    write_outputs(out)
    print(out.to_string(index=False))
