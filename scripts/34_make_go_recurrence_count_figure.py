from pathlib import Path
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).resolve().parents[1]
GO_DIR = BASE / 'results' / 'tables' / 'GO_terms'
FIG_DIR = BASE / 'Figs' / 'Supplementary'
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR = BASE / 'results' / 'tables' / 'GO_terms'

SPECIES = ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']
TRAITS = ['PC1', 'pH', 'soilN']
ONTOS = ['GO_BP', 'GO_MF']
ONTO_LABEL = {'GO_BP': 'GO biological process', 'GO_MF': 'GO molecular function'}
TRAIT_LABEL = {'PC1': 'PC1', 'pH': 'soil pH', 'soilN': 'soil nitrogen'}
PANEL = {('PC1','GO_BP'):'A',('PC1','GO_MF'):'B',('pH','GO_BP'):'C',('pH','GO_MF'):'D',('soilN','GO_BP'):'E',('soilN','GO_MF'):'F'}


def wrap_term(term, width=38):
    return '\n'.join(textwrap.wrap(term, width=width, break_long_words=False))


def parse_panther(path: Path, species: str, trait: str, ontology: str) -> pd.DataFrame:
    with open(path) as fh:
        lines = [line.rstrip('\n') for line in fh]
    start = None
    for i, line in enumerate(lines):
        if line.startswith('GO biological process complete') or line.startswith('GO molecular function complete'):
            start = i + 1
            break
    rows = []
    for line in lines[start:]:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 7:
            continue
        try:
            p = float(parts[6].replace('E', 'e'))
            fold = float(parts[5])
        except ValueError:
            continue
        rows.append({
            'term': parts[0],
            'direction': parts[4],
            'fold_enrichment': fold,
            'p_value': p,
            'species': species,
            'trait': trait,
            'ontology': ontology,
        })
    return pd.DataFrame(rows)


def load_all():
    dfs = []
    for species in SPECIES:
        for trait in TRAITS:
            for onto in ONTOS:
                path = GO_DIR / f'{species}_{trait}_{onto}.txt'
                dfs.append(parse_panther(path, species, trait, onto))
    df = pd.concat(dfs, ignore_index=True)
    df = df[(df['direction'] == '+') & (df['p_value'] < 0.05)].copy()
    return df


def recurrence_table(df, trait, onto):
    sub = df[(df['trait'] == trait) & (df['ontology'] == onto)].copy()
    if sub.empty:
        return pd.DataFrame()
    rec = sub.groupby('term').agg(
        n_species=('species', 'nunique'),
        min_p=('p_value', 'min'),
        max_fold=('fold_enrichment', 'max'),
        species_list=('species', lambda x: ';'.join(sorted(set(x))))
    ).reset_index()
    rec = rec[rec['n_species'] >= 2].copy()
    rec = rec.sort_values(['n_species', 'min_p'], ascending=[False, True]).head(10)
    return rec


def draw_panel(ax, rec, trait, onto, cmap):
    rec = rec.copy()
    rec['wrapped_term'] = rec['term'].map(lambda x: wrap_term(x, 42))
    rec['neglog10p'] = -np.log10(rec['min_p'])

    if rec.empty:
        ax.text(0.5, 0.5, 'No recurring terms\n(recurrence >= 2 species)', ha='center', va='center')
        ax.axis('off')
        return

    rec = rec.sort_values(['n_species', 'min_p'], ascending=[True, False])
    colors = cmap((rec['neglog10p'] - rec['neglog10p'].min()) / (rec['neglog10p'].max() - rec['neglog10p'].min() + 1e-9))
    bars = ax.barh(rec['wrapped_term'], rec['n_species'], color=colors, edgecolor='black', linewidth=0.4)
    for bar, sp in zip(bars, rec['species_list']):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2, sp, va='center', fontsize=8)

    ax.set_xlim(0, 5.6)
    ax.set_xticks(range(0, 6))
    ax.set_xlabel('Number of species with significant enrichment')
    ax.set_ylabel('Recurring GO terms')
    ax.set_title(f"{PANEL[(trait, onto)]}. {TRAIT_LABEL[trait]}: {ONTO_LABEL[onto]}", loc='left', fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=rec['neglog10p'].min(), vmax=rec['neglog10p'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.01, fraction=0.035)
    cbar.set_label('-log10(min raw P-value)', fontsize=9)


def make_trait_figure(df, trait, cmap):
    fig, axes = plt.subplots(1, 2, figsize=(18, 14), constrained_layout=True)
    for j, onto in enumerate(ONTOS):
        rec = recurrence_table(df, trait, onto)
        draw_panel(axes[j], rec, trait, onto, cmap)

    fig.suptitle(f"Recurring enriched GO terms across species: {TRAIT_LABEL[trait]}", fontsize=17, fontweight='bold')
    fig.text(0.5, 0.01, 'Only overrepresented terms with raw P < 0.05 were considered. Panels show the top 10 recurring exact GO terms present in at least two species for each ontology.', ha='center', fontsize=10)

    png = FIG_DIR / f'SuppFig_GO_recurrence_counts_{trait}_top10.png'
    pdf = FIG_DIR / f'SuppFig_GO_recurrence_counts_{trait}_top10.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    print('WROTE', png)
    print('WROTE', pdf)


def main():
    df = load_all()
    fig, axes = plt.subplots(3, 2, figsize=(15, 18), constrained_layout=True)
    cmap = sns.color_palette('YlOrRd', as_cmap=True)
    all_tables = []

    for i, trait in enumerate(TRAITS):
        for j, onto in enumerate(ONTOS):
            ax = axes[i, j]
            rec = recurrence_table(df, trait, onto)
            all_tables.append(rec.assign(trait=trait, ontology=onto))
            draw_panel(ax, rec, trait, onto, cmap)

    fig.suptitle('Recurring enriched GO terms across species', fontsize=17, fontweight='bold')
    fig.text(0.5, 0.01, 'Only overrepresented terms with raw P < 0.05 were considered. Panels show the top 10 recurring exact GO terms present in at least two species for each trait and ontology.', ha='center', fontsize=10)
    png = FIG_DIR / 'SuppFig_GO_recurrence_counts_PC1_pH_soilN.png'
    pdf = FIG_DIR / 'SuppFig_GO_recurrence_counts_PC1_pH_soilN.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    print('WROTE', png)
    print('WROTE', pdf)

    out = pd.concat(all_tables, ignore_index=True)
    out.to_csv(OUT_DIR / 'GO_recurrence_counts_PC1_pH_soilN.csv', index=False)
    print('WROTE', OUT_DIR / 'GO_recurrence_counts_PC1_pH_soilN.csv')

    for trait in TRAITS:
        make_trait_figure(df, trait, cmap)


if __name__ == '__main__':
    main()
