from pathlib import Path
import re
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
ONTO = ['GO_BP', 'GO_MF']
PANEL_LABEL = {'GO_BP': 'A', 'GO_MF': 'B'}
ONTO_LABEL = {'GO_BP': 'GO biological process', 'GO_MF': 'GO molecular function'}
TRAIT_LABEL = {'PC1': 'PC1', 'pH': 'soil pH', 'soilN': 'soil nitrogen'}


def wrap_term(term, width=40):
    return '\n'.join(textwrap.wrap(term, width=width, break_long_words=False))


def parse_panther(path: Path, species: str, trait: str, ontology: str) -> pd.DataFrame:
    with open(path) as fh:
        lines = [line.rstrip('\n') for line in fh]
    start = None
    for i, line in enumerate(lines):
        if line.startswith('GO biological process complete') or line.startswith('GO molecular function complete'):
            start = i + 1
            break
    if start is None:
        return pd.DataFrame(columns=['term', 'ref_count', 'test_count', 'expected', 'direction', 'fold_enrichment', 'p_value', 'species', 'trait', 'ontology'])

    rows = []
    for line in lines[start:]:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 7:
            continue
        try:
            rows.append({
                'term': parts[0],
                'ref_count': float(parts[1]),
                'test_count': float(parts[2]),
                'expected': float(str(parts[3]).replace('.', '0.', 1) if str(parts[3]).startswith('.') else parts[3]),
                'direction': parts[4],
                'fold_enrichment': float(parts[5]),
                'p_value': float(parts[6].replace('E', 'e')),
                'species': species,
                'trait': trait,
                'ontology': ontology,
            })
        except ValueError:
            continue
    return pd.DataFrame(rows)


def load_all():
    dfs = []
    for species in SPECIES:
        for trait in TRAITS:
            for ontology in ONTO:
                path = GO_DIR / f'{species}_{trait}_{ontology}.txt'
                if path.exists():
                    dfs.append(parse_panther(path, species, trait, ontology))
    return pd.concat(dfs, ignore_index=True)


def make_trait_figure(df: pd.DataFrame, trait: str):
    trait_df = df[df['trait'] == trait].copy()
    if trait_df.empty:
        return

    fig, axes = plt.subplots(2, 1, figsize=(13, 16), constrained_layout=True)
    for ax, ontology in zip(axes, ONTO):
        sub = trait_df[trait_df['ontology'] == ontology].copy()
        top_terms = (
            sub.groupby('term', as_index=False)['p_value']
               .min()
               .sort_values('p_value', ascending=True)
               .head(20)['term']
               .tolist()
        )
        plot_df = sub[sub['term'].isin(top_terms)].copy()
        plot_df['neglog10p'] = -np.log10(plot_df['p_value'])

        # Keep ordering by strongest term overall.
        term_order = (
            plot_df.groupby('term', as_index=False)['p_value']
                   .min()
                   .sort_values('p_value', ascending=False)['term']
                   .tolist()
        )
        plot_df['term_wrapped'] = plot_df['term'].map(lambda x: wrap_term(x, width=42))
        wrapped_order = [wrap_term(t, width=42) for t in term_order]
        plot_df['term_wrapped'] = pd.Categorical(plot_df['term_wrapped'], categories=wrapped_order, ordered=True)
        plot_df['species'] = pd.Categorical(plot_df['species'], categories=SPECIES, ordered=True)

        sns.scatterplot(
            data=plot_df,
            x='species',
            y='term_wrapped',
            size='fold_enrichment',
            sizes=(40, 280),
            hue='neglog10p',
            palette='YlOrRd',
            edgecolor='black',
            linewidth=0.3,
            ax=ax,
            legend='brief',
        )
        ax.set_title(f"{PANEL_LABEL[ontology]}. {ONTO_LABEL[ontology]}", loc='left', fontweight='bold')
        ax.set_xlabel('Species')
        ax.set_ylabel('Top 20 enriched terms')
        ax.grid(True, axis='x', linewidth=0.3)
        ax.grid(True, axis='y', linewidth=0.2, alpha=0.4)
        ax.spines[['top', 'right']].set_visible(False)

        # Clean legend labels a bit.
        leg = ax.get_legend()
        if leg is not None:
            leg.set_title('')
            texts = leg.get_texts()
            for t in texts:
                if t.get_text() == 'neglog10p':
                    t.set_text('-log10(P)')
                elif t.get_text() == 'fold_enrichment':
                    t.set_text('Fold enrichment')

        # Save plotted table for reproducibility.
        plot_out = OUT_DIR / f'{trait}_{ontology}_top20_terms_for_plot.csv'
        plot_df[['term', 'species', 'fold_enrichment', 'p_value', 'neglog10p']].sort_values(['term', 'species']).to_csv(plot_out, index=False)

    fig.suptitle(f"Top enriched GO terms across species for {TRAIT_LABEL[trait]} candidate genes", fontsize=16, fontweight='bold')
    fig.text(0.5, 0.01, 'Point size shows fold enrichment; color shows -log10(raw P-value). Terms are selected by the 20 smallest raw P-values across species within each ontology.', ha='center', fontsize=10)
    png = FIG_DIR / f'SuppFig_GO_{trait}_top20_dotplot.png'
    pdf = FIG_DIR / f'SuppFig_GO_{trait}_top20_dotplot.pdf'
    fig.savefig(png, dpi=300, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    plt.close(fig)
    print('WROTE', png)
    print('WROTE', pdf)


def main():
    df = load_all()
    for trait in TRAITS:
        make_trait_figure(df, trait)


if __name__ == '__main__':
    main()
