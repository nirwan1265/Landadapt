import math
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).resolve().parents[1]
FIGDIR = BASE / 'Figs' / 'FigX'
FIGDIR.mkdir(parents=True, exist_ok=True)
OUTCSV = BASE / 'results' / 'tables' / 'Supplementary' / 'FigX_cross_species_summary_values.csv'

species_order = ['Arabidopsis', 'Barley', 'Rice', 'Maize', 'Sorghum']
trait_order = ['PC1', 'PC2', 'PC3', 'pH', 'soilN', 'AM_rel', 'AM_roots', 'aridity']
trait_labels = {
    'PC1': 'PC1', 'PC2': 'PC2', 'PC3': 'PC3', 'pH': 'pH', 'soilN': 'SoilN',
    'AM_rel': 'AM rel', 'AM_roots': 'AM roots', 'aridity': 'Aridity'
}

# Panel A: Bonferroni-significant SNP counts from the manuscript/current raw GWAS summaries.
bonf = {
    'PC1': {'Arabidopsis': 38, 'Barley': 8, 'Rice': 602, 'Maize': 1, 'Sorghum': 6},
    'PC2': {'Arabidopsis': 13, 'Barley': 1, 'Rice': 0, 'Maize': 1, 'Sorghum': 0},
    'PC3': {'Arabidopsis': 29, 'Barley': 24, 'Rice': 10, 'Maize': 0, 'Sorghum': 80},
    'pH': {'Arabidopsis': 6, 'Barley': 6, 'Rice': 383, 'Maize': 0, 'Sorghum': 3},
    'soilN': {'Arabidopsis': 64, 'Barley': 4, 'Rice': 202, 'Maize': 1, 'Sorghum': 2},
    'AM_rel': {'Arabidopsis': 824, 'Barley': 9, 'Rice': 268, 'Maize': 1, 'Sorghum': 0},
    'AM_roots': {'Arabidopsis': 1160, 'Barley': 125, 'Rice': 1536, 'Maize': 5, 'Sorghum': 4},
    'aridity': {'Arabidopsis': 436, 'Barley': 1, 'Rice': 0, 'Maize': 4, 'Sorghum': 90},
}

# Panel B: top-0.5% mapped gene counts (current final mixed-window setup)
gene_table_files = {
    'PC1': BASE / 'results' / 'mixed_window_10kb_csv' / 'pc1_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    'PC2': BASE / 'results' / 'mixed_window_10kb_csv' / 'pc2_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    'PC3': BASE / 'results' / 'mixed_window_10kb_csv' / 'pc3_worldclim_top0p5_gene_table_mixed_window_10kb.csv',
    'pH': BASE / 'results' / 'mixed_window_10kb_csv' / 'ph_top0p5_gene_table_mixed_window_10kb.csv',
    'soilN': BASE / 'results' / 'mixed_window_10kb_csv' / 'soilN_top0p5_gene_table_mixed_window_10kb.csv',
    'AM_rel': BASE / 'results' / 'mixed_window_10kb_csv' / 'am_rel_abundance_colonization_top0p5_gene_table_mixed_window_10kb.csv',
    'AM_roots': BASE / 'results' / 'mixed_window_10kb_csv' / 'am_roots_colonized_top0p5_gene_table_mixed_window_10kb.csv',
    'aridity': BASE / 'results' / 'mixed_window_10kb_csv' / 'aridity_index_top0p5_gene_table_mixed_window_10kb.csv',
}

gene_counts = {}
for trait, path in gene_table_files.items():
    df = pd.read_csv(path)
    gene_counts[trait] = df.groupby('species')['gene'].nunique().to_dict()

# Panels C/D: shared orthogroup counts from current supplementary tables.
pc = pd.read_csv(BASE / 'results' / 'tables' / 'Supplementary' / 'PC' / 'SuppTable3A_PC_orthogroups_atleast_2_desc.csv')
no = pd.read_csv(BASE / 'results' / 'tables' / 'Supplementary' / 'SuppTable3A_orthogroups_atleast_2_desc_noPC.csv')
no = no[~no['Phenotype'].isin(['PC1','PC2','PC3'])].copy()
all_orth = pd.concat([pc, no], ignore_index=True)
all_orth['TraitSimple'] = all_orth['Phenotype'].replace({
    'am_rel_abundance_colonization': 'AM_rel',
    'am_roots_colonized': 'AM_roots',
    'aridity_index': 'aridity'
})
all5_counts = all_orth[all_orth['Gene_count'] == 5].groupby('TraitSimple').size().to_dict()
ge4_counts = all_orth[all_orth['Gene_count'] >= 4].groupby('TraitSimple').size().to_dict()

# Build matrices
A = pd.DataFrame(index=species_order, columns=trait_order, data=0)
B = pd.DataFrame(index=species_order, columns=trait_order, data=0)
for trait in trait_order:
    for sp in species_order:
        A.loc[sp, trait] = bonf[trait].get(sp, 0)
        B.loc[sp, trait] = gene_counts[trait].get(sp, 0)
A = A.astype(int)
B = B.astype(int)
C = pd.Series({t: int(all5_counts.get(t, 0)) for t in trait_order})
D = pd.Series({t: int(ge4_counts.get(t, 0)) for t in trait_order})

# Save source values
rows = []
for sp in species_order:
    for tr in trait_order:
        rows.append({'panel': 'A', 'species': sp, 'trait': tr, 'value': int(A.loc[sp, tr])})
        rows.append({'panel': 'B', 'species': sp, 'trait': tr, 'value': int(B.loc[sp, tr])})
for tr in trait_order:
    rows.append({'panel': 'C', 'species': '', 'trait': tr, 'value': int(C.loc[tr])})
    rows.append({'panel': 'D', 'species': '', 'trait': tr, 'value': int(D.loc[tr])})
pd.DataFrame(rows).to_csv(OUTCSV, index=False)

sns.set_theme(style='whitegrid')
fig = plt.figure(figsize=(15, 11))
gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1], wspace=0.25, hspace=0.28)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[0, 1])
axD = fig.add_subplot(gs[1, 1])

A_plot = np.log10(A + 1)
B_plot = np.log10(B + 1)
annotA = A.applymap(lambda x: f'{x:,}')
annotB = B.applymap(lambda x: f'{x:,}')


def add_heatmap_labels(ax, values, labels, fontsize=10):
    """Add explicit labels so every exported cell remains readable."""
    arr = values.to_numpy(dtype=float)
    label_arr = labels.to_numpy()
    threshold = arr.min() + 0.55 * (arr.max() - arr.min())
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            color = 'white' if arr[i, j] >= threshold else 'black'
            ax.text(
                j + 0.5,
                i + 0.5,
                label_arr[i, j],
                ha='center',
                va='center',
                fontsize=fontsize,
                color=color,
            )


sns.heatmap(A_plot, ax=axA, cmap='YlOrRd', annot=False, cbar_kws={'label': 'log10(count + 1)'}, linewidths=0.5, linecolor='white')
add_heatmap_labels(axA, A_plot, annotA, fontsize=10)
axA.set_title('A. Bonferroni-significant SNP counts', loc='left', fontweight='bold')
axA.set_xlabel('Trait')
axA.set_ylabel('Species')
axA.set_xticklabels([trait_labels[t] for t in trait_order], rotation=45, ha='right')

sns.heatmap(B_plot, ax=axB, cmap='YlGnBu', annot=False, cbar_kws={'label': 'log10(count + 1)'}, linewidths=0.5, linecolor='white')
add_heatmap_labels(axB, B_plot, annotB, fontsize=10)
axB.set_title('B. Top-0.5% mapped gene counts', loc='left', fontweight='bold')
axB.set_xlabel('Trait')
axB.set_ylabel('Species')
axB.set_xticklabels([trait_labels[t] for t in trait_order], rotation=45, ha='right')

x = np.arange(len(trait_order))
bar_colors = sns.color_palette('Set2', n_colors=len(trait_order))
axC.bar(x, [C[t] for t in trait_order], color=bar_colors, edgecolor='black', linewidth=0.5)
axC.set_title('C. All-five shared orthogroups', loc='left', fontweight='bold')
axC.set_xticks(x)
axC.set_xticklabels([trait_labels[t] for t in trait_order], rotation=45, ha='right')
axC.set_ylabel('Orthogroup count')
for i, t in enumerate(trait_order):
    axC.text(i, C[t] + 0.15, str(int(C[t])), ha='center', va='bottom', fontsize=10, color='black')
axC.spines[['top', 'right']].set_visible(False)

axD.bar(x, [D[t] for t in trait_order], color=bar_colors, edgecolor='black', linewidth=0.5)
axD.set_title('D. Shared orthogroups in >=4 species', loc='left', fontweight='bold')
axD.set_xticks(x)
axD.set_xticklabels([trait_labels[t] for t in trait_order], rotation=45, ha='right')
axD.set_ylabel('Orthogroup count')
for i, t in enumerate(trait_order):
    axD.text(i, D[t] + 0.25, str(int(D[t])), ha='center', va='bottom', fontsize=10, color='black')
axD.spines[['top', 'right']].set_visible(False)

fig.suptitle('Cross-species environmental GWAS summary across traits', fontsize=16, fontweight='bold', y=0.98)
fig.text(0.5, 0.005, 'AM rel = AM fungal relative abundance colonization; AM roots = AM fungal roots colonized', ha='center', fontsize=10)
fig.tight_layout(rect=[0, 0.02, 1, 0.96])

png = FIGDIR / 'FigX_cross_species_gwas_summary.png'
pdf = FIGDIR / 'FigX_cross_species_gwas_summary.pdf'
fig.savefig(png, dpi=300, bbox_inches='tight')
fig.savefig(pdf, bbox_inches='tight')
print('WROTE', png)
print('WROTE', pdf)
print('WROTE', OUTCSV)
