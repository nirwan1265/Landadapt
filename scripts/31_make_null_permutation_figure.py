from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).resolve().parents[1]
SUM = BASE / 'results' / 'tables' / 'Null_model' / 'orthogroup_null_permutation_summary.csv'
DIST = BASE / 'results' / 'tables' / 'Null_model' / 'orthogroup_null_permutation_distributions.csv'
OUTDIR = BASE / 'Figs' / 'Supplementary'
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTCSV = BASE / 'results' / 'tables' / 'Null_model' / 'orthogroup_null_permutation_plot_values.csv'

order = ['PC1', 'PC2', 'PC3', 'pH', 'soilN', 'AM_rel', 'AM_roots', 'aridity']
labels = {
    'PC1': 'PC1', 'PC2': 'PC2', 'PC3': 'PC3', 'pH': 'pH', 'soilN': 'SoilN',
    'AM_rel': 'AM rel', 'AM_roots': 'AM roots', 'aridity': 'Aridity'
}

summary = pd.read_csv(SUM)
dist = pd.read_csv(DIST)
summary = summary.set_index('Trait_label').loc[order].reset_index()
dist['Trait_label'] = pd.Categorical(dist['Trait_label'], categories=order, ordered=True)
dist = dist.sort_values(['Trait_label', 'perm_index'])

plot_rows = []
for _, r in summary.iterrows():
    plot_rows.append({
        'Trait_label': r['Trait_label'],
        'Observed_all5': r['Observed_all5'],
        'Null_mean_all5': r['Null_mean_all5'],
        'Null_q025_all5': r['Null_q025_all5'],
        'Null_q975_all5': r['Null_q975_all5'],
        'Empirical_p_all5': r['Empirical_p_all5'],
        'Observed_ge4': r['Observed_ge4'],
        'Null_mean_ge4': r['Null_mean_ge4'],
        'Null_q025_ge4': r['Null_q025_ge4'],
        'Null_q975_ge4': r['Null_q975_ge4'],
        'Empirical_p_ge4': r['Empirical_p_ge4'],
    })
pd.DataFrame(plot_rows).to_csv(OUTCSV, index=False)

sns.set_theme(style='whitegrid')
fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True, height_ratios=[1, 1.25])
palette = sns.color_palette('Set2', n_colors=len(order))

panels = [
    ('all5_count', 'Observed_all5', 'Null_mean_all5', 'Null_q025_all5', 'Null_q975_all5', 'Empirical_p_all5', 'A. All-five shared orthogroup recurrence vs null'),
    ('ge4_count', 'Observed_ge4', 'Null_mean_ge4', 'Null_q025_ge4', 'Null_q975_ge4', 'Empirical_p_ge4', 'B. Shared orthogroups in >=4 species vs null'),
]

for ax, (dist_col, obs_col, mean_col, low_col, high_col, p_col, title) in zip(axes, panels):
    sns.violinplot(
        data=dist,
        x='Trait_label',
        y=dist_col,
        order=order,
        inner=None,
        cut=0,
        linewidth=0.8,
        palette=palette,
        ax=ax,
    )
    x = np.arange(len(order))
    means = summary[mean_col].to_numpy(dtype=float)
    lows = summary[low_col].to_numpy(dtype=float)
    highs = summary[high_col].to_numpy(dtype=float)
    obs = summary[obs_col].to_numpy(dtype=float)
    pvals = summary[p_col].to_numpy(dtype=float)
    lower_err = np.maximum(0, means - lows)
    upper_err = np.maximum(0, highs - means)

    ax.errorbar(
        x,
        means,
        yerr=[lower_err, upper_err],
        fmt='o',
        color='black',
        ecolor='black',
        elinewidth=1.2,
        capsize=4,
        markersize=5,
        zorder=4,
        label='Null mean with 95% interval',
    )
    ax.scatter(x, obs, marker='D', s=55, color='crimson', edgecolor='black', linewidth=0.4, zorder=5, label='Observed')

    for i, (o, p) in enumerate(zip(obs, pvals)):
        y = max(o, highs[i]) + (0.18 if dist_col == 'all5_count' else 0.35)
        ax.text(i, y, f'obs={int(o)}\nP={p:.4f}', ha='center', va='bottom', fontsize=9)

    ax.set_title(title, loc='left', fontweight='bold')
    ax.set_ylabel('Orthogroup count')
    ax.set_xlabel('')
    ax.spines[['top', 'right']].set_visible(False)

axes[-1].set_xlabel('Trait')
axes[-1].set_xticklabels([labels[t] for t in order], rotation=40, ha='right')
axes[0].tick_params(labelbottom=False)
axes[0].legend(loc='upper right', frameon=True)
fig.suptitle('Observed shared orthogroup recurrence exceeds matched-null expectations', fontsize=15, fontweight='bold', y=0.98)
fig.text(0.5, 0.01, 'AM rel = AM fungal relative abundance colonization; AM roots = AM fungal roots colonized', ha='center', fontsize=10)
fig.tight_layout(rect=[0, 0.03, 1, 0.96])

png = OUTDIR / 'SuppFig_null_permutation_orthogroup_recurrence.png'
pdf = OUTDIR / 'SuppFig_null_permutation_orthogroup_recurrence.pdf'
fig.savefig(png, dpi=300, bbox_inches='tight')
fig.savefig(pdf, bbox_inches='tight')
print('WROTE', png)
print('WROTE', pdf)
print('WROTE', OUTCSV)
