#!/usr/bin/env python3
'''Plots for the Fab-side flex ddG alanine scan. All figures are saved with transparent backgrounds.

Reads analysis_output/fab_alascan_summary.csv and fab_complete-struct_scores_results.csv,
writes PNG (300 dpi) + PDF to plots/.
'''

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
from matplotlib.cm import ScalarMappable

OUT = 'plots'
TEXT, TEXT2, GRID = '#0b0b0b', '#52514e', '#d6d5d0'
CHAIN = {'A': ('#2a78d6', 'Heavy chain (A)'), 'B': ('#eb6834', 'Light chain (B)')}
HOTSPOT = 1.0  # kcal/mol

plt.rcParams.update({
    'font.size': 9, 'text.color': TEXT, 'axes.labelcolor': TEXT2,
    'xtick.color': TEXT2, 'ytick.color': TEXT, 'axes.edgecolor': GRID,
    'axes.linewidth': 1, 'xtick.major.size': 0, 'ytick.major.size': 0,
    'figure.facecolor': 'none', 'axes.facecolor': 'none', 'savefig.transparent': True,
})

def label(case):
    chain, mut = case.split('_')
    return '%s %s' % (chain, mut)

def clean(ax, grid_axis='x'):
    for s in ['top', 'right', 'left']:
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, lw=1, ls='-')
    ax.set_axisbelow(True)

def save(fig, name):
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT, name + '.' + ext), dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)

def ranked_bars(summ):
    s = summ.sort_values('ddG_GAM')
    y = np.arange(len(s))
    colors = [CHAIN[c[0]][0] for c in s.case]
    fig, ax = plt.subplots(figsize=(5.2, 6.4))
    ax.barh(y, s.ddG_GAM, height=0.72, color=colors)
    ax.errorbar(s.ddG_GAM, y, xerr=[s.ddG_GAM - s.GAM_lo, s.GAM_hi - s.ddG_GAM],
                fmt='none', ecolor=TEXT2, elinewidth=1, capsize=0)
    ax.axvline(0, color=TEXT2, lw=1)
    ax.axvline(HOTSPOT, color=TEXT2, lw=1, ls=(0, (1, 2)))
    ax.text(HOTSPOT + 0.05, len(s) - 0.4, 'hotspot ≥ 1 kcal/mol', color=TEXT2, fontsize=8, va='bottom')
    ax.set_yticks(y)
    ax.set_yticklabels([label(c) for c in s.case])
    for yi, (v, hi) in enumerate(zip(s.ddG_GAM, s.GAM_hi)):
        if v >= HOTSPOT:
            ax.text(hi + 0.08, yi, '%.1f' % v, va='center', fontsize=8, color=TEXT)
    ax.set_ylim(-0.7, len(s) - 0.1)
    ax.set_xlabel('ΔΔG binding, flex ddG GAM (kcal/mol)\nerror bars: 95% bootstrap over 35 trajectories')
    ax.set_title('Fab alanine scan: binding hotspots', loc='left', fontsize=11, color=TEXT)
    ax.legend(handles=[Patch(color=v[0], label=v[1]) for v in CHAIN.values()],
              loc='lower right', frameon=False)
    clean(ax)
    save(fig, 'fab_alascan_ranked')

def trajectory_spread(summ, structs):
    order = summ.sort_values('ddG_GAM', ascending=False).case.tolist()
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    for i, c in enumerate(order):
        v = structs.loc[structs.case_name == c, 'total_score'].values
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=9, color=CHAIN[c[0]][0],
                   alpha=0.55, linewidths=0)
        ax.plot([i - 0.3, i + 0.3], [v.mean()] * 2, color=TEXT, lw=2, solid_capstyle='round')
    ax.axhline(0, color=TEXT2, lw=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label(c) for c in order], rotation=60, ha='right', fontsize=8)
    ax.tick_params(axis='x', colors=TEXT)
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel('ΔΔG per trajectory (talaris2014, REU)')
    ax.set_title('Per-trajectory ΔΔG (dots) and mean (bar), ordered by GAM ΔΔG',
                 loc='left', fontsize=11, color=TEXT)
    ax.legend(handles=[Patch(color=v[0], label=v[1]) for v in CHAIN.values()],
              loc='upper right', frameon=False)
    clean(ax, 'y')
    ax.spines['left'].set_visible(False)
    save(fig, 'fab_alascan_trajectories')

def term_heatmap(summ, structs):
    terms = ['fa_atr', 'fa_sol', 'fa_rep', 'fa_elec', 'hbond_sc', 'hbond_bb_sc']
    names = ['Attraction\n(fa_atr)', 'Solvation\n(fa_sol)', 'Repulsion\n(fa_rep)',
             'Electrostatics\n(fa_elec)', 'Side-chain\nH-bond', 'Backbone–SC\nH-bond']
    hot = summ[summ.ddG_GAM >= HOTSPOT].sort_values('ddG_GAM', ascending=False).case.tolist()
    m = structs[structs.case_name.isin(hot)].groupby('case_name')[terms].mean().loc[hot]
    lim = np.abs(m.values).max()
    cmap = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#f0efec', '#f0a0a0', '#b8292a'])
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    norm = TwoSlopeNorm(0, -lim, lim)
    # Separate rectangles rather than imshow, so the gaps between cells stay transparent.
    for (i, j), v in np.ndenumerate(m.values):
        ax.add_patch(Rectangle((j - 0.47, i - 0.46), 0.94, 0.92, color=cmap(norm(v)), lw=0))
    ax.set_xlim(-0.5, len(terms) - 0.5)
    ax.set_ylim(len(hot) - 0.5, -0.5)
    im = ScalarMappable(norm=norm, cmap=cmap)
    for (i, j), v in np.ndenumerate(m.values):
        ax.text(j, i, '%.1f' % v, ha='center', va='center', fontsize=8,
                color='#ffffff' if abs(v) > 0.6 * lim else TEXT)
    ax.set_xticks(range(len(terms)))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_yticks(range(len(hot)))
    ax.set_yticklabels([label(c) for c in hot])
    ax.tick_params(axis='x', colors=TEXT, top=True, labeltop=True, bottom=False, labelbottom=False)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label('ΔΔG contribution (REU)\nred: alanine loses this', color=TEXT2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TEXT2, size=0)
    ax.set_title('What each hotspot loses: per-term ΔΔG (talaris, before GAM)',
                 loc='left', fontsize=11, color=TEXT, pad=34)
    save(fig, 'fab_alascan_hotspot_terms')

if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    summ = pd.read_csv('analysis_output/fab_alascan_summary.csv')
    structs = pd.read_csv('analysis_output/fab_complete-struct_scores_results.csv')
    ranked_bars(summ)
    trajectory_spread(summ, structs)
    term_heatmap(summ, structs)
    print('wrote', sorted(os.listdir(OUT)))
