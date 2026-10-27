'''Figures for an alanine scan, and a side-by-side comparison of targets. Transparent PNG (300 dpi) + PDF.

Chains are colored by categorical slot in the order they appear in target.json `groups`, so a chain
keeps its color whichever mutations are plotted. Only mutations on the `scan.side` groups are shown.
'''
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from . import config  # noqa: E402

TEXT, TEXT2, GRID = '#0b0b0b', '#52514e', '#d6d5d0'
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']

plt.rcParams.update({
    'font.size': 9, 'text.color': TEXT, 'axes.labelcolor': TEXT2,
    'xtick.color': TEXT2, 'ytick.color': TEXT, 'axes.edgecolor': GRID,
    'axes.linewidth': 1, 'xtick.major.size': 0, 'ytick.major.size': 0,
    'figure.facecolor': 'none', 'axes.facecolor': 'none', 'savefig.transparent': True,
})


def chain_styles(cfg):
    '''{chain: (color, legend label)} for every chain, slots in config order.'''
    order = config.chain_order(cfg)
    if len(order) > len(SLOTS):
        raise ValueError('more than %d chains; fold some into one group label' % len(SLOTS))
    return {c: (SLOTS[k], '%s (%s)' % (cfg['chain_labels'].get(c, 'Chain'), c)) for k, c in enumerate(order)}


def side_chains(cfg):
    return [c for g in cfg['scan']['side'] for c in cfg['groups'][g]]


def label(case):
    return case.replace('_', ' ').replace('-', ' + ')


def clean(ax, grid_axis='x'):
    for s in ['top', 'right', 'left']:
        ax.spines[s].set_visible(False)
    ax.grid(axis=grid_axis, color=GRID, lw=1, ls='-')
    ax.set_axisbelow(True)


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(out_dir, name + '.' + ext), dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)


def legend(ax, styles, cases, loc):
    present = {c[0] for c in cases}
    ax.legend(handles=[Patch(color=v[0], label=v[1]) for c, v in styles.items() if c in present],
              loc=loc, frameon=False)


def _bars(ax, s, styles, hotspot, xlim=None):
    y = np.arange(len(s))
    ax.barh(y, s.ddG_GAM, height=0.72, color=[styles[c[0]][0] for c in s.case])
    ax.errorbar(s.ddG_GAM, y, xerr=[s.ddG_GAM - s.GAM_lo, s.GAM_hi - s.ddG_GAM],
                fmt='none', ecolor=TEXT2, elinewidth=1, capsize=0)
    ax.axvline(0, color=TEXT2, lw=1)
    ax.axvline(hotspot, color=TEXT2, lw=1, ls=(0, (1, 2)))
    ax.text(hotspot + 0.05, len(s) - 0.4, 'hotspot ≥ %g kcal/mol' % hotspot, color=TEXT2, fontsize=8, va='bottom')
    ax.set_yticks(y)
    ax.set_yticklabels([label(c) for c in s.case])
    for yi, (v, hi) in enumerate(zip(s.ddG_GAM, s.GAM_hi)):
        if v >= hotspot:
            ax.text(hi + 0.08, yi, '%.1f' % v, va='center', fontsize=8, color=TEXT)
    ax.set_ylim(-0.7, len(s) - 0.1)
    if xlim:
        ax.set_xlim(*xlim)
    clean(ax)


def ranked_bars(cfg, summ, styles, out_dir):
    s = summ.sort_values('ddG_GAM')
    hotspot = cfg['analysis']['hotspot']
    fig, ax = plt.subplots(figsize=(5.2, max(2.4, 0.22 * len(s) + 0.9)))
    _bars(ax, s, styles, hotspot)
    ax.set_xlabel('ΔΔG binding, flex ddG GAM (kcal/mol)\nerror bars: 95%% bootstrap over %d trajectories'
                  % summ.n.max())
    ax.set_title('%s alanine scan: binding hotspots' % cfg['title'], loc='left', fontsize=11, color=TEXT)
    legend(ax, styles, s.case, 'lower right')
    save(fig, out_dir, 'alascan_ranked')


def trajectory_spread(cfg, summ, structs, styles, out_dir):
    order = summ.sort_values('ddG_GAM', ascending=False).case.tolist()
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(max(3.5, 0.34 * len(order)), 3.6))
    for i, c in enumerate(order):
        v = structs.loc[structs.case_name == c, 'total_score'].values
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(v)), v, s=9, color=styles[c[0]][0], alpha=0.55, linewidths=0)
        ax.plot([i - 0.3, i + 0.3], [v.mean()] * 2, color=TEXT, lw=2, solid_capstyle='round')
    ax.axhline(0, color=TEXT2, lw=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label(c) for c in order], rotation=60, ha='right', fontsize=8)
    ax.tick_params(axis='x', colors=TEXT)
    ax.set_xlim(-0.7, len(order) - 0.3)
    ax.set_ylabel('ΔΔG per trajectory (talaris2014, REU)')
    ax.set_title('Per-trajectory ΔΔG (dots) and mean (bar), ordered by GAM ΔΔG', loc='left', fontsize=11, color=TEXT)
    legend(ax, styles, order, 'upper right')
    clean(ax, 'y')
    ax.spines['left'].set_visible(False)
    save(fig, out_dir, 'alascan_trajectories')


def term_heatmap(cfg, summ, structs, out_dir):
    terms = ['fa_atr', 'fa_sol', 'fa_rep', 'fa_elec', 'hbond_sc', 'hbond_bb_sc']
    names = ['Attraction\n(fa_atr)', 'Solvation\n(fa_sol)', 'Repulsion\n(fa_rep)',
             'Electrostatics\n(fa_elec)', 'Side-chain\nH-bond', 'Backbone–SC\nH-bond']
    hot = summ[summ.ddG_GAM >= cfg['analysis']['hotspot']].sort_values('ddG_GAM', ascending=False).case.tolist()
    if not hot:
        print('no hotspots, skipping term heatmap')
        return
    m = structs[structs.case_name.isin(hot)].groupby('case_name')[terms].mean().loc[hot]
    lim = np.abs(m.values).max()
    cmap = LinearSegmentedColormap.from_list('div', ['#1c5cab', '#86b6ef', '#f0efec', '#f0a0a0', '#b8292a'])
    fig, ax = plt.subplots(figsize=(6.2, max(1.8, 0.37 * len(hot) + 0.9)))
    norm = TwoSlopeNorm(0, -lim, lim)
    # Separate rectangles rather than imshow, so the gaps between cells stay transparent.
    for (i, j), v in np.ndenumerate(m.values):
        ax.add_patch(Rectangle((j - 0.47, i - 0.46), 0.94, 0.92, color=cmap(norm(v)), lw=0))
        ax.text(j, i, '%.1f' % v, ha='center', va='center', fontsize=8,
                color='#ffffff' if abs(v) > 0.6 * lim else TEXT)
    ax.set_xlim(-0.5, len(terms) - 0.5)
    ax.set_ylim(len(hot) - 0.5, -0.5)
    ax.set_xticks(range(len(terms)))
    ax.set_xticklabels(names, fontsize=8)
    ax.set_yticks(range(len(hot)))
    ax.set_yticklabels([label(c) for c in hot])
    ax.tick_params(axis='x', colors=TEXT, top=True, labeltop=True, bottom=False, labelbottom=False)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.04, pad=0.03)
    cb.set_label('ΔΔG contribution (REU)\nred: alanine loses this', color=TEXT2)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=TEXT2, size=0)
    ax.set_title('What each hotspot loses: per-term ΔΔG (talaris, before GAM)', loc='left', fontsize=11,
                 color=TEXT, pad=34)
    save(fig, out_dir, 'alascan_hotspot_terms')


def load_side(cfg, analysis_dir=None):
    d = analysis_dir or os.path.join(cfg['dir'], 'analysis_output')
    summ = pd.read_csv(os.path.join(d, 'alascan_summary.csv'))
    structs = pd.read_csv(os.path.join(d, 'struct_scores.csv'))
    side = side_chains(cfg)
    return summ[summ.case.str[0].isin(side)], structs[structs.case_name.str[0].isin(side)]


def plot_target(cfg, analysis_dir=None, out_dir=None):
    out_dir = out_dir or os.path.join(cfg['dir'], 'plots')
    summ, structs = load_side(cfg, analysis_dir)
    styles = chain_styles(cfg)
    ranked_bars(cfg, summ, styles, out_dir)
    trajectory_spread(cfg, summ, structs, styles, out_dir)
    term_heatmap(cfg, summ, structs, out_dir)
    print('wrote', out_dir)


def compare(cfgs, out_dir):
    '''One ranked-bar panel per target on a shared x scale, plus hotspots.csv.'''
    data = [(cfg, load_side(cfg)[0].sort_values('ddG_GAM')) for cfg in cfgs]
    lo = min(0, min(s.GAM_lo.min() for _, s in data)) - 0.2
    hi = max(s.GAM_hi.max() for _, s in data) + 0.6
    n = max(len(s) for _, s in data)
    fig, axes = plt.subplots(1, len(data), figsize=(4.4 * len(data), max(2.4, 0.22 * n + 0.9)), squeeze=False)
    for ax, (cfg, s) in zip(axes[0], data):
        styles = chain_styles(cfg)
        _bars(ax, s, styles, cfg['analysis']['hotspot'], (lo, hi))
        ax.set_ylim(-0.7, n - 0.1)
        ax.set_title(cfg['title'], loc='left', fontsize=11, color=TEXT)
        ax.set_xlabel('ΔΔG binding, GAM (kcal/mol)')
        legend(ax, styles, s.case, 'lower right')
    fig.suptitle('Alanine scan hotspots by target (95% bootstrap CI)', x=0.01, ha='left', fontsize=12, color=TEXT)
    save(fig, out_dir, 'compare_ranked')
    rows = [dict(target=cfg['title'], case=r.case, ddG_GAM=r.ddG_GAM, GAM_lo=r.GAM_lo, GAM_hi=r.GAM_hi)
            for cfg, s in data for r in s.itertuples() if r.ddG_GAM >= cfg['analysis']['hotspot']]
    hot = pd.DataFrame(rows).sort_values(['target', 'ddG_GAM'], ascending=[True, False])
    hot.to_csv(os.path.join(out_dir, 'hotspots.csv'), index=False)
    print(hot.round(2).to_string(index=False))
    print('wrote', out_dir)
