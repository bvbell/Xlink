'''Plot the degeneracy of the passing macrocycle set, from scripts/degeneracy.py's CSVs.

Two panels, one story: the variants with the best pass rates are the ones that found a single
conformer 50 times, so "passed" badly overstates how many distinct designs exist.

  top    per variant -- passing structures vs distinct conformers (<1.0 A backbone RMSD).
         The diagonal is "every structure is its own conformer"; everything sits far below it.
  bottom per scaffold -- the funnel from passing structures, to conformers distinct within a
         variant, to conformers still distinct once the five turn seeds are pooled.

Table view is the source CSVs; every plotted value is a column there.
'''
import collections
import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

SURFACE, GRID, AXIS = '#fcfcfb', '#e1e0d9', '#c3c2b7'
INK, INK_2, MUTED = '#0b0b0b', '#52514e', '#898781'
SCAFFOLD = {'short': '#2a78d6', 'medium': '#eb6834', 'long': '#1baf7a'}   # categorical slots 1-3
STAGE = ('#86b6ef', '#2a78d6', '#104281')                                 # ordinal ramp, light->dark
STAGE_LABEL = ('Passing structures', 'Distinct within a variant', 'Distinct across turn seeds')
TARGET_LABEL = {'apj': 'APJ (apelin receptor)', 'uor': 'uOR (mu-opioid receptor)'}


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=MUTED, labelsize=9, length=0)


def scatter_panel(ax, rows, target):
    style(ax)
    ax.grid(axis='both', color=GRID, linewidth=1.0, linestyle='-', zorder=0)
    ax.set_axisbelow(True)
    lim = 53
    ax.plot([0, lim], [0, lim], color=AXIS, linewidth=1.0, zorder=1)
    # The diagonal leaves the top of the panel; label it where it does, not where it would end.
    ax.annotate('every structure\na distinct conformer', xy=(15.4, 15.4), xytext=(18.5, 14.2),
                color=MUTED, fontsize=8, ha='left', va='top',
                arrowprops=dict(arrowstyle='-', color=AXIS, linewidth=1.0, shrinkA=2, shrinkB=2))
    for scaffold, colour in SCAFFOLD.items():
        pts = [r for r in rows if r['scaffold'] == scaffold]
        ax.scatter([int(r['passed']) for r in pts], [int(r['clusters_1.00']) for r in pts],
                   s=90, c=colour, edgecolors=SURFACE, linewidths=2.0, zorder=3,
                   label='%s motif' % scaffold)
    ax.set_xlim(-2, lim)
    ax.set_ylim(-0.7, lim * 0.32)
    ax.set_xlabel('Passing structures (of 50)', color=INK_2, fontsize=9.5)
    ax.set_title(TARGET_LABEL[target], color=INK, fontsize=11, loc='left', pad=10)


def funnel_panel(ax, variants, cross):
    style(ax)
    ax.grid(axis='x', color=GRID, linewidth=1.0, linestyle='-', zorder=0)
    ax.set_axisbelow(True)

    passed, within, across = (collections.Counter() for _ in range(3))
    for r in variants:
        key = (r['target'], r['scaffold'])
        passed[key] += int(r['passed'])
        within[key] += int(r['clusters_1.00'])
    for r in cross:
        across[(r['target'], r['scaffold'])] += int(r['clusters_1.00'])

    keys = [(t, s) for t in ('apj', 'uor') for s in ('short', 'medium', 'long')][::-1]
    for y, key in enumerate(keys):
        values = (passed[key], within[key], across[key])
        ax.plot(values, [y, y, y], color=STAGE[1], linewidth=2.0, solid_capstyle='round', zorder=2)
        for value, colour in zip(values, STAGE):
            ax.scatter([value], [y], s=110, c=colour, edgecolors=SURFACE, linewidths=2.0, zorder=3)
        ax.annotate('%d' % values[0], (values[0], y), xytext=(9, 0), textcoords='offset points',
                    va='center', color=INK_2, fontsize=9)
        ax.annotate('%d' % values[2], (values[2], y), xytext=(-9, 0), textcoords='offset points',
                    va='center', ha='right', color=INK_2, fontsize=9)

    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels(['%s  %s' % (t.upper(), s) for t, s in keys], color=INK_2, fontsize=9.5)
    ax.set_xlim(-28, 440)
    ax.set_ylim(-0.7, len(keys) - 0.3)
    ax.set_xlabel('Structures', color=INK_2, fontsize=9.5)
    ax.set_title('The same collapse per scaffold: 1,200 passing structures are ~60 conformers',
                 color=INK, fontsize=11, loc='left', pad=10)
    handles = [plt.Line2D([], [], marker='o', linestyle='', markersize=9, markerfacecolor=c,
                          markeredgecolor=SURFACE, markeredgewidth=2, label=l)
               for c, l in zip(STAGE, STAGE_LABEL)]
    ax.legend(handles=handles, loc='lower right', frameon=False, fontsize=9, labelcolor=INK_2,
              handletextpad=0.4, borderaxespad=0.8)


def main(out):
    variants, cross = load('degeneracy_variants.csv'), load('degeneracy_cross.csv')
    fig = plt.figure(figsize=(11, 8.2), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.45, wspace=0.16,
                          left=0.115, right=0.97, top=0.83, bottom=0.08)

    axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
    for ax, target in zip(axes, ('apj', 'uor')):
        scatter_panel(ax, [r for r in variants if r['target'] == target], target)
    axes[0].set_ylabel('Distinct conformers\n(<1.0 A backbone RMSD)', color=INK_2, fontsize=9.5)
    axes[0].legend(loc='upper left', frameon=False, fontsize=9, labelcolor=INK_2,
                   handletextpad=0.2, borderaxespad=0.2)
    axes[1].sharey(axes[0])

    funnel_panel(fig.add_subplot(gs[1, :]), variants, cross)

    fig.suptitle('Passing structures are not distinct designs', x=0.115, y=0.955, ha='left',
                 color=INK, fontsize=14.5)
    fig.text(0.115, 0.905, 'Each variant runs 50 GenKIC trajectories on one fixed sequence. '
             'The best-passing variants converge on a single backbone.',
             ha='left', color=INK_2, fontsize=10)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print('wrote %s' % out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'degeneracy.png'))
