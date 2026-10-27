#!/usr/bin/env python
"""Embed the all-vs-all TM matrix in 2D with UMAP and t-SNE.

Both run on the precomputed distance D = 1 - TM, so no feature vectors are
invented: the pairwise structural distance is the only input.

Run with the protein-shapes env, which is the only one here with umap:
  /home/groups/possu/miniconda/envs/protein-shapes/bin/python embed_diversity.py
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE

# dataviz skill palette (validated: ordinal blue --ordinal, categorical
# first-3 slots --pairs all, both on the light surface #fcfcfb).
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
# 5 ordinal steps, blue 250->650; light end clears 2:1 on the surface.
LINKER_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
CLUSTER_HUES = ["#2a78d6", "#eb6834", "#1baf7a"]  # slots 1-3
OTHER = "#c3c2b7"


def embed(D, seed, perplexity, n_neighbors, min_dist):
    """Return (tsne_xy, umap_xy); umap_xy is None if umap is unavailable."""
    tsne = TSNE(
        n_components=2, metric="precomputed", init="random",
        perplexity=perplexity, random_state=seed,
    ).fit_transform(D)

    try:
        import umap
    except ImportError:
        print("umap not importable - t-SNE only")
        return tsne, None
    # densmap=False, euclidean output; precomputed input metric.
    um = umap.UMAP(
        n_components=2, metric="precomputed", n_neighbors=n_neighbors,
        min_dist=min_dist, random_state=seed,
    ).fit_transform(D)
    return tsne, um


def _panel(ax, xy, colors, title):
    ax.set_facecolor(SURFACE)
    # 8px+ markers with a surface ring, so overlapping points stay countable.
    ax.scatter(
        xy[:, 0], xy[:, 1], s=46, c=colors,
        edgecolors=SURFACE, linewidths=0.8, zorder=3,
    )
    ax.set_title(title, fontsize=10, color=INK, loc="left", pad=7)
    # The embedding axes carry no units; ticks would invite false reading.
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color("#e1e0d9")


def plot(names, linker, cluster, tsne, um, out_png):
    uniq_link = sorted(set(linker.tolist()))
    link_color = {v: LINKER_RAMP[i % len(LINKER_RAMP)]
                  for i, v in enumerate(uniq_link)}

    # 42 clusters is far past the 3-slot cap for an all-pairs form, so the
    # tail folds into one gray "Other" rather than cycling hues.
    sizes = {c: int((cluster == c).sum()) for c in set(cluster.tolist())}
    top = sorted(sizes, key=lambda c: -sizes[c])[:3]
    clus_color = {c: CLUSTER_HUES[i] for i, c in enumerate(top)}

    methods = [("t-SNE", tsne)] + ([("UMAP", um)] if um is not None else [])
    fig, axes = plt.subplots(
        len(methods), 2, figsize=(10.5, 4.9 * len(methods)), squeeze=False,
    )
    fig.patch.set_facecolor(SURFACE)

    for row, (label, xy) in enumerate(methods):
        _panel(
            axes[row][0], xy, [link_color[v] for v in linker],
            f"{label} - colored by linker length",
        )
        _panel(
            axes[row][1], xy,
            [clus_color.get(c, OTHER) for c in cluster],
            f"{label} - colored by TM cluster (top 3 of {len(sizes)})",
        )
        # Direct labels on the big blobs: the aqua slot sits below 3:1 on this
        # surface, so identity cannot rest on the swatch alone.
        for c in top:
            m = cluster == c
            axes[row][1].text(
                xy[m, 0].mean(), xy[m, 1].mean(), f"n={sizes[c]}",
                fontsize=8, color=INK, ha="center", va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec="none",
                          alpha=0.85),
            )

    handles = [
        plt.Line2D([], [], marker="o", ls="", markersize=7,
                   markerfacecolor=link_color[v], markeredgecolor=SURFACE,
                   label=f"linker {v}")
        for v in uniq_link
    ]
    axes[0][0].legend(
        handles=handles, loc="upper left", bbox_to_anchor=(0, -0.02), ncol=5,
        frameon=False, fontsize=8, labelcolor=INK_2, handletextpad=0.3,
        columnspacing=1.0,
    )
    handles = [
        plt.Line2D([], [], marker="o", ls="", markersize=7,
                   markerfacecolor=clus_color[c], markeredgecolor=SURFACE,
                   label=f"cluster {c} (n={sizes[c]})")
        for c in top
    ] + [plt.Line2D([], [], marker="o", ls="", markersize=7,
                    markerfacecolor=OTHER, markeredgecolor=SURFACE,
                    label=f"other ({len(sizes) - 3} clusters)")]
    axes[0][1].legend(
        handles=handles, loc="upper left", bbox_to_anchor=(0, -0.02), ncol=2,
        frameon=False, fontsize=8, labelcolor=INK_2, handletextpad=0.3,
        columnspacing=1.0,
    )

    fig.suptitle(
        "Design diversity: 2D embedding of the all-vs-all TM matrix",
        fontsize=12.5, color=INK, x=0.02, ha="left", y=0.985,
    )
    fig.text(
        0.02, 0.955,
        f"{len(names)} designs, distance = 1 - TM-score. "
        "Axes are arbitrary; only neighborhoods carry meaning.",
        fontsize=9, color=MUTED, ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(out_png, dpi=200, facecolor=SURFACE)
    print("wrote", out_png)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--npz", default="runs_prod/diversity_pairs.npz")
    p.add_argument("--out_dir", default="runs_prod")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--perplexity", type=float, default=30.0)
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--min_dist", type=float, default=0.1)
    a = p.parse_args()

    d = np.load(a.npz, allow_pickle=True)
    tm = np.asarray(d["tm"], dtype=float)
    tm = (tm + tm.T) / 2.0
    D = 1.0 - tm
    np.fill_diagonal(D, 0.0)
    D[D < 0] = 0.0

    names, linker, cluster = d["names"], d["linker"], d["cluster"]
    print(f"{len(names)} designs, TM {tm.min():.3f}-{tm.max():.3f}")

    tsne, um = embed(D, a.seed, a.perplexity, a.n_neighbors, a.min_dist)

    csv = os.path.join(a.out_dir, "diversity_embedding.csv")
    with open(csv, "w") as fh:
        fh.write("design,linker,cluster,tsne1,tsne2,umap1,umap2\n")
        for i, n in enumerate(names):
            u = (f"{um[i,0]:.4f},{um[i,1]:.4f}" if um is not None else ",")
            fh.write(f"{n},{linker[i]},{cluster[i]},"
                     f"{tsne[i,0]:.4f},{tsne[i,1]:.4f},{u}\n")
    print("wrote", csv)

    plot(names, linker, cluster, tsne, um,
         os.path.join(a.out_dir, "diversity_embedding.png"))


if __name__ == "__main__":
    main()
