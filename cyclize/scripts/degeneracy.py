'''How many genuinely distinct designs are in the passing set?

Every variant runs 50 independent GenKIC trajectories on ONE fixed sequence, so the 1200 passers are
at most 55 distinct molecules, and within a variant they are conformers that may or may not differ.
This reports both axes:

  * sequence degeneracy -- unique one-letter sequences across variants;
  * conformational degeneracy -- greedy clustering of passing structures by backbone RMSD, within a
    variant and then across variants of the same length (same scaffold+route, different turn seed).

Structures are compared after optimal superposition on their own backbone atoms; Rosetta does not
preserve the input frame, so raw coordinates are meaningless (same reason validate.py superimposes).
'''
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cyclize import struct

BB = ('N', 'CA', 'C')
THRESHOLDS = (0.25, 0.5, 1.0)


def backbone(path):
    res = struct.chain_residues(struct.read(path), 'A')
    coords = [r.atoms[a] for r in res for a in BB if a in r.atoms]
    if len(coords) != 3 * len(res):
        raise struct.StructError('%s: incomplete backbone' % path)
    return struct.sequence(res), np.asarray(coords, float)


def pairwise_rmsd(coords):
    '''Symmetric RMSD matrix for a stack of (natoms, 3) arrays, each optimally superimposed.'''
    n = len(coords)
    X = np.stack([c - c.mean(axis=0) for c in coords])
    natoms = X.shape[1]
    out = np.zeros((n, n))
    for i in range(n):
        H = np.einsum('ab,nac->nbc', X[i], X[i + 1:])
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(np.einsum('nab,nbc->nac', Vt.transpose(0, 2, 1), U.transpose(0, 2, 1))))
        S[:, 2] *= d
        sq = (X[i] ** 2).sum() + (X[i + 1:] ** 2).sum(axis=(1, 2)) - 2 * S.sum(axis=1)
        out[i, i + 1:] = out[i + 1:, i] = np.sqrt(np.maximum(sq, 0) / natoms)
    return out


def cluster(d, cutoff):
    '''Greedy leader clustering -> [(centre, [members])], largest cluster first.

    Repeatedly takes the most-connected unassigned point as a centre and absorbs everything within
    `cutoff` of it. Deterministic: ties break on the lowest index.
    '''
    unassigned = set(range(len(d)))
    out = []
    while unassigned:
        idx = sorted(unassigned)
        pop = max(idx, key=lambda i: (sum(1 for j in idx if d[i, j] <= cutoff), -i))
        members = [j for j in idx if d[pop, j] <= cutoff]
        out.append((pop, members))
        unassigned -= set(members)
    return sorted(out, key=lambda cm: -len(cm[1]))


def main(targets, csv_dir=None):
    per_variant, cross = [], []
    for target in targets:
        root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), target)
        rows = [r for r in csv.DictReader(open(os.path.join(root, 'analysis_output/results.csv')))
                if r['pass'] == 'True']
        by_variant = {}
        for r in rows:
            by_variant.setdefault(r['variant'], []).append(r['structure'])

        print('\n=== %s: %d passing structures in %d variants ===' % (target, len(rows), len(by_variant)))
        seqs, loaded, totals = {}, {}, {c: 0 for c in THRESHOLDS}
        print('%-32s %5s %s' % ('variant', 'pass', '  '.join('<%.2fA' % c for c in THRESHOLDS)))
        for variant in sorted(by_variant):
            coords = []
            for name in sorted(by_variant[variant]):
                seq, xyz = backbone(os.path.join(root, 'output', variant, name))
                seqs.setdefault(seq, []).append(variant)
                coords.append(xyz)
            loaded[variant] = (seq, coords)
            d = pairwise_rmsd(coords)
            counts = [cluster(d, c) for c in THRESHOLDS]
            for c, cl in zip(THRESHOLDS, counts):
                totals[c] += len(cl)
            print('%-32s %5d %s   top cluster %d/%d' % (
                variant, len(coords), '  '.join('%5d' % len(s) for s in counts),
                len(counts[-1][0][1]), len(coords)))
            row = {'target': target, 'variant': variant, 'scaffold': variant.split('__')[0],
                   'route': variant.split('__')[1], 'turn': variant.split('__')[2],
                   'residues': len(coords[0]) // 3, 'passed': len(coords),
                   'top_cluster': counts[-1][0]}
            row.update(('clusters_%.2f' % c, len(s)) for c, s in zip(THRESHOLDS, counts))
            per_variant.append(row)

        print('-- unique sequences: %d across %d variants' % (len(seqs), len(by_variant)))
        for seq, vs in sorted(seqs.items(), key=lambda kv: -len(set(kv[1]))):
            if len(set(vs)) > 1:
                print('   SHARED by %d variants: %s  %s' % (len(set(vs)), seq, sorted(set(vs))))
        print('-- within-variant clusters, summed: %s'
              % '  '.join('%.2fA=%d' % (c, totals[c]) for c in THRESHOLDS))

        # Cross-variant: variants of the same length are the same scaffold with a different turn seed.
        groups = {}
        for variant, (seq, coords) in loaded.items():
            groups.setdefault(len(coords[0]), []).append((variant, coords))
        print('-- cross-variant (same backbone length = same scaffold, different turn):')
        for natoms in sorted(groups):
            flat, labels = [], []
            for variant, coords in sorted(groups[natoms]):
                flat.extend(coords)
                labels.extend([variant] * len(coords))
            d = pairwise_rmsd(flat)
            sizes = [len(m) for _, m in cluster(d, 1.0)]
            print('   %3d residues  %2d variants  %4d structures -> %4d clusters <1.0A (largest %d)'
                  % (natoms // 3, len(groups[natoms]), len(flat), len(sizes), sizes[0]))
            cross.append({'target': target, 'residues': natoms // 3,
                          'scaffold': sorted(groups[natoms])[0][0].split('__')[0],
                          'variants': len(groups[natoms]), 'passed': len(flat),
                          'clusters_1.00': len(sizes), 'largest': sizes[0]})

    if csv_dir:
        for name, rows in (('degeneracy_variants.csv', per_variant), ('degeneracy_cross.csv', cross)):
            path = os.path.join(csv_dir, name)
            with open(path, 'w') as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
            print('wrote %s' % path)


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    out = [a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--csv-dir=')]
    main(args or ['apj', 'uor'], out[0] if out else None)
