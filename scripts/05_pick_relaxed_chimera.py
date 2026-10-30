'''Pick the best constrained-FastRelax decoy of the CLDN6/9v32-Fab chimera.

The raw superposition chimera has 56 heavy-atom pairs under 2.0 A against the Fab, several of them
backbone-backbone. relax.sbatch produced 20 decoys at two coordinate-constraint strengths (sd 0.5 and
sd 1.0). This scores each decoy on what actually matters for a ddG run:

  - did the clashes clear?                  (heavy-atom pairs A vs H+L under 2.0/2.5/3.0 A)
  - did the structure stay put?             (CA RMSD to the input chimera, all and interface-only)
  - are the chain-H gaps still gaps?        (they must not have been zipped shut by the jump)
  - total talaris2014 score

Selection: among decoys with zero sub-2.0 A pairs, take the lowest total score. If none are clean,
report that and take the fewest clashes, breaking ties on score.

Writes the winner to cldn6_9v32fab_relaxed_complex.pdb with hydrogens stripped and Rosetta's trailing
score block dropped, so it reads like the other targets' inputs.
'''
import glob, os, sys

import numpy as np

ROOT = '/scratch/users/bvbell/xlink'
RELAX = os.path.join(ROOT, 'flexddg/cldn6_9v32fab/relax/decoys')
START = os.path.join(ROOT, 'cldn6_9v32fab_complex.pdb')
OUT = os.path.join(ROOT, 'cldn6_9v32fab_relaxed_complex.pdb')


def read_atoms(path, heavy_only=True):
    out = []
    for l in open(path):
        if not l.startswith(('ATOM', 'HETATM')):
            continue
        if l[16] not in (' ', 'A'):
            continue
        if heavy_only and (l[76:78].strip() == 'H' or l[12:16].strip().startswith('H')):
            continue
        out.append(l.rstrip('\n'))
    return out


def xyz(l):
    return (float(l[30:38]), float(l[38:46]), float(l[46:54]))


def key(l):
    return (l[21], int(l[22:26]), l[26], l[12:16])


def clashes(atoms, cutoffs=(2.0, 2.5, 3.0)):
    a = [l for l in atoms if l[21] == 'A']
    f = [l for l in atoms if l[21] in 'HL']
    XA, XF = np.array([xyz(l) for l in a]), np.array([xyz(l) for l in f])
    d = np.sqrt(((XA[:, None, :] - XF[None, :, :]) ** 2).sum(-1))
    worst = []
    if d.size:
        i, j = np.unravel_index(np.argmin(d), d.shape)
        worst = [a[i][17:26].strip(), f[j][17:26].strip(), float(d[i, j])]
    return [int((d < c).sum()) for c in cutoffs], worst


def ca_map(atoms):
    return {key(l)[:3]: np.array(xyz(l)) for l in atoms if l[12:16] == ' CA '}


def rmsd(m1, m2, subset=None):
    ks = [k for k in m1 if k in m2 and (subset is None or k in subset)]
    if not ks:
        return float('nan'), 0
    d = np.array([m1[k] - m2[k] for k in ks])
    return float(np.sqrt((d ** 2).sum(1).mean())), len(ks)


def breaks(atoms):
    C, N, order = {}, {}, []
    for l in atoms:
        k = key(l)[:3]
        if not order or order[-1] != k:
            order.append(k)
        if l[12:16] == ' C  ':
            C[k] = np.array(xyz(l))
        elif l[12:16] == ' N  ':
            N[k] = np.array(xyz(l))
    out = []
    for a, b in zip(order, order[1:]):
        if a[0] == b[0] and a in C and b in N:
            d = float(np.linalg.norm(C[a] - N[b]))
            if d > 2.0:
                out.append((a[0], a[1], b[1], d))
    return out


def total_score(pdb):
    '''Rosetta appends a score table after END; the pose line starts with "pose".'''
    for l in open(pdb):
        if l.startswith('pose '):
            return float(l.split()[-1])
    return float('nan')


def interface_keys(atoms, cutoff=5.0):
    a = [l for l in atoms if l[21] == 'A']
    f = [l for l in atoms if l[21] in 'HL']
    XA, XF = np.array([xyz(l) for l in a]), np.array([xyz(l) for l in f])
    d = np.sqrt(((XA[:, None, :] - XF[None, :, :]) ** 2).sum(-1))
    ks = set()
    for i, j in zip(*np.where(d < cutoff)):
        ks.add(key(a[i])[:3])
        ks.add(key(f[j])[:3])
    return ks


def main():
    start = read_atoms(START)
    s_ca = ca_map(start)
    iface = interface_keys(start)
    s_clash, s_worst = clashes(start)
    print('input chimera: clashes <2.0/<2.5/<3.0 = %s, worst %s-%s %.2f A, breaks %s\n'
          % (s_clash, s_worst[0], s_worst[1], s_worst[2], breaks(start)))

    rows = []
    for pdb in sorted(glob.glob(os.path.join(RELAX, '*', '*', '*.pdb'))):
        var = pdb.split(os.sep)[-3]
        atoms = read_atoms(pdb)
        chains = {l[21] for l in atoms}
        if chains != {'A', 'H', 'L'}:
            print('SKIP %s: chains %s' % (pdb, sorted(chains)))
            continue
        c, worst = clashes(atoms)
        m = ca_map(atoms)
        all_rms, n_all = rmsd(s_ca, m)
        if_rms, n_if = rmsd(s_ca, m, iface)
        rows.append(dict(pdb=pdb, var=var, c20=c[0], c25=c[1], c30=c[2], worst=worst,
                         score=total_score(pdb), rms=all_rms, n=n_all, if_rms=if_rms, n_if=n_if,
                         brk=breaks(atoms)))
    if not rows:
        sys.exit('no decoys found under %s' % RELAX)

    print('%-6s %-28s %6s %6s %6s %10s %7s %8s  %s'
          % ('var', 'decoy', '<2.0', '<2.5', '<3.0', 'score', 'CA rms', 'if rms', 'H gaps'))
    for r in sorted(rows, key=lambda r: (r['c20'], r['score'])):
        print('%-6s %-28s %6d %6d %6d %10.1f %7.2f %8.2f  %s'
              % (r['var'], os.path.basename(r['pdb'])[:28], r['c20'], r['c25'], r['c30'],
                 r['score'], r['rms'], r['if_rms'],
                 ','.join('%s%d-%d:%.1f' % b for b in r['brk']) or 'NONE(!)'))

    clean = [r for r in rows if r['c20'] == 0]
    if clean:
        best = min(clean, key=lambda r: r['score'])
        print('\n%d/%d decoys cleared every sub-2.0 A pair.' % (len(clean), len(rows)))
    else:
        best = min(rows, key=lambda r: (r['c20'], r['score']))
        print('\nWARNING: no decoy cleared every sub-2.0 A pair; best still has %d.' % best['c20'])
    print('chosen: %s (%s)  score %.1f  clashes %d/%d/%d  CA rms %.2f A  interface CA rms %.2f A'
          % (best['pdb'], best['var'], best['score'], best['c20'], best['c25'], best['c30'],
             best['rms'], best['if_rms']))

    with open(OUT, 'w') as f:
        last = None
        for l in read_atoms(best['pdb']):
            if last and l[21] != last:
                f.write('TER\n')
            f.write(l + '\n')
            last = l[21]
        f.write('TER\nEND\n')
    print('wrote %s (%d heavy atoms, hydrogens stripped)' % (OUT, len(read_atoms(best['pdb']))))


if __name__ == '__main__':
    main()
