'''Did the rm_lower phantom actually move geometry, or only inflate the reported score?

Until 2026-10-29 build.xml stripped LOWER_TERMINUS_VARIANT from residue 1 on the disulfide routes,
where nothing is ever bonded to that N. Rosetta scored its rama at phi = 0 -- a forbidden region --
for about +10 REU. The claim that this biased SAMPLING rather than just reporting rests on an
inference: phi(1) is fixed garbage that minimisation cannot touch, but psi(1) is real and free, and
residue 1 sits in the GenKIC loop with backbone enabled in the FastRelax MoveMap. So relax should
have spent effort optimising psi(1) against a bogus phi. Residue 1 is the prepended Cys, and its
backbone orientation places SG -- the disulfide the route exists to form.

This compares the pre-fix runs (kept in <target>/output_prebug/) against the reruns, on the three
things the inference predicts should move: psi(1), the S-S distance, and the closure/pass yield.

Run after the reruns finish. Reports only -- writes nothing.
'''
import argparse
import collections
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from cyclize import struct      # noqa: E402


def psi(residues, num):
    '''psi of residue `num`: N-CA-C-N(next). Returns None at the chain end.'''
    by = {r.num: r for r in residues}
    a, b = by.get(num), by.get(num + 1)
    if not a or not b or not all(x in a.atoms for x in ('N', 'CA', 'C')) or 'N' not in b.atoms:
        return None
    return struct.dihedral(a.atoms['N'], a.atoms['CA'], a.atoms['C'], b.atoms['N'])


def ss_distance(residues):
    '''Closest SG-SG pair, which for these designs is the disulfide.'''
    sgs = [r.atoms['SG'] for r in residues if 'SG' in r.atoms]
    if len(sgs) < 2:
        return None
    return min(struct.dist(a, b) for i, a in enumerate(sgs) for b in sgs[i + 1:])


def survey(d):
    out = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(d, '*.pdb'))):
        try:
            res = struct.chain_residues(struct.read(f), 'A')
        except struct.StructError:
            continue
        p, s = psi(res, res[0].num), ss_distance(res)
        if p is not None:
            out['psi1'].append(p)
        if s is not None:
            out['ss'].append(s)
        out['n'].append(1)
    return out


def stats(xs):
    if not xs:
        return '      n/a       '
    m = sum(xs) / len(xs)
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / len(xs)) if len(xs) > 1 else 0.0
    return '%7.2f +/- %-5.2f' % (m, sd)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('targets', nargs='*', default=['apj', 'uor'])
    args = ap.parse_args(argv)

    for target in args.targets:
        old_root = os.path.join(ROOT, target, 'output_prebug')
        if not os.path.isdir(old_root):
            print('%s: no output_prebug/, nothing to compare' % target)
            continue
        print('\n=== %s ===' % target)
        print('%-30s %6s %6s  %-19s %-19s  %-15s %-15s'
              % ('variant', 'n_pre', 'n_post', 'psi(1) pre', 'psi(1) post',
                 'S-S pre', 'S-S post'))
        agg = {k: ([], []) for k in ('psi1', 'ss')}
        for vid in sorted(os.listdir(old_root)):
            new = os.path.join(ROOT, target, 'output', vid)
            if not os.path.isdir(new):
                continue
            pre, post = survey(os.path.join(old_root, vid)), survey(new)
            for k in agg:
                agg[k][0].extend(pre[k])
                agg[k][1].extend(post[k])
            print('%-30s %6d %6d  %-19s %-19s  %-15s %-15s'
                  % (vid, len(pre['n']), len(post['n']),
                     stats(pre['psi1']), stats(post['psi1']),
                     stats(pre['ss']), stats(post['ss'])))
        print('-' * 118)
        for k, label in (('psi1', 'psi(residue 1)'), ('ss', 'S-S distance')):
            pre, post = agg[k]
            if pre and post:
                print('%-30s pre %s   post %s   shift %+.2f'
                      % (label, stats(pre), stats(post),
                         sum(post) / len(post) - sum(pre) / len(pre)))


if __name__ == '__main__':
    main()
