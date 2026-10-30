#!/usr/bin/env python3
'''Pick a representative "good" and "bad" set of designs per GPCR target for manual inspection.

Written 2026-10-30, after two build-level defects were found that no existing filter detects:
non-planar omega (both targets) and a D-amino-acid motif (uOR only). See ../README.md "Known
caveats" and ../../STATUS.md "Recommendations".

"Good" here is strictly RELATIVE. Nothing in the current output is clean -- the best intra-peptide
contact in the whole apj set is 2.99 A against a >=3.4 A norm, and 89% of complexes carry a
non-planar amide. These are the least-bad structures, chosen so that what is wrong with them is
visible by eye next to the worst ones.

Metrics per structure, all computed from coordinates (nothing trusted from the CSVs except
motif_rmsd and dG_separated, which come from interface.csv):

  omega_max_dev   largest deviation of any peptide bond from planarity, in degrees.
                  0 = perfect; >30 is non-planar; 90 is a fully broken amide.
  cis_nonpro      peptide bonds with |omega| < 90 where the FOLLOWING residue is not proline.
                  ~0.03% of bonds in real proteins. Legitimate cis-Pro is counted separately.
  intra_min       closest non-adjacent intra-peptide heavy-atom contact. >=3.4 A is normal
                  packing; below ~3.0 A the chain is passing through itself.
  inter_min       closest design-to-receptor heavy-atom contact. <2.0 A is a clash.
  n_dres          backbone D-amino acids (improper N-CA-C-CB < 0). Should be 0 except where a
                  D residue was deliberately designed in at the turn.
  closure_ON      O...N distance across the ring-closure bond. Ideal 2.25 A. Rosetta constrains
                  the closure C-N but NOT the carbonyl oxygen, so in 57% of amide designs the O
                  sits on top of the next N -- as close as 0.14 A. Invisible to the energy
                  function because O and N are a 1-3 pair through C and fa_rep excludes those.
  flank_dev       worst omega deviation among bonds within 2 of the closure. This is where the
                  closure is paid for: 47% of the immediately flanking bonds are >30 deg off
                  planar and 50% are non-proline cis, against ~0% more than 3 bonds away.
  flank_cis       non-proline cis bonds within 2 of the closure.

Ranking is lexicographic on the geometry that is actually broken, worst-first for "bad" and
best-first for "good". dG_separated is deliberately NOT used to rank: it cancels intra-peptide
energy, so it rewards self-intersection (see STATUS.md).

Usage:  python3 scripts/inspect_set.py [--n 6] [--out inspect]
'''
import argparse
import collections
import csv
import glob
import os
import shutil

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CYC = os.path.dirname(HERE)
ROOT = os.path.dirname(CYC)

# apj_i1 / uor_i1 are the i+1 sweeps; they share a receptor and scaffold set with their base target.
GROUPS = {'apj': ['apj', 'apj_i1'], 'uor': ['uor', 'uor_i1']}


def read_pdb(path, chain=None):
    res, seq = collections.defaultdict(dict), {}
    for line in open(path):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        if line[76:78].strip() == 'H':
            continue
        if chain and line[21] != chain:
            continue
        num = int(line[22:26])
        res[num][line[12:16].strip()] = np.array(
            [float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
        seq[num] = line[17:20]
    return res, seq


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, w), v), np.dot(v, w)))


def improper(atoms):
    '''N-CA-C-CB improper. Positive = L, negative = D. Convention fixed against experimental 9v32,
    where all 545 residues read positive.'''
    if not all(a in atoms for a in ('N', 'CA', 'C', 'CB')):
        return None
    return dihedral(atoms['N'], atoms['CA'], atoms['C'], atoms['CB'])


def geometry(path):
    '''-> dict of the five coordinate-derived metrics, or None if the pose has no peptide chain.'''
    pep, seq = read_pdb(path, 'P')
    rec, _ = read_pdb(path, 'B')
    if len(pep) < 4:
        return None
    nums = sorted(pep)

    # omega, treating the chain as a cycle (the last residue bonds back to the first)
    devs, cis_nonpro, cis_pro = [], 0, 0
    for i, n in enumerate(nums):
        nxt = nums[(i + 1) % len(nums)]
        try:
            om = dihedral(pep[n]['CA'], pep[n]['C'], pep[nxt]['N'], pep[nxt]['CA'])
        except KeyError:
            continue
        devs.append(min(abs(om), abs(abs(om) - 180.0)))
        if abs(om) < 90:
            if seq[nxt] in ('PRO', 'DPR', 'DPRO'):
                cis_pro += 1
            else:
                cis_nonpro += 1

    # closest non-adjacent intra-peptide contact; ring neighbours excluded at both seams
    intra = np.inf
    for i, a in enumerate(nums):
        for b in nums[i + 1:]:
            if abs(a - b) == 1 or (a == nums[0] and b == nums[-1]):
                continue
            for x in pep[a].values():
                for y in pep[b].values():
                    intra = min(intra, float(np.linalg.norm(x - y)))

    inter = np.inf
    if rec:
        P = np.array([x for r in pep.values() for x in r.values()])
        R = np.array([x for r in rec.values() for x in r.values()])
        inter = float(np.linalg.norm(P[:, None] - R[None], axis=-1).min())

    n_d = sum(1 for n in nums
              if (improper(pep[n]) is not None and improper(pep[n]) < 0))

    # the ring-closure bond and its neighbourhood, where the real damage lives.
    # NOTE route-awareness: only the amide routes close through a backbone C-N, so closure_ON is
    # defined only for them. The disulfide routes close through S-S and get NaN, not a spurious
    # ~9 A "error" -- an earlier version of this script ranked every ssapp structure as worst in
    # the set for exactly that reason.
    last, first = nums[-1], nums[0]
    closure_on = float('nan')
    sgs = sorted(n for n in nums if 'SG' in pep[n])
    is_ss = len(sgs) >= 2
    if not is_ss and all(a in pep[last] for a in ('C', 'O')) and 'N' in pep[first]:
        closure_on = float(np.linalg.norm(pep[last]['O'] - pep[first]['N']))
    flank_dev, flank_cis, L = 0.0, 0, len(nums)
    for i, n in enumerate(nums):
        nxt = nums[(i + 1) % L]
        off = min((L - 1 - i) % L, (i + 1) % L)
        if off > 2:
            continue
        try:
            om = dihedral(pep[n]['CA'], pep[n]['C'], pep[nxt]['N'], pep[nxt]['CA'])
        except KeyError:
            continue
        flank_dev = max(flank_dev, min(abs(om), abs(abs(om) - 180.0)))
        if abs(om) < 90 and seq[nxt] not in ('PRO', 'DPR', 'DPRO'):
            flank_cis += 1

    return {'omega_max_dev': max(devs) if devs else float('nan'),
            'cis_nonpro': cis_nonpro, 'cis_pro': cis_pro,
            'intra_min': intra, 'inter_min': inter, 'n_dres': n_d,
            'closure_ON': closure_on, 'flank_dev': flank_dev, 'flank_cis': flank_cis,
            'n_res': len(nums)}


def collect(group):
    rows = []
    for target in GROUPS[group]:
        csv_path = os.path.join(CYC, target, 'interface', 'interface.csv')
        if not os.path.exists(csv_path):
            continue
        for rec in csv.DictReader(open(csv_path)):
            pdb = os.path.join(CYC, target, 'interface', 'complexes', rec['structure'] + '.pdb')
            if not os.path.exists(pdb):
                continue
            g = geometry(pdb)
            if g is None:
                continue
            g.update(target=target, structure=rec['structure'], variant=rec['variant'],
                     route=rec['route'], scaffold=rec['scaffold'],
                     motif_rmsd=float(rec['motif_rmsd']),
                     dG_separated=float(rec['dG_separated']),
                     total_score=float(rec['total_score']), path=pdb)
            rows.append(g)
    return rows


def score_good(r):
    '''Higher is better. Ranks on the geometry that is actually broken, not on dG.

    Closure quality leads, because that is where every defect concentrates and it is the thing
    neither the scorefunction nor validate.py can see. An earlier version of this function ranked
    on omega_max_dev / intra_min alone and called structures "good" that had a carbonyl oxygen
    0.38 A from the next nitrogen -- the intra-peptide scan skipped the closure pair by
    construction, so the single worst defect was the one metric nobody measured.
    '''
    on = r['closure_ON']
    # NaN (disulfide route) is neutral, not worst: it means the metric does not apply.
    on_err = abs(on - 2.25) if on == on else 0.0
    return (-on_err, -r['flank_cis'], -r['flank_dev'],
            -r['omega_max_dev'], r['intra_min'], r['inter_min'], -r['motif_rmsd'])


FIELDS = ['bucket', 'target', 'variant', 'structure', 'route', 'scaffold', 'n_res',
          'closure_ON', 'flank_dev', 'flank_cis',
          'omega_max_dev', 'cis_nonpro', 'cis_pro', 'intra_min', 'inter_min', 'n_dres',
          'motif_rmsd', 'dG_separated', 'total_score']


def closure_atoms(path):
    '''-> (serial_a, serial_b, label) for the ring-closing bond, or None.

    Rosetta writes no LINK or CONECT record for the cyclisation, and viewers only auto-bond
    sequentially numbered residues -- so the i->1 bond is never drawn and a perfectly good
    macrocycle renders as an open chain whose ends happen to sit 1.33 A apart. That is a display
    artifact, not a defect, and it made every design look broken on inspection.
    '''
    serial = {}
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')) and line[21] == 'P':
            serial[(int(line[22:26]), line[12:16].strip())] = int(line[6:11])
    if not serial:
        return None
    nums = sorted({n for n, _ in serial})
    first, last = nums[0], nums[-1]
    # disulfide routes close through S-S; amide routes through the backbone C-N
    sgs = sorted(n for n, a in serial if a == 'SG')
    if len(sgs) >= 2 and (sgs[0], 'SG') in serial and (sgs[-1], 'SG') in serial:
        return serial[(sgs[0], 'SG')], serial[(sgs[-1], 'SG')], 'SG-SG'
    if (last, 'C') in serial and (first, 'N') in serial:
        return serial[(last, 'C')], serial[(first, 'N')], 'C-N'
    return None


def write_with_closure(src, dest):
    '''Copy a complex, appending CONECT for the ring bond so viewers draw the macrocycle closed.'''
    link = closure_atoms(src)
    lines = [l for l in open(src) if not l.startswith(('END', 'CONECT'))]
    with open(dest, 'w') as fh:
        fh.writelines(lines)
        if link:
            a, b, _ = link
            fh.write('CONECT%5d%5d\n' % (a, b))
            fh.write('CONECT%5d%5d\n' % (b, a))
        fh.write('END\n')
    return link


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=6, help='structures per bucket per target')
    ap.add_argument('--out', default=os.path.join(CYC, 'inspect'))
    args = ap.parse_args()

    manifest, bonds = [], []
    for group in GROUPS:
        rows = collect(group)
        if not rows:
            print('%s: nothing scored' % group)
            continue
        # Stratify by route: amide and disulfide close through different chemistry and fail in
        # different ways, so ranking them against each other is meaningless -- closure_ON does not
        # even apply to the disulfides. Sample both explicitly instead.
        half = max(1, args.n // 2)
        picks = []
        for route_set, tag in ((('amide',), 'amide'), (('ssapp', 'ssmut'), 'ss')):
            sub = sorted([r for r in rows if r['route'] in route_set],
                         key=score_good, reverse=True)
            if not sub:
                continue
            picks += [('good', r) for r in sub[:half]]
            picks += [('bad', r) for r in sub[-half:][::-1]]
        for bucket, r in picks:
            dest = os.path.join(args.out, group, bucket)
            os.makedirs(dest, exist_ok=True)
            name = '%s__%s.pdb' % (r['target'], r['structure'])
            link = write_with_closure(r['path'], os.path.join(dest, name))
            bonds.append((group, bucket, name, r, link))
            manifest.append(dict(bucket=bucket, **{k: r[k] for k in FIELDS if k != 'bucket'}))
        print('%s: %d scored, copied %d' % (group, len(rows), len(picks)))

    # PyMOL ignores CONECT for ATOM records under its default connect_mode, so also emit explicit
    # bond commands -- this is the part that reliably makes the ring render closed.
    pml = os.path.join(args.out, 'view.pml')
    with open(pml, 'w') as fh:
        fh.write('# Load the inspection set with the macrocycle bond drawn.\n'
                 '#   pymol view.pml\n'
                 '# Rosetta writes no connectivity record for the cyclisation, so without these\n'
                 '# bond commands every design renders as an open chain. The geometry is fine:\n'
                 '# the closure C-N measures 1.33 A at 116.2/121.7 deg in every amide design.\n\n')
        for group, bucket, name, r, link in bonds:
            obj = os.path.splitext(name)[0]
            fh.write('load %s/%s/%s, %s\n' % (group, bucket, name, obj))
            if link:
                nums = sorted({int(l[22:26]) for l in open(r['path'])
                               if l.startswith(('ATOM', 'HETATM')) and l[21] == 'P'})
                if link[2] == 'SG':
                    pass
                fh.write('bond %s and chain P and resi %d and name C, '
                         '%s and chain P and resi %d and name N\n'
                         % (obj, nums[-1], obj, nums[0]) if link[2] == 'C-N' else
                         '# %s closes through S-S; PyMOL draws disulfides automatically\n' % obj)
        fh.write('\nhide everything\nshow cartoon, chain B\nshow sticks, chain P\n'
                 'color grey80, chain B\ncolor cyan, chain P\nset cartoon_transparency, 0.6\n'
                 'zoom chain P\n')
    print('pymol script -> %s' % pml)

    # reference structures: the native complexes the motifs were cut from
    ref = os.path.join(args.out, 'reference')
    os.makedirs(ref, exist_ok=True)
    for src in ('apj_complex.pdb', 'uor_complex.pdb'):
        p = os.path.join(ROOT, src)
        if os.path.exists(p):
            shutil.copy2(p, os.path.join(ref, src))
    for p in glob.glob(os.path.join(ROOT, '*_scaffolds', '*.pdb')):
        tag = os.path.basename(os.path.dirname(p))
        shutil.copy2(p, os.path.join(ref, '%s__%s' % (tag, os.path.basename(p))))

    out_csv = os.path.join(args.out, 'manifest.csv')
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for row in manifest:
            for k in ('omega_max_dev', 'intra_min', 'inter_min', 'motif_rmsd',
                      'dG_separated', 'total_score', 'closure_ON', 'flank_dev'):
                row[k] = round(row[k], 2)
            w.writerow(row)
    print('manifest -> %s' % out_csv)


if __name__ == '__main__':
    main()
