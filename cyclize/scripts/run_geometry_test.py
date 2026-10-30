#!/usr/bin/env python3
'''Build and score a handful of variants under one test arm, then report closure geometry.

Driven by scripts/geometry_test.sbatch; see that file for the experimental design. This runs the
real `build.xml` path so the test measures the production protocol, not a reimplementation of it.

Arms:
Validity requires NO cis peptide bonds at all (user's rule, 2026-10-30) -- not merely no
non-proline cis. Both input nanobodies are 100% trans, so this forbids nothing that was native.

  base    relaxed input, omega unconstrained   (control -- expected to reproduce the failure)
  omega   omega_trans=True                     (pin every sampled peptide bond trans)
  loop7   omega_trans=True, target_loop=7      (also release two more cold terminal motif
                                                residues; identity is preserved, nothing is
                                                designed, and still only 2 residues are added)

Both closure chemistries are tested, because both are independently broken in the archived set:

  amide      closure C-N is fine but the carbonyl O lands on the next N in 57% of structures,
             and the two bonds flanking the closure are 47% non-planar / 50% non-proline cis
  ssapp      S-S distance median 2.28 A against an ideal 2.05, and as high as 5.82 A -- which is
             not a bond at all; |chi3| (CB-SG-SG-CB) outside the 60-120 deg gauche window in
             43 of 74 structures
  ssmut      the single archived structure has S-S 1.70 A (compressed) and chi3 164 deg (trans),
             which is the pinned-motif strain already noted in ../README.md

omega_trans pins peptide bonds only and deliberately skips the S-S, so for the disulfide routes it
fixes the FLANKING bonds while the S-S geometry itself is governed by the ForceDisulfides +
FastRelax step. Both are measured here.

Writes <out>/<target>__<scaffold>__<arm>__<routetag>/{*.pdb, result.csv}. result.csv carries one
row per structure with the geometry that decides validity.
'''
import argparse
import collections
import csv
import glob
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CYC = os.path.dirname(HERE)
ROOT = os.path.dirname(CYC)
sys.path.insert(0, CYC)

from cyclize import build, config, struct, variants  # noqa: E402

# the relaxed inputs produced by relax_inputs.sbatch + reexcise_motifs.py
RELAXED = {
    'apj': {'complex': os.path.join(ROOT, 'apj_complex_relaxed.pdb'),
            'scaffolds': os.path.join(ROOT, 'apj_scaffolds_relaxed'),
            'scan': os.path.join(ROOT, 'flexddg/apj_nb/analysis_output/alascan_summary.csv')},
    'uor': {'complex': os.path.join(ROOT, 'uor_complex_relaxed.pdb'),
            'scaffolds': os.path.join(ROOT, 'uOR_scaffolds_relaxed'),
            'scan': os.path.join(ROOT, 'flexddg/uor_nb/analysis_output/alascan_summary.csv')},
}
ARMS = {
    'base':  {'omega_trans': False},
    'omega': {'omega_trans': True},
    'loop7': {'omega_trans': True, 'target_loop': 7},
}


def load_scan(path):
    scan = {}
    if not os.path.exists(path):
        return scan
    for row in csv.DictReader(open(path)):
        case = row.get('case', '')
        try:
            scan[int(case[2:-1])] = float(row['ddG_GAM'])
        except (ValueError, KeyError):
            continue
    return scan


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, w), v), np.dot(v, w)))


def geometry(path):
    res, seq = collections.defaultdict(dict), {}
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')) and line[76:78].strip() != 'H':
            n = int(line[22:26])
            res[n][line[12:16].strip()] = np.array(
                [float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
            seq[n] = line[17:20]
    nums = sorted(res)
    if len(nums) < 4:
        return None
    worst_dev, cis_nonpro, n_cis, L = 0.0, 0, 0, len(nums)
    for i, n in enumerate(nums):
        nxt = nums[(i + 1) % L]
        try:
            om = dihedral(res[n]['CA'], res[n]['C'], res[nxt]['N'], res[nxt]['CA'])
        except KeyError:
            continue
        worst_dev = max(worst_dev, min(abs(om), abs(abs(om) - 180.0)))
        # No cis bonds AT ALL -- cis-Pro is legitimate in nature but a liability in a synthesised
        # macrocycle, and both input nanobodies are 100% trans, so nothing native is lost.
        # cis_nonpro is still reported separately because it is what the archived set was measured
        # against, but validity below uses the total.
        if abs(om) < 90:
            n_cis += 1
            if seq[nxt] not in ('PRO', 'DPR', 'DPRO'):
                cis_nonpro += 1
    # amide closure: the carbonyl O must sit ~2.25 A from the next N (O=C-N ~123 deg).
    # disulfide closure: SG-SG ~2.05 A with |CB-SG-SG-CB| in the 60-120 deg gauche window.
    last, first = nums[-1], nums[0]
    sgs = sorted(n for n in nums if 'SG' in res[n])
    is_ss = len(sgs) >= 2
    on, ss_dist, ss_chi3 = float('nan'), float('nan'), float('nan')
    if is_ss:
        a, b = sgs[0], sgs[-1]
        ss_dist = float(np.linalg.norm(res[a]['SG'] - res[b]['SG']))
        if 'CB' in res[a] and 'CB' in res[b]:
            ss_chi3 = float(dihedral(res[a]['CB'], res[a]['SG'], res[b]['SG'], res[b]['CB']))
    elif all(x in res[last] for x in ('C', 'O')) and 'N' in res[first]:
        on = float(np.linalg.norm(res[last]['O'] - res[first]['N']))
    intra = np.inf
    for i, a in enumerate(nums):
        for b in nums[i + 1:]:
            if abs(a - b) == 1 or (a == nums[0] and b == nums[-1]):
                continue
            for x in res[a].values():
                for y in res[b].values():
                    intra = min(intra, float(np.linalg.norm(x - y)))
    return {'omega_max_dev': round(worst_dev, 2), 'cis_nonpro': cis_nonpro, 'n_cis': n_cis,
            'closure_ON': round(on, 3), 'ss_dist': round(ss_dist, 3),
            'ss_chi3': round(ss_chi3, 1), 'is_ss': int(is_ss),
            'intra_min': round(intra, 2), 'n_res': L}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', required=True, choices=sorted(RELAXED))
    ap.add_argument('--scaffold', required=True)
    ap.add_argument('--arm', required=True, choices=sorted(ARMS))
    ap.add_argument('--nstruct', type=int, default=5)
    ap.add_argument('--out', required=True)
    ap.add_argument('--routes', default='amide',
                    help='comma list: amide and/or ssapp,ssmut')
    ap.add_argument('--seeds', default='DPRO-GLY,ASN-GLY,GLY-GLY',
                    help='amide turn seeds to test; ignored for the disulfide routes, whose '
                         '"seed" field is the Cys pair rather than a turn')
    ap.add_argument('--max-variants', type=int, default=3)
    a = ap.parse_args()

    cfg = RELAXED[a.target]
    for p in (cfg['complex'], os.path.join(cfg['scaffolds'], a.scaffold + '.pdb')):
        if not os.path.exists(p):
            sys.exit('missing %s -- run relax_inputs.sbatch then reexcise_motifs.py --write' % p)

    res = struct.chain_residues(struct.read(cfg['complex']), 'A')
    keep = {r.num for r in struct.chain_residues(
        struct.read(os.path.join(cfg['scaffolds'], a.scaffold + '.pdb')), 'A')}
    motif = [r for r in res if r.num in keep]
    scan = load_scan(cfg['scan'])

    proto = dict(config.DEFAULTS['protocol'], nstruct=a.nstruct, **ARMS[a.arm])
    kw = {}
    if 'target_loop' in ARMS[a.arm]:
        kw['target_loop'] = ARMS[a.arm]['target_loop']
    vs = variants.enumerate_variants(a.scaffold, motif, scan, **kw)
    routes = [r for r in a.routes.split(',') if r]
    vs = [v for v in vs if v['route'] in routes]
    if 'amide' in routes:
        # the amide id's trailing field is the turn seed, so select on it
        want = set(a.seeds.split(','))
        vs = ([v for v in vs if v['route'] == 'amide' and v['id'].split('__')[-1] in want]
              + [v for v in vs if v['route'] != 'amide'])
    # disulfide variants are few per scaffold (and ssmut exists only where Cb-Cb is in the
    # 3.4-4.6 A window -- uOR/short alone), so take the first few of each route rather than
    # filtering on a seed name that means something different for them.
    per_route, picked = collections.Counter(), []
    for v in vs:
        if per_route[v['route']] >= a.max_variants:
            continue
        per_route[v['route']] += 1
        picked.append(v)
    vs = picked
    if not vs:
        print('no variants for routes=%s scaffold=%s target=%s' % (a.routes, a.scaffold, a.target))
        return

    tag_routes = 'amide' if routes == ['amide'] else 'ss'
    work = os.path.join(a.out, '%s__%s__%s__%s' % (a.target, a.scaffold, a.arm, tag_routes))
    os.makedirs(work, exist_ok=True)
    rows = []
    for v in vs:
        plan = build.plan(v, motif)
        xml = build.xml(v, plan, proto)
        tag = v['id'].replace('/', '_')
        d = os.path.join(work, tag)
        os.makedirs(d, exist_ok=True)
        xp = os.path.join(d, 'build.xml')
        open(xp, 'w').write(xml)
        cmd = ['rosetta_scripts.default.linuxgccrelease',
               '-parser:protocol', xp, '-s', cfg['complex'],
               '-nstruct', str(a.nstruct), '-out:path:pdb', d,
               '-out:file:scorefile', os.path.join(d, 'score.sc'),
               '-overwrite', '-ignore_unrecognized_res', '-mute', 'all']
        print('  %s loop=%s' % (tag, plan['loop']), flush=True)
        subprocess.run(cmd, cwd=d, check=False)
        for pdb in sorted(glob.glob(os.path.join(d, '*.pdb'))):
            g = geometry(pdb)
            if g is None:
                continue
            g.update(target=a.target, scaffold=a.scaffold, arm=a.arm, variant=v['id'],
                     route=v['route'], loop=len(plan['loop']),
                     structure=os.path.basename(pdb))
            # Validity is route-aware: the two chemistries close through different bonds, so the
            # closure criterion differs. Everything else (flanking planarity, no non-proline cis)
            # applies to both -- omega_trans pins peptide bonds for the disulfide routes too.
            backbone_ok = g['n_cis'] == 0 and g['omega_max_dev'] <= 30
            if g['is_ss']:
                closure_ok = (1.95 <= g['ss_dist'] <= 2.15
                              and 60 <= abs(g['ss_chi3']) <= 120)
            else:
                closure_ok = 2.15 < g['closure_ON'] < 2.35
            g['valid'] = int(backbone_ok and closure_ok)
            rows.append(g)

    if rows:
        fields = ['target', 'scaffold', 'arm', 'route', 'variant', 'loop', 'structure',
                  'closure_ON', 'ss_dist', 'ss_chi3', 'omega_max_dev', 'n_cis',
                  'cis_nonpro', 'intra_min', 'n_res', 'valid']
        with open(os.path.join(work, 'result.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({k: r[k] for k in fields})
        ok = sum(r['valid'] for r in rows)
        print('\n%s %s %s [%s]: %d/%d structures valid'
              % (a.target, a.scaffold, a.arm, a.routes, ok, len(rows)))
        by = collections.Counter()
        for r in rows:
            by[(r['route'], r['valid'])] += 1
        for route in sorted({r['route'] for r in rows}):
            print('   %-6s %d/%d valid' % (route, by[(route, 1)],
                                           by[(route, 1)] + by[(route, 0)]))
    else:
        print('\n%s %s %s [%s]: no structures produced'
              % (a.target, a.scaffold, a.arm, a.routes))


if __name__ == '__main__':
    main()
