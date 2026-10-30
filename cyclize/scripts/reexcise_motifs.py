#!/usr/bin/env python3
'''Re-excise the CDR3 scaffolds from relaxed input complexes, preserving residue numbering.

Written 2026-10-30 as step 2 of the rebuild (step 1 is scripts/relax_inputs.sbatch).

Why: the previous design set was invalid partly because the motif source and the alignment target
had drifted apart. On uOR they were actively inconsistent -- residues 102-107 are D-amino acids in
`uor_complex.pdb` but Rosetta rebuilds them as L on read, so designs were built L and then
superimposed onto D coordinates. Taking both the scaffold AND the alignment target from the SAME
relaxed pose removes that mismatch by construction.

The residue ranges are read from the EXISTING scaffolds rather than hard-coded, so the new
scaffolds cover exactly the same positions as the old ones and the `i + j = 209` register check
still applies unchanged:

    apj  short 101-108   medium 100-109   long  99-110
    uOR  short 101-108   medium 100-109   long  98-111

Picks the relaxed decoy with the lowest total_score whose motif region stays closest to the input
(so the relax regularises geometry without walking the binding conformation away).

Usage:
    python3 scripts/reexcise_motifs.py                 # report only
    python3 scripts/reexcise_motifs.py --write         # write scaffolds + relaxed complexes
'''
import argparse
import collections
import glob
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CYC = os.path.dirname(HERE)
ROOT = os.path.dirname(CYC)
RELAX = os.path.join(HERE, 'relax_inputs')

TARGETS = {
    'apj': {'input': os.path.join(ROOT, 'apj_complex.pdb'),
            'scaffold_dir': os.path.join(ROOT, 'apj_scaffolds'),
            'out_complex': os.path.join(ROOT, 'apj_complex_relaxed.pdb'),
            'out_scaffold_dir': os.path.join(ROOT, 'apj_scaffolds_relaxed')},
    'uor': {'input': os.path.join(ROOT, 'uor_complex.pdb'),
            'scaffold_dir': os.path.join(ROOT, 'uOR_scaffolds'),
            'out_complex': os.path.join(ROOT, 'uor_complex_relaxed.pdb'),
            'out_scaffold_dir': os.path.join(ROOT, 'uOR_scaffolds_relaxed')},
}


def residues_of(path, chain=None):
    out = collections.defaultdict(dict)
    for line in open(path):
        if line.startswith(('ATOM', 'HETATM')) and line[76:78].strip() != 'H':
            if chain and line[21] != chain:
                continue
            out[int(line[22:26])][line[12:16].strip()] = np.array(
                [float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
    return out


def scaffold_ranges(scaffold_dir):
    '''-> {name: sorted residue numbers}, taken from the existing scaffolds.'''
    out = {}
    for p in sorted(glob.glob(os.path.join(scaffold_dir, '*.pdb'))):
        name = os.path.splitext(os.path.basename(p))[0]
        out[name] = sorted(residues_of(p))
    return out


def ca_rmsd(a, b, nums):
    pairs = [(a[n]['CA'], b[n]['CA']) for n in nums
             if n in a and n in b and 'CA' in a[n] and 'CA' in b[n]]
    if not pairs:
        return float('nan')
    return float(np.sqrt(np.mean([np.sum((x - y) ** 2) for x, y in pairs])))


def read_scores(sc):
    rows = []
    if not os.path.exists(sc):
        return rows
    lines = [l for l in open(sc) if l.startswith('SCORE:')]
    if len(lines) < 2:
        return rows
    hdr = lines[0].split()[1:]
    for l in lines[1:]:
        v = l.split()[1:]
        if len(v) != len(hdr):
            continue
        d = dict(zip(hdr, v))
        try:
            d['total_score'] = float(d['total_score'])
        except (KeyError, ValueError):
            continue
        rows.append(d)
    return rows


def pick(target, cfg, motif_nums):
    decoys = sorted(glob.glob(os.path.join(RELAX, target, '*_000*.pdb')))
    if not decoys:
        return None, []
    scores = {d['description']: d['total_score'] for d in read_scores(
        os.path.join(RELAX, target, 'relax.sc')) if 'description' in d}
    native = residues_of(cfg['input'], 'A')
    rows = []
    for p in decoys:
        tag = os.path.splitext(os.path.basename(p))[0]
        res = residues_of(p, 'A')
        rows.append({'path': p, 'tag': tag,
                     'score': scores.get(tag, float('nan')),
                     'motif_rmsd': ca_rmsd(native, res, motif_nums)})
    # lowest score among decoys whose motif stayed within 0.75 A; fall back to best motif fidelity
    ok = [r for r in rows if r['motif_rmsd'] == r['motif_rmsd'] and r['motif_rmsd'] <= 0.75]
    best = min(ok, key=lambda r: r['score']) if ok else min(rows, key=lambda r: r['motif_rmsd'])
    return best, rows


def write_scaffold(src_complex, nums, dest):
    '''Copy chain A residues `nums` verbatim out of a relaxed complex, numbering preserved.'''
    kept = []
    for line in open(src_complex):
        if line.startswith(('ATOM', 'HETATM')) and line[21] == 'A':
            if int(line[22:26]) in nums:
                kept.append(line)
    with open(dest, 'w') as fh:
        fh.writelines(kept)
        fh.write('TER\nEND\n')
    return len(kept)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()

    # Chained afterany on the relax, not afterok: the relax is paced to finish ~10 decoys in about
    # the time it is allotted, so a timeout on the last one is likely and is harmless -- 9 decoys
    # are as good as 10 for picking a relaxed pose. But if NOTHING was produced we must fail loudly
    # so the geometry test (afterok on this job) does not run against missing inputs.
    wrote = 0
    for target, cfg in TARGETS.items():
        ranges = scaffold_ranges(cfg['scaffold_dir'])
        allnums = sorted({n for v in ranges.values() for n in v})
        best, rows = pick(target, cfg, allnums)
        print('\n=== %s ===' % target)
        if best is None:
            print('  no relaxed decoys yet in %s' % os.path.join(RELAX, target))
            continue
        for r in sorted(rows, key=lambda r: r['score']):
            mark = '  <-- chosen' if r['tag'] == best['tag'] else ''
            print('  %-40s score %10.1f  motif CA RMSD %5.2f%s'
                  % (r['tag'], r['score'], r['motif_rmsd'], mark))
        if not args.write:
            continue
        shutil.copy2(best['path'], cfg['out_complex'])
        os.makedirs(cfg['out_scaffold_dir'], exist_ok=True)
        for name, nums in ranges.items():
            dest = os.path.join(cfg['out_scaffold_dir'], name + '.pdb')
            n = write_scaffold(best['path'], set(nums), dest)
            print('  wrote %-44s %2d residues, %3d atoms'
                  % (os.path.relpath(dest, ROOT), len(nums), n))
        print('  wrote %s' % os.path.relpath(cfg['out_complex'], ROOT))
        wrote += 1

    if args.write and wrote < len(TARGETS):
        sys.exit('ERROR: wrote %d of %d targets -- relax produced no usable decoys for the rest'
                 % (wrote, len(TARGETS)))


if __name__ == '__main__':
    main()
