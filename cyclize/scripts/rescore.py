'''Score closed designs on rama_prepro -- total, and summed over the GenKIC loop.

The naive way to do this is wrong. rama_prepro is a two-body term on consecutive pairs, and for a
head-to-tail ring the wrap pair (1, N) only exists if the pose actually knows it is cyclic. Reading
an amide design straight from its PDB takes phi(1) from the atom-tree root stub instead of the
geometry, so the seam residues -- which are the turn residues -- come out wrong, silently, by an
amount that depends on where the numbering happens to start. PeptideCyclizeMover fixes it.

The disulfide routes take the OPPOSITE rule, and this is the part that is easy to get backwards:

  amide  head-to-tail ring       -> PeptideCyclizeMover REQUIRED
  ssapp  linear, closed by S-S   -> read as-is; the mover would force an amide across a 4.1 A gap
  ssmut  linear, closed by S-S   -> read as-is, same reason

There is also a live bug to steer around (see "Known caveats" in the README): build.xml strips
LOWER_TERMINUS_VARIANT from residue 1 for every route, including ssapp, where nothing is ever bonded
to that N. That made Rosetta score residue 1's rama at phi = 0 -- a forbidden region -- and invent
about +10 REU during every ssapp run. Rescoring must NOT reproduce it: we read the deposited PDB
as-is and let Rosetta apply the normal terminus patch, which correctly gives residue 1 no rama at all.

That last point makes the loop sums route-dependent in a way worth remembering when comparing them:
on a disulfide route the two chain-terminal residues have no defined rama_prepro, and on ssapp those
are exactly the two Cys the route adds. `n_loop_scored` records how many of the loop's residues
actually carried a value, so a short sum is visible rather than silently flattering.

Writes <target>/rescore/rescore.csv. Nothing is written into output/<vid>/ -- runner.state counts
*.pdb there to decide whether a variant is done.
'''
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ROSETTA = 'rosetta_scripts.default.linuxgccrelease'
BREAKDOWN = 'residue_energy_breakdown.default.linuxgccrelease'
CYCLIC_ROUTES = ('amide',)          # everything else closes through a disulfide and reads as linear


def xml(cyclize_it):
    '''RosettaScripts protocol for one route class. `cyclize_it` is the whole difference.'''
    mover = '    <PeptideCyclizeMover name="cyc" residue_selector="pep"/>\n' if cyclize_it else ''
    step = '    <Add mover="cyc"/>\n' if cyclize_it else ''
    return ('<ROSETTASCRIPTS>\n'
            '  <SCOREFXNS><ScoreFunction name="sfxn" weights="ref2015"/></SCOREFXNS>\n'
            '  <RESIDUE_SELECTORS><Chain name="pep" chains="A"/></RESIDUE_SELECTORS>\n'
            '  <SIMPLE_METRICS>\n'
            '    <TotalEnergyMetric name="rpp_tot" scoretype="rama_prepro" scorefxn="sfxn"/>\n'
            '  </SIMPLE_METRICS>\n'
            '  <MOVERS>\n' + mover +
            '    <RunSimpleMetrics name="report" metrics="rpp_tot" prefix="rpp_"/>\n'
            '  </MOVERS>\n'
            '  <PROTOCOLS>\n' + step +
            '    <Add mover="report"/>\n'
            '  </PROTOCOLS>\n'
            '  <OUTPUT scorefxn="sfxn"/>\n'
            '</ROSETTASCRIPTS>\n')


def run(cmd, cwd):
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit('%s failed (exit %d):\n%s' % (cmd[0], out.returncode, out.stderr[-3000:]))
    return out


def parse_breakdown(path):
    '''-> {tag: {(resi1, resi2): rama_prepro}}, two-body rows only.

    residue_energy_breakdown writes a whitespace table whose header names resi1/resi2 and one column
    per score term. The structure tag is `pose_id`; the trailing `description` is per-ROW
    ("<tag>_<resi1>_<resi2>") and keying on it silently yields one pair per structure. One-body rows
    carry "--" for resi2 and are skipped, rama_prepro being two-body.
    '''
    out, cols = {}, None
    for line in open(path):
        f = line.split()
        if not f or f[0] != 'SCORE:':
            continue
        if cols is None and 'resi1' in f:
            cols = f[1:]
            for need in ('pose_id', 'resi1', 'resi2', 'rama_prepro'):
                if need not in cols:
                    sys.exit('%s: breakdown table has no %r column' % (path, need))
            continue
        if cols is None:
            continue
        row = dict(zip(cols, f[1:]))
        if row.get('resi2') in (None, '--'):
            continue
        out.setdefault(row['pose_id'], {})[(int(row['resi1']), int(row['resi2']))] = \
            float(row['rama_prepro'])
    if cols is None:
        sys.exit('%s: no header row with resi1; residue_energy_breakdown output not understood' % path)
    return out


def owner(pair, n):
    '''Which residue a consecutive-pair rama_prepro term belongs to.

    The term for residue i lives on (i, i+1). The cyclic wrap row is (1, N) and belongs to N, not 1 --
    getting this backwards moves the turn's largest contribution onto a pinned motif residue.
    '''
    a, b = pair
    if {a, b} == {1, n} and n > 2:
        return n
    return min(a, b)


def collect(target, use_all):
    '''-> [(pdb path, variant, route)] for one target.'''
    root = os.path.join(ROOT, target)
    manifest = os.path.join(root, 'representatives', 'representatives.csv')
    if not use_all:
        if not os.path.exists(manifest):
            sys.exit('%s: no representatives yet; run scripts/representatives.py first (or pass '
                     '--all to score every passing structure)' % target)
        rows = list(csv.DictReader(open(manifest)))
        return [(os.path.join(root, 'representatives', r['pdb']), r['variant'], r['route'])
                for r in rows]
    rows = [r for r in csv.DictReader(open(os.path.join(root, 'analysis_output/results.csv')))
            if r['pass'] == 'True']
    return [(os.path.join(root, 'output', r['variant'], r['structure']), r['variant'], r['route'])
            for r in rows]


def score_target(target, use_all=False, keep=False):
    root = os.path.join(ROOT, target)
    work = os.path.join(root, 'rescore')
    os.makedirs(work, exist_ok=True)
    jobs = collect(target, use_all)
    manifests = {}

    rows = []
    for cyclic, label in ((True, 'cyclic'), (False, 'linear')):
        group = [j for j in jobs if (j[2] in CYCLIC_ROUTES) == cyclic]
        if not group:
            continue
        print('%s/%s: %d structures (%s)' % (target, label, len(group),
                                             'PeptideCyclizeMover' if cyclic else 'read as-is'))
        open(os.path.join(work, '%s.xml' % label), 'w').write(xml(cyclic))
        listfile = os.path.join(work, '%s.list' % label)
        open(listfile, 'w').write(''.join(p + '\n' for p, _, _ in group))
        silent = '%s.silent' % label
        run([ROSETTA, '-l', listfile, '-parser:protocol', '%s.xml' % label,
             '-out:file:silent', silent, '-out:file:silent_struct_type', 'binary',
             '-out:file:scorefile', '%s.sc' % label, '-no_nstruct_label',
             '-overwrite', '-mute', 'all'], cwd=work)
        # Hand the CYCLIZED pose on as a binary silent: a PDB does not round-trip the declared bond
        # and would put the breakdown straight back into the naive state.
        reb = '%s.reb' % label
        run([BREAKDOWN, '-in:file:silent', silent, '-in:file:silent_struct_type', 'binary',
             '-score:weights', 'ref2015', '-out:file:silent', reb,
             '-overwrite', '-mute', 'all'], cwd=work)
        pairs_by_tag = parse_breakdown(os.path.join(work, reb))

        for path, variant, route in group:
            tag = os.path.basename(path)[:-4]
            pairs = pairs_by_tag.get(tag)
            if pairs is None:
                print('  WARNING: no breakdown rows for %s' % tag)
                continue
            if variant not in manifests:
                manifests[variant] = json.load(
                    open(os.path.join(root, 'inputs', variant, 'input.json')))
            man = manifests[variant]
            n = man['rosetta_last']
            loop = set(man['rosetta_loop'])
            per = {}
            for pair, value in pairs.items():
                per[owner(pair, n)] = per.get(owner(pair, n), 0.0) + value
            scored = sorted(loop & set(per))
            rows.append({
                'target': target, 'structure': tag, 'variant': variant, 'route': route,
                'scaffold': variant.split('__')[0], 'turn': variant.split('__')[2],
                'cyclized': int(route in CYCLIC_ROUTES),
                'n_res': n,
                'rama_total': round(sum(per.values()), 3),
                # Blank, never 0.0, when nothing scored: ssmut has an empty GenKIC loop, and a zero
                # would sort to the top of a ranking as though it were the best turn in the set.
                'rama_loop': round(sum(per[i] for i in scored), 3) if scored else '',
                'rama_loop_per_res': round(sum(per[i] for i in scored) / len(scored), 3) if scored else '',
                'n_loop': len(loop), 'n_loop_scored': len(scored),
                'loop_unscored': ' '.join(str(i) for i in sorted(loop - set(per))),
            })

    if not keep:
        for f in os.listdir(work):
            if f.endswith(('.silent', '.reb', '.list')):
                os.remove(os.path.join(work, f))
    out = os.path.join(work, 'rescore.csv')
    rows.sort(key=lambda r: (r["scaffold"], r["rama_loop"] if r["rama_loop"] != "" else float("inf")))
    with open(out, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('  -> %s (%d rows)' % (out, len(rows)))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('targets', nargs='*', default=['apj', 'uor'])
    ap.add_argument('--all', action='store_true',
                    help='score every passing structure, not just the cluster representatives')
    ap.add_argument('--keep', action='store_true', help='keep the silent/breakdown intermediates')
    args = ap.parse_args(argv)

    for tool in (ROSETTA, BREAKDOWN):
        if not shutil.which(tool):
            sys.exit('%s not on PATH; ml biology rosetta/3.14 first' % tool)

    for target in args.targets:
        rows = score_target(target, args.all, args.keep)
        best = sorted([r for r in rows if r["rama_loop"] != ""], key=lambda r: r["rama_loop"])[:5]
        for r in best:
            print('      %-34s loop %7.2f over %d/%d res   total %8.2f'
                  % (r['variant'], r['rama_loop'], r['n_loop_scored'], r['n_loop'], r['rama_total']))


if __name__ == '__main__':
    main()
