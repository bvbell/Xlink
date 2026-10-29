'''python3 tests/test_cyclize.py  (~2 s, no Rosetta needed)

Pins the parts that encode the science rather than the plumbing: the hairpin register, the rules that
decide which routes a scaffold is eligible for, the phi-sign turn-seed rules, and the index arithmetic
in the generated XML. A Rosetta run is NOT exercised here -- use `cyclize smoke` for that.
'''
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REPO = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)

from cyclize import build, config, hairpin, struct, variants   # noqa: E402

FAILS = []


def check(label, got, want):
    if got != want:
        FAILS.append('%s: got %r, want %r' % (label, got, want))
        print('  FAIL %s: got %r, want %r' % (label, got, want))
    else:
        print('  ok   %s' % label)


def approx(label, got, want, tol=0.05):
    if got is None or abs(got - want) > tol:
        FAILS.append('%s: got %s, want %s +/- %s' % (label, got, want, tol))
        print('  FAIL %s: got %s, want %s' % (label, got, want))
    else:
        print('  ok   %s (%.2f)' % (label, got))


def motif(path, chain='A'):
    return struct.chain_residues(struct.read(os.path.join(REPO, path), chain), chain)


def test_register():
    '''Every backbone H-bond in both CDR3s obeys i + j = 209, with no exceptions.'''
    print('register')
    for pdb, chain in [('apj_complex.pdb', 'A'), ('uor_complex.pdb', 'A')]:
        full = [r for r in motif(pdb, chain) if 94 <= r.num <= 115]
        k, n_ok, n_tot = hairpin.register(full)
        check('%s register K' % pdb, k, 209)
        check('%s all H-bonds agree' % pdb, n_ok, n_tot)
    # Every scaffold the user selected is on-register, which forces an even residue count.
    for path in ['apj_scaffolds/short.pdb', 'apj_scaffolds/medium.pdb', 'apj_scaffolds/long.pdb',
                 'uOR_scaffolds/short.pdb', 'uOR_scaffolds/medium.pdb', 'uOR_scaffolds/long.pdb']:
        res = motif(path)
        check('%s termini sum' % path, res[0].num + res[-1].num, 209)
        check('%s even length' % path, len(res) % 2, 0)


def test_staple_eligibility():
    '''Only uOR/short's terminal pair is close enough to staple by mutation (Cb-Cb 4.4 A).'''
    print('staple eligibility')
    want = {'apj_scaffolds/short.pdb': False, 'apj_scaffolds/medium.pdb': False,
            'apj_scaffolds/long.pdb': False, 'uOR_scaffolds/short.pdb': True,
            'uOR_scaffolds/medium.pdb': False, 'uOR_scaffolds/long.pdb': False}
    for path, ok in want.items():
        t = hairpin.terminal_pair(motif(path), 209)
        check('%s staple_ok' % path, t['staple_ok'], ok)
    approx('uOR/short Cb-Cb', hairpin.terminal_pair(motif('uOR_scaffolds/short.pdb'), 209)['cb'], 4.4)


def test_turn_rules():
    '''L-Pro can never sit at i+1; the seed table must be self-consistent.'''
    print('turn seed rules')
    for a, b, turn, _ in variants.TURN_SEEDS:
        check('seed %s-%s valid' % (a, b), variants.check_seed(a, b, turn), True)
    try:
        variants.check_seed('PRO', 'GLY', "II'")
        FAILS.append('L-Pro at i+1 was accepted')
        print('  FAIL L-Pro at i+1 was accepted')
    except variants.VariantError:
        print('  ok   L-Pro at i+1 rejected')
    check("turn_type(+,-)", hairpin.turn_type(60, -80), "II'")
    check("turn_type(+,+)", hairpin.turn_type(60, 90), "I'")
    check("turn_type(-,-) not a hairpin turn", hairpin.turn_type(-60, -30) in hairpin.HAIRPIN_TURNS, False)


def test_offregister_refused():
    '''An off-register motif must be a hard error, not a silently bad design.'''
    print('off-register refusal')
    res = motif('apj_scaffolds/long.pdb')[:-1]     # 99-109, sums to 208
    try:
        variants.enumerate_variants('bad', res, {})
        FAILS.append('off-register motif was accepted')
        print('  FAIL off-register motif was accepted')
    except variants.VariantError as e:
        check('mentions the register', 'register' in str(e), True)


def test_plan_indices():
    '''The generated pose layout and GenKIC loop must match what build.py documents.'''
    print('plan indices')
    res = motif('apj_scaffolds/long.pdb')
    scan = config.read_scan(os.path.join(REPO, 'flexddg/apj_nb/analysis_output/alascan_summary.csv'))
    vs = {v['id']: v for v in variants.enumerate_variants('long', res, scan)}

    amide = build.plan(vs['long__amide__DPRO-GLY'], res)
    check('amide last residue', amide['last'], 14)
    check('amide loop', amide['loop'], [13, 14, 1, 2, 3])
    check('amide loop >= 4 so sampling is not fully determined', len(amide['loop']) >= 4, True)
    check('amide closes C->N', (amide['close']['atom1'], amide['close']['atom2']), ('C', 'N'))
    check('amide seeded residue', build.seeded(vs['long__amide__DPRO-GLY'], amide['n']), 13)

    ss = build.plan(vs['long__ssapp__CYD-CYD'], res)
    check('ssapp motif shifted by the prepend', ss['motif_range'], (2, 13))
    check('ssapp closes SG-SG', (ss['close']['atom1'], ss['close']['atom2']), ('SG', 'SG'))
    check('ssapp loop crosses the disulfide', ss['loop'][:2], [14, 1])
    # The rama perturber aborts on residues lacking a polymer connection; in this topology the ring
    # closes through the S-S, so residues 1 and 14 are polymer-terminal.
    check('ssapp polymer termini excluded', ss['polymer_termini'], {1, 14})
    # The motif offset must reach the manifest: validation pairs design to native positionally, and
    # an offset of 0 here misaligned every pair, so the first apj run scored 0 of 150 ssapp structures.
    check('ssapp motif offset is 1', ss['motif_range'][0] - 1, 1)
    check('amide motif offset is 0', amide['motif_range'][0] - 1, 0)

    x = build.xml(vs['long__amide__DPRO-GLY'], amide, config.DEFAULTS['protocol'])
    check('DPRO is built', 'resname="DPRO"' in x, True)
    check('rama filter off by default', 'rama_prepro_check' in x, False)
    check('relax is MoveMap-restricted', '<MoveMap' in x, True)
    check('motif is pinned', 'AddConstraintsToCurrentConformationMover' in x, True)
    check('amide uses PeptideCyclizeMover', 'PeptideCyclizeMover' in x, True)

    xs = build.xml(vs['long__ssapp__CYD-CYD'], ss, config.DEFAULTS['protocol'])
    check('ssapp uses ForceDisulfides', 'ForceDisulfides' in xs, True)
    check('ssapp prepends', '<Prepend' in xs, True)
    # Without these the sulfurs are never moved off the extended build geometry and ssapp closes 0/50.
    check('ssapp randomizes the S-S torsion',
          '<AddAtoms res1="14" atom1="SG" res2="1" atom2="SG"/>' in xs, True)
    addatoms = [l for l in xs.splitlines() if '<AddAtoms' in l]
    check('ssapp randomizes both Cys side chains',
          sum('atom1="CB"' in l and 'atom2="SG"' in l for l in addatoms), 2)
    check('ssapp perturbs 7 torsions in total', len(addatoms), 7)
    check('ssapp randomizes the terminal Cys backbone',
          '<AddAtoms res1="14" atom1="N" res2="14" atom2="CA"/>' in xs, True)
    check('ssapp has 3 randomize_dihedral perturbers',
          xs.count('effect="randomize_dihedral"'), 3)
    check('amide has none', x.count('effect="randomize_dihedral"'), 0)


def test_ssmut_needs_no_sampling():
    '''uOR/short staples the native pair, so there is no loop and no bond to declare.'''
    print('ssmut')
    res = motif('uOR_scaffolds/short.pdb')
    scan = config.read_scan(os.path.join(REPO, 'flexddg/uor_nb/analysis_output/alascan_summary.csv'))
    vs = {v['id']: v for v in variants.enumerate_variants('short', res, scan)}
    check('ssmut offered for uOR/short', 'short__ssmut__CYS-CYS' in vs, True)
    p = build.plan(vs['short__ssmut__CYS-CYS'], res)
    check('no residues added', p['added'], [])
    check('no loop', p['loop'], [])
    check('no sampling', p['sample'], False)
    # MutateResidue makes a plain CYS with no SG connection, so DeclareBond would fail here;
    # ForceDisulfides applies the CYD patch and forms the bond itself.
    check('bond not declared explicitly', p['declare'], False)
    check('both termini mutated', sorted(n for n, _ in p['mutate']), [1, 8])


def test_cys_names():
    """There is no "DCYD": the D disulfide cysteine must be spelled "DCYS:disulfide"."""
    print('cysteine naming')
    labels = [l for l, _ in variants.CYS_FORMS]
    resnames = [r for _, r in variants.CYS_FORMS]
    check('L form resname', resnames[0], 'CYD')
    check('D form resname', resnames[1], 'DCYS:disulfide')
    check('no bare DCYD is ever emitted as a resname', 'DCYD' in resnames, False)
    # A colon cannot appear in a variant id, because ids become directory names.
    for lab in labels:
        check('label %s is path-safe' % lab, ':' in lab or '/' in lab, False)
    res = motif('apj_scaffolds/long.pdb')
    scan = config.read_scan(os.path.join(REPO, 'flexddg/apj_nb/analysis_output/alascan_summary.csv'))
    for v in variants.enumerate_variants('long', res, scan):
        if ':' in v['id'] or '/' in v['id']:
            FAILS.append('variant id %r is not path-safe' % v['id'])
            print('  FAIL variant id %r is not path-safe' % v['id'])
    print('  ok   every variant id is path-safe')


def test_runs_helper():
    print('misc helpers')
    check('runs splits non-contiguous loops', build.runs([13, 14, 1, 2, 3]), [(1, 3), (13, 14)])
    check('pivots spread across the loop', [x[0] for x in build.pick_pivots([13, 14, 1, 2, 3])],
          [13, 1, 3])
    check('D residues read from HETATM', struct.AA3to1['DPR'], 'p')
    check('ss_torsions groups: bond, sidechains, backbone',
          [len(g) for g in build.ss_torsions(14, 1)], [1, 4, 2])


if __name__ == '__main__':
    for fn in [test_register, test_staple_eligibility, test_turn_rules, test_offregister_refused,
               test_plan_indices, test_ssmut_needs_no_sampling, test_cys_names, test_runs_helper]:
        fn()
    print('\n%d failure(s)' % len(FAILS))
    for f in FAILS:
        print('  ' + f)
    sys.exit(1 if FAILS else 0)
