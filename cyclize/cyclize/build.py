'''Generate the per-variant RosettaScripts XML and input manifest.

Rosetta renumbers an input PDB to 1..N, so everything below works in Rosetta indices and records the
PDB-index mapping in input.json for the validation step to undo.

Pose layout for a motif of N residues:

    amide / ssmut          1 .. N      the motif, native coordinates
                           N+1, N+2    the two added turn residues (amide only)
    ssapp                  1           Cys prepended before the motif N-terminus
                           2 .. N+1    the motif
                           N+2         Cys appended after the motif C-terminus

Prepending renumbers the pose, so ssapp does it first and plan() shifts every index by one.

The GenKIC loop is [released motif C-term] + [added] + [released motif N-term], ordered so that it
spans the closing bond, with the first non-released motif residue at either end acting as a fixed
anchor. How many cold motif residues get released is decided by variants.pick_release, which targets
more than GenKIC's bare 3-residue minimum -- see the TARGET_LOOP note there for why 3 does not work.
'''
import json, os

from . import struct

PEPBOND = {'length': 1.328685, 'angle1': 116.2, 'angle2': 121.7}   # radians in the source, degrees here
SSBOND = {'length': 2.05, 'angle1': 103.0, 'angle2': 103.0}
# rama_prepro_check is the single biggest throttle on closure here. Measured on apj/long, 4 structures
# x 5000 attempts: filter off 4/4 closed, cutoff 5.0 gave 1/4, cutoff 1.0 on the seeded residue gave
# 0/4, and moving the pivot off the seeded residue gave 0/4 at any cutoff. So GenKIC is left to do only
# what it is good at -- closing the ring -- and backbone quality is judged afterwards: validate.py
# checks the seeded residue's phi sign, and the score function sees rama_prepro anyway. Set
# protocol.rama_cutoff to a number to turn the in-sampling filter back on.
RAMA_CUTOFF = None


def seeded(variant, n):
    '''The chirality-critical residue: i+1 of a seeded turn, the first appended residue.

    GenKIC solves the torsions flanking its pivots, so a pivot on this residue overrides the seed --
    measured on apj/long, D-Pro as pivot 1 came out at phi = -55 and +145, which its ring cannot adopt.
    Keeping a pivot OFF it instead starves the closure of leverage (0/6 closed), so the pivots stay
    spread and the seed is enforced the idiomatic way: a tight rama_prepro_check on this residue
    rejects wrong-handed solutions during sampling, and validate.py re-checks the phi sign after.
    '''
    if variant['route'] == 'amide' and variant['added']:
        return n + 1
    return None


def pick_pivots(loop):
    '''Three (residue, atom) pivots spread across `loop`.

    Pivot atoms must be at least 3 apart in the atom chain; consecutive CA atoms are exactly 3 apart,
    so any three distinct loop residues work, and spreading them gives the closure the most leverage.
    The closed bond must lie between pivot 1 and pivot 3, which holds since it falls inside the loop.
    '''
    return [(loop[0], 'CA'), (loop[len(loop) // 2], 'CA'), (loop[-1], 'CA')]


def ss_torsions(upper, lower):
    """Torsions a disulfide closure has to sample, as (res, atom) pairs naming each central bond.

    Mirrors simple_cycpep_predict's terminal_disulfide setup: the S-S torsion itself, both Cys side
    chains marching CA-CB-SG, and what backbone each Cys has. `upper` is the appended Cys (polymer
    upper terminus, so only phi is defined); `lower` is the prepended one (only psi is defined).
    """
    return [
        [((upper, 'SG'), (lower, 'SG'))],
        [((upper, 'CB'), (upper, 'SG')), ((upper, 'CA'), (upper, 'CB')),
         ((lower, 'CB'), (lower, 'SG')), ((lower, 'CA'), (lower, 'CB'))],
        [((upper, 'N'), (upper, 'CA')), ((lower, 'CA'), (lower, 'C'))],
    ]


def _indices(residues):
    '''-> (pdb_to_rosetta, rosetta_to_pdb) for a motif Rosetta will renumber 1..N.'''
    pdb2r = {r.num: i + 1 for i, r in enumerate(residues)}
    return pdb2r, {v: k for k, v in pdb2r.items()}


def plan(variant, residues):
    '''Resolve a variant into concrete Rosetta indices. -> dict used by xml() and the manifest.'''
    pdb2r, r2pdb = _indices(residues)
    n = len(residues)
    added = list(variant['added'])
    last = n + len(added)
    if variant['route'] == 'ssmut':
        # Nothing is appended and nothing is sampled: the native terminal pair already sits in the
        # disulfide window, so we mutate, declare the bond and minimise.
        return {'n': n, 'added': [], 'last': n, 'loop': [], 'pivots': [], 'release': None,
                'close': {'res1': n, 'atom1': 'SG', 'res2': 1, 'atom2': 'SG',
                          'prior': 'CB', 'following': 'CB', 'geom': SSBOND},
                'mutate': [(pdb2r[num], name) for num, name in variant['mutate']],
                'pdb2r': pdb2r, 'r2pdb': r2pdb, 'sample': False,
                'prepend': None, 'append': None, 'motif_range': (1, n),
                'polymer_termini': {1, n}, 'declare': False, 'ss_torsions': []}
    if variant['route'] == 'ssapp':
        # One Cys beyond each terminus, bonded to each other. Prepending renumbers the pose, so it is
        # done first and everything downstream shifts by one:
        #     1        Cys prepended before the motif N-terminus
        #     2..N+1   the motif
        #     N+2      Cys appended after the motif C-terminus
        # The chain to close runs motif-Cterm -> Cys(N+2) -> S-S -> Cys(1) -> motif-Nterm, so the loop
        # crosses the disulfide rather than a peptide bond.
        pdb2r = {num: i + 2 for i, num in enumerate(r.num for r in residues)}
        r2pdb = {v: k for k, v in pdb2r.items()}
        last = n + 2
        rel_n = [pdb2r[x] for x in variant['release_n']]
        rel_c = [pdb2r[x] for x in variant['release_c']]
        # Chain order: [released motif C-term] -> Cys(last) -> S-S -> Cys(1) -> [released motif N-term]
        loop = sorted(rel_c) + [last, 1] + sorted(rel_n)
        close = {'res1': last, 'atom1': 'SG', 'res2': 1, 'atom2': 'SG',
                 'prior': 'CB', 'following': 'CB', 'geom': SSBOND}
        prepend, append = added[0], added[1]
    else:
        rel_n = [pdb2r[x] for x in variant['release_n']]
        rel_c = [pdb2r[x] for x in variant['release_c']]
        # Chain order: [released motif C-term] -> added -> new amide -> [released motif N-term]
        loop = sorted(rel_c) + list(range(n + 1, last + 1)) + sorted(rel_n)
        close = {'res1': last, 'atom1': 'C', 'res2': 1, 'atom2': 'N',
                 'prior': 'CA', 'following': 'CA', 'geom': PEPBOND}
        prepend, append = None, None
    if len(loop) < 3:
        raise ValueError('%s: GenKIC loop would be %d residues, below its 3-residue minimum'
                         % (variant['id'], len(loop)))
    pivots = pick_pivots(loop)
    return {'n': n, 'added': added, 'last': last, 'loop': loop, 'pivots': pivots,
            'release_n': rel_n, 'release_c': rel_c, 'close': close, 'mutate': [],
            'pdb2r': pdb2r, 'r2pdb': r2pdb, 'sample': True, 'prepend': prepend, 'append': append,
            'motif_range': (2, n + 1) if variant['route'] == 'ssapp' else (1, n),
            'polymer_termini': {1, last} if variant['route'] == 'ssapp' else set(),
            'declare': True,
            'ss_torsions': ss_torsions(last, 1) if variant['route'] == 'ssapp' else []}


def runs(indices):
    """[13,14,1,2,3] -> [(1,3),(13,14)]: contiguous spans, for MoveMap, which takes ranges only."""
    out = []
    for i in sorted(indices):
        if out and i == out[-1][1] + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return [tuple(x) for x in out]


def xml(variant, p, protocol):
    '''The RosettaScripts protocol for one variant, as a string.'''
    sf = protocol['scorefxn']
    out = ['<ROSETTASCRIPTS>', '  <SCOREFXNS>',
           '    <ScoreFunction name="sfxn" weights="%s"/>' % sf,
           '    <ScoreFunction name="sfxn_cst" weights="%s">' % sf]
    for t in ('atom_pair_constraint', 'angle_constraint', 'dihedral_constraint'):
        out.append('      <Reweight scoretype="%s" weight="1.0"/>' % t)
    out += ['    </ScoreFunction>', '  </SCOREFXNS>', '  <RESIDUE_SELECTORS>',
            '    <Index name="first" resnums="1"/>',
            '    <Index name="last" resnums="%d"/>' % p['last'],
            '    <Index name="motif" resnums="%d-%d"/>' % p['motif_range'],
            '  </RESIDUE_SELECTORS>', '  <MOVERS>']
    if p['added']:
        out.append('    <PeptideStubMover name="build">')
        if p['prepend']:
            # Prepend first: it renumbers the pose, and plan() already accounted for the shift.
            out.append('      <Prepend resname="%s" anchor_rsd="1"/>' % p['prepend'])
            out.append('      <Append resname="%s" anchor_rsd="%d"/>' % (p['append'], p['n'] + 1))
        else:
            for i, name in enumerate(p['added']):
                anchor = ' anchor_rsd="%d"' % p['n'] if i == 0 else ''
                out.append('      <Append resname="%s"%s/>' % (name, anchor))
        out.append('    </PeptideStubMover>')
    for num, name in p['mutate']:
        out.append('    <MutateResidue name="mut%d" target="%d" new_res="%s" preserve_atom_coords="0"/>'
                   % (num, num, name))
    # A head-to-tail macrocycle has no termini; the patches must come off before a bond can be
    # declared across them. A disulfide-closed ring still HAS termini -- the backbone really is
    # linear and the ring closes through S-S -- so stripping them there leaves residue 1 with nothing
    # bonded to its N, and Rosetta then scores its rama at phi = 0, a forbidden region, inventing
    # about +10 REU of phantom energy. That was emitted unconditionally until 2026-10-29 and biased
    # the FastRelax of every ssapp and ssmut run before then; see "Known caveats" in the README.
    strip_termini = p['close']['atom1'] != 'SG'
    if strip_termini:
        out += ['    <ModifyVariantType name="rm_upper" remove_type="UPPER_TERMINUS_VARIANT" residue_selector="last"/>',
                '    <ModifyVariantType name="rm_lower" remove_type="LOWER_TERMINUS_VARIANT" residue_selector="first"/>']
    c = p['close']
    if p['declare']:
        out.append('    <DeclareBond name="declare" res1="%d" atom1="%s" res2="%d" atom2="%s" '
                   'add_termini="0"/>' % (c['res1'], c['atom1'], c['res2'], c['atom2']))
    if p['sample']:
        out += ['    <GeneralizedKIC name="genkic" closure_attempts="%d" stop_when_n_solutions_found="%d"'
                % (protocol['closure_attempts'], protocol['min_solutions']),
                '                    selector="lowest_energy_selector" selector_scorefunction="sfxn">']
        for r in p['loop']:
            out.append('      <AddResidue res_index="%d"/>' % r)
        (p1, a1), (p2, a2), (p3, a3) = p['pivots']
        out.append('      <SetPivots res1="%d" atom1="%s" res2="%d" atom2="%s" res3="%d" atom3="%s"/>'
                   % (p1, a1, p2, a2, p3, a3))
        g = c['geom']
        out.append('      <CloseBond prioratom_res="%d" prioratom="%s" res1="%d" atom1="%s" '
                   'res2="%d" atom2="%s" followingatom_res="%d" followingatom="%s" '
                   'bondlength="%s" angle1="%s" angle2="%s" randomize_flanking_torsions="true"/>'
                   % (c['res1'], c['prior'], c['res1'], c['atom1'], c['res2'], c['atom2'],
                      c['res2'], c['following'], g['length'], g['angle1'], g['angle2']))
        # Each residue's own Ramachandran table does the seeding: D-Pro cannot propose negative phi.
        # randomize_backbone_by_rama_prepro aborts on a residue lacking either polymer connection, so
        # in the ssapp topology (ring closes through the S-S) the two Cys are excluded from it.
        perturb = [r for r in p['loop'] if r not in p['polymer_termini']]
        if perturb:
            out.append('      <AddPerturber effect="randomize_backbone_by_rama_prepro">')
            for r in perturb:
                out.append('        <AddResidue index="%d"/>' % r)
            out.append('      </AddPerturber>')
        # ...but leaving them unsampled entirely is why ssapp closed 0/50 on the first submission:
        # PeptideStubMover builds them pointing off into space and nothing ever moves them, so the two
        # sulfurs can never meet (100,000 attempts, zero solutions, no error). simple_cycpep_predict's
        # terminal_disulfide mode randomizes them with explicit randomize_dihedral perturbers instead --
        # the S-S torsion, each Cys side chain, and the Cys backbone - which is what this mirrors.
        for tors in p['ss_torsions']:
            out.append('      <AddPerturber effect="randomize_dihedral">')
            for (r1, a1), (r2, a2) in tors:
                out.append('        <AddAtoms res1="%d" atom1="%s" res2="%d" atom2="%s"/>'
                           % (r1, a1, r2, a2))
            out.append('      </AddPerturber>')
        # Pin every sampled peptide bond to trans.
        #
        # This is the fix for the defect that invalidated the first design set (2026-10-30). With
        # CA pivots the closure algorithm solves phi/psi at the pivots and leaves **omega free** --
        # the note above TARGET_LOOP in variants.py says so explicitly -- and nothing downstream
        # constrains it. The result: 47% of the bonds flanking the closure came out more than 30
        # deg from planar and 50% came out non-proline cis, against ~0% more than three bonds away.
        # PeptideCyclizeMover constrains the closure bond itself, which is why that one bond was
        # always clean, but it does not touch the neighbours.
        #
        # Because omega is NOT a solved DOF, setting it here holds through the closure solve.
        # Perturbers are emitted before the filters deliberately: GenKIC applies them in order.
        if protocol.get('omega_trans'):
            ss_pairs = {(min(a, b), max(a, b))
                        for tors in p['ss_torsions'] for (a, _), (b, _) in tors}
            loop = p['loop']
            for i in range(len(loop) - 1):
                r1, r2 = loop[i], loop[i + 1]
                if (min(r1, r2), max(r1, r2)) in ss_pairs:
                    continue          # disulfide, not a peptide bond
                out += ['      <AddPerturber effect="set_dihedral">',
                        '        <AddAtoms res1="%d" atom1="C" res2="%d" atom2="N"/>' % (r1, r2),
                        '        <AddValue value="180"/>',
                        '      </AddPerturber>']
        out.append('      <AddFilter type="loop_bump_check"/>')
        seed = seeded(variant, p['n'])
        for r in p['loop']:
            cut = protocol['seed_rama_cutoff'] if r == seed else protocol['rama_cutoff']
            if cut is not None:
                out.append('      <AddFilter type="rama_prepro_check" residue="%d" '
                           'rama_cutoff_energy="%s"/>' % (r, cut))
        out.append('    </GeneralizedKIC>')
    # Hold the motif rigid. An unrestrained FastRelax moves it ~5 A and throws away the whole point of
    # grafting a ddG-validated conformation, so relax is confined to the GenKIC loop by MoveMap and the
    # motif is additionally pinned by CA coordinate constraints.
    out.append('    <AddConstraintsToCurrentConformationMover name="pin" use_distance_cst="0" '
               'coord_dev="%s" residue_selector="motif" CA_only="1"/>' % protocol['coord_dev'])
    c2 = p['close']
    if c2['atom1'] == 'SG':
        out.append('    <ForceDisulfides name="cyclize" disulfides="%d:%d"/>'
                   % (min(c2['res1'], c2['res2']), max(c2['res1'], c2['res2'])))
    else:
        out.append('    <PeptideCyclizeMover name="cyclize"/>')
    out.append('    <FastRelax name="relax" scorefxn="sfxn_cst" repeats="%d">' % protocol['relax_rounds'])
    out.append('      <MoveMap name="mm" bb="0" chi="0" jump="0">')
    for lo, hi in runs(p['loop'] or []):
        out.append('        <Span begin="%d" end="%d" bb="1" chi="1"/>' % (lo, hi))
    if not p['loop']:   # ssmut: nothing was sampled, so let the two new Cys repack but keep bb fixed
        for num, _ in p['mutate']:
            out.append('        <Span begin="%d" end="%d" bb="0" chi="1"/>' % (num, num))
    out.append('      </MoveMap>')
    out.append('    </FastRelax>')
    out += ['  </MOVERS>', '  <PROTOCOLS>']
    if p['added']:
        out.append('    <Add mover="build"/>')
    for num, _ in p['mutate']:
        out.append('    <Add mover="mut%d"/>' % num)
    if strip_termini:
        out += ['    <Add mover="rm_upper"/>', '    <Add mover="rm_lower"/>']
    if p['declare']:
        out.append('    <Add mover="declare"/>')
    if p['sample']:
        out.append('    <Add mover="genkic"/>')
    out += ['    <Add mover="pin"/>', '    <Add mover="cyclize"/>', '    <Add mover="relax"/>', '  </PROTOCOLS>',
            '  <OUTPUT scorefxn="sfxn"/>', '</ROSETTASCRIPTS>']
    return '\n'.join(out) + '\n'


def write_variant(cfg, variant, residues, dest):
    '''Write inputs/<variant id>/{protocol.xml,input.json} and return the directory.'''
    os.makedirs(dest, exist_ok=True)
    p = plan(variant, residues)
    open(os.path.join(dest, 'protocol.xml'), 'w').write(xml(variant, p, cfg['protocol']))
    manifest = dict(variant)
    manifest.update({'pdb': os.path.relpath(cfg['scaffolds'][variant['scaffold']], dest),
                     'motif_first_pdb': residues[0].num, 'motif_last_pdb': residues[-1].num,
                     'motif_seq': struct.sequence(residues), 'n_motif': p['n'],
                     'rosetta_loop': p['loop'], 'rosetta_last': p['last'],
                     'motif_range': list(p['motif_range']),
                     'rosetta_to_pdb': {str(k): v for k, v in p['r2pdb'].items()}})
    json.dump(manifest, open(os.path.join(dest, 'input.json'), 'w'), indent=1, sort_keys=True)
    return dest
