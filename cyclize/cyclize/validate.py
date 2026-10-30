'''Superimpose closed designs back onto the parent complex and filter them.

Two things can go wrong with a macrocycle that closed perfectly:

  * the motif drifted, so the ddG-validated conformation is gone;
  * the ring closed on the receptor side, so the added residues clash with the target.

Both are invisible without the receptor, and the scaffolds are excised motifs with no receptor chain,
so every design is superimposed onto the complex by its motif CA atoms first. Note that Rosetta does
not preserve the input coordinate frame, so RMSD must be measured AFTER superposition -- comparing raw
coordinates reports several angstroms of pure rigid-body displacement and looks like failure.

GEOMETRY GATES (added 2026-10-30). The two checks above are necessary and were not sufficient: the
entire first design set -- 311 scored complexes -- passed them while being chemically invalid. The
gates below close the holes that let that happen. Each exists because something specific got through:

  * `cis`          every peptide bond must be trans. GenKIC with CA pivots solves phi/psi at the
                   pivots and leaves omega FREE (see variants.TARGET_LOOP), so 50% of the bonds
                   flanking the closure came out cis. The user's rule is NO cis bonds at all, not
                   merely no non-proline cis -- cis-Pro is legitimate in nature but is a liability
                   in a synthesised macrocycle, and both input nanobodies are 100% trans, so
                   nothing native is lost by forbidding it.
  * `omega`        every peptide bond must be within `omega_tol` of planar. 47% of the flanking
                   bonds were >30 deg off, up to a fully broken 90.
  * `closure`      route-aware. Amide: the carbonyl O must sit ~2.25 A from the next N -- it was as
                   close as 0.14 A in 57% of amide designs, and NOTHING could see it because O and
                   N are a 1-3 pair through C, which fa_rep excludes. Disulfide: SG-SG ~2.05 A with
                   |CB-SG-SG-CB| in the gauche window -- ssapp ran a median 2.28 A and as high as
                   5.82 A, which is not a bond.
  * `intra`        the design must not pass through itself. The old receptor-clash test never
                   compared the design to itself. Polar pairs get a slacker cutoff because a
                   2.8-3.0 A N...O contact is a hydrogen bond, not a clash.
  * `chirality`    no unintended D residues. uOR's motif residues 102-107 are D in the input while
                   carrying standard L names, so Rosetta silently rebuilt them L and every design
                   was superimposed onto mirror-imaged coordinates. Expected D positions are passed
                   in, since D residues at the turn are deliberate.

The receptor-clash test also no longer skips the motif. It did, for a defensible reason -- without
the skip, native contacts like Y102 are flagged and every structure fails -- but the skip assumed
the motif keeps its native sidechain geometry, which is exactly what broke. `motif_clash` is the
separate, slacker threshold for motif residues.
'''
import math

from . import struct

# Peptide-bond geometry. omega is +-180 for trans; anything within `omega_tol` of 0 is cis.
OMEGA_TOL = 30.0            # degrees from planarity (0 or 180) before a bond is called distorted
CIS_CUTOFF = 90.0           # |omega| below this is cis
CLOSURE_ON = (2.05, 2.45)   # amide: carbonyl O to the next N, ideal 2.25
SS_DIST = (1.90, 2.20)      # disulfide: SG-SG, ideal 2.05
SS_CHI3 = (60.0, 120.0)     # disulfide: |CB-SG-SG-CB|, gauche
INTRA_MIN = 3.0             # non-adjacent intra-peptide heavy atoms, non-polar pair
INTRA_MIN_POLAR = 2.6       # N/O pairs may hydrogen bond


def _dihedral(p0, p1, p2, p3):
    import numpy as np
    b0, b1, b2 = (np.asarray(p0, float) - p1, np.asarray(p2, float) - p1,
                  np.asarray(p3, float) - p2)
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - b0.dot(b1) * b1
    w = b2 - b2.dot(b1) * b1
    return math.degrees(math.atan2(np.cross(b1, w).dot(v), v.dot(w)))


def _planarity(omega):
    '''-> degrees from the nearest planar value (0 for cis, 180 for trans).'''
    return min(abs(omega), abs(abs(omega) - 180.0))


def ring_bonds(design, cyclic=True):
    '''-> [(res_i, res_j, omega)] for each peptide bond, following the ring if `cyclic`.

    `design` is an ordered residue list. The last->first bond is included only for amide closures;
    a disulfide ring has no peptide bond there and would report a meaningless dihedral.
    '''
    out = []
    n = len(design)
    span = n if cyclic else n - 1
    for i in range(span):
        a, b = design[i], design[(i + 1) % n]
        try:
            out.append((a, b, _dihedral(a.atoms['CA'], a.atoms['C'],
                                        b.atoms['N'], b.atoms['CA'])))
        except KeyError:
            continue
    return out


def closure_geometry(design, route):
    '''-> dict of the route-appropriate closure measurements, with None where not applicable.'''
    out = {'closure_ON': None, 'ss_dist': None, 'ss_chi3': None}
    import numpy as np
    if route in ('ssapp', 'ssmut'):
        cys = [r for r in design if 'SG' in r.atoms]
        if len(cys) >= 2:
            a, b = cys[0], cys[-1]
            out['ss_dist'] = float(np.linalg.norm(
                np.asarray(a.atoms['SG']) - np.asarray(b.atoms['SG'])))
            if 'CB' in a.atoms and 'CB' in b.atoms:
                out['ss_chi3'] = _dihedral(a.atoms['CB'], a.atoms['SG'],
                                           b.atoms['SG'], b.atoms['CB'])
    else:
        last, first = design[-1], design[0]
        if 'O' in last.atoms and 'C' in last.atoms and 'N' in first.atoms:
            out['closure_ON'] = float(np.linalg.norm(
                np.asarray(last.atoms['O']) - np.asarray(first.atoms['N'])))
    return out


def min_intra_distance(design):
    '''-> (worst_gap, is_polar_pair) over non-adjacent intra-peptide heavy atoms.

    Ring neighbours are excluded at both seams. Polar (N/O to N/O) pairs are reported separately
    because they legitimately approach to ~2.8 A as hydrogen bonds.
    '''
    import numpy as np
    best, best_polar = float('inf'), float('inf')
    n = len(design)
    for i in range(n):
        for j in range(i + 1, n):
            if j - i == 1 or (i == 0 and j == n - 1):
                continue
            for an, x in design[i].atoms.items():
                for bn, y in design[j].atoms.items():
                    d = float(np.linalg.norm(np.asarray(x) - np.asarray(y)))
                    if an[0] in 'NO' and bn[0] in 'NO':
                        best_polar = min(best_polar, d)
                    else:
                        best = min(best, d)
    return best, best_polar


def d_residues(design):
    '''-> [res] whose backbone chirality is D (improper N-CA-C-CB < 0).

    Sign convention fixed against experimental 9v32, where all 545 residues read positive.
    '''
    out = []
    for r in design:
        if not all(a in r.atoms for a in ('N', 'CA', 'C', 'CB')):
            continue
        if _dihedral(r.atoms['N'], r.atoms['CA'], r.atoms['C'], r.atoms['CB']) < 0:
            out.append(r)
    return out


def kabsch(mobile, target):
    '''-> (rotation 3x3, mobile centroid, target centroid) minimising RMSD of paired point lists.'''
    import numpy as np
    P, Q = np.asarray(mobile, float), np.asarray(target, float)
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T.dot(Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T.dot(U.T)))
    R = Vt.T.dot(np.diag([1.0, 1.0, d])).dot(U.T)
    return R, pc, qc


def apply_transform(points, R, pc, qc):
    import numpy as np
    return (np.asarray(points, float) - pc).dot(R.T) + qc


def rmsd(a, b):
    return math.sqrt(sum(struct.dist(x, y) ** 2 for x, y in zip(a, b)) / len(a))


def superimpose_on_motif(design, native_motif, atom='CA', offset=0):
    '''Fit `design` onto the native motif using the motif residues they share.

    `design` and `native_motif` are residue lists. Pairing is positional, but `offset` is the number of
    design residues that precede the motif -- 1 for the ssapp route, which PREPENDS a Cys, so the motif
    starts at design position 2. Getting this wrong silently misaligns every pair and the name check
    then rejects all of them: the first apj run scored 0 of 150 ssapp structures for exactly that
    reason. Returns (transform, motif_rmsd, paired_count).
    '''
    pairs = [(d, n) for d, n in zip(design[offset:], native_motif)
             if atom in d.atoms and atom in n.atoms and d.name == n.name]
    if len(pairs) < 3:
        raise struct.StructError('only %d motif residues matched by name and position; the design and '
                                 'the scaffold are not the same motif' % len(pairs))
    R, pc, qc = kabsch([d.atoms[atom] for d, _ in pairs], [n.atoms[atom] for _, n in pairs])
    moved = apply_transform([d.atoms[atom] for d, _ in pairs], R, pc, qc)
    return (R, pc, qc), rmsd(moved, [n.atoms[atom] for _, n in pairs]), len(pairs)


def min_receptor_distance(design, receptor, transform, skip_nums=()):
    '''Closest heavy-atom approach between transformed design residues and the receptor.

    -> (distance, design residue, receptor residue). `skip_nums` drops the motif, so the number
    reported is about the residues the design ADDED, not the native contacts we want to keep.
    '''
    best = (float('inf'), None, None)
    rec_atoms = [(r, xyz) for r in receptor for xyz in r.atoms.values()]
    for d in design:
        if d.num in skip_nums:
            continue
        for xyz in apply_transform(list(d.atoms.values()), *transform):
            for r, rxyz in rec_atoms:
                v = struct.dist(tuple(xyz), rxyz)
                if v < best[0]:
                    best = (v, d, r)
    return best


def seed_phi(design, seed_num):
    '''Measured phi of the seeded turn residue, or None if it is not resolvable.'''
    pp = struct.phi_psi(design)
    return pp.get(seed_num, (None, None))[0]


def check(design_path, native_motif, receptor, motif_nums=None, clash=3.0, motif_rmsd=1.0,
          seed_num=None, offset=0, n_motif=None, route='amide', expect_d=(),
          geometry=True, omega_tol=OMEGA_TOL, motif_clash=2.2):
    '''-> dict verdict for one closed design.

    The motif is identified by POSITION in the design, not by residue number: routes that add residues
    make Rosetta renumber the pose 1..N, but ssmut adds nothing so the original PDB numbering survives.
    A numeric range therefore matches one and not the other -- which made uOR ssmut report its own
    motif hotspot Y102 as an "added residue" clashing with the receptor, failing all 50 structures.

    `route` selects the closure check. `expect_d` is the set of design residue numbers that are
    legitimately D (the turn seed may be a D residue) -- any OTHER D residue is a chirality error.
    `geometry=False` restores the pre-2026-10-30 behaviour, for reproducing archived results only.
    '''
    des = struct.chain_residues(struct.read(design_path), 'A')
    transform, rms, npair = superimpose_on_motif(des, native_motif, offset=offset)
    n_motif = n_motif if n_motif is not None else len(native_motif)
    motif_set = {r.num for r in des[offset:offset + n_motif]}
    dist, dres, rres = min_receptor_distance(des, receptor, transform, skip_nums=motif_set)
    # The motif is no longer exempt from the receptor check, only held to a slacker threshold: the
    # native interface genuinely packs motif residues to 2.44-2.80 A, but a transplanted sidechain
    # at 1.08 A is an error, and the old blanket skip could not tell those apart.
    m_dist, m_dres, m_rres = min_receptor_distance(
        des, receptor, transform, skip_nums={r.num for r in des} - motif_set)
    reasons = []
    if rms > motif_rmsd:
        reasons.append('motif drifted %.2f A (limit %.2f)' % (rms, motif_rmsd))
    if dist < clash:
        reasons.append('added residue %s within %.2f A of receptor %s (limit %.2f)'
                       % (dres, dist, rres, clash))
    if m_dist < motif_clash:
        reasons.append('motif residue %s within %.2f A of receptor %s (limit %.2f)'
                       % (m_dres, m_dist, m_rres, motif_clash))
    # The seed only counts if the closure algebra did not flip it: a D residue at i+1 must end up at
    # positive phi, or the turn is not the type that was seeded and the register will not continue.
    phi = seed_phi(des, seed_num) if seed_num else None
    if seed_num and phi is not None and phi <= 0:
        reasons.append('seeded residue %d came out at phi %.0f, wrong handedness for a I\'/II\' turn'
                       % (seed_num, phi))

    geo = {'closure_ON': None, 'ss_dist': None, 'ss_chi3': None,
           'omega_max_dev': None, 'n_cis': None, 'intra_min': None,
           'intra_min_polar': None, 'n_bad_chirality': None}
    if geometry:
        is_ss = route in ('ssapp', 'ssmut')
        bonds = ring_bonds(des, cyclic=not is_ss)
        if bonds:
            devs = [_planarity(o) for _, _, o in bonds]
            cis = [(a, b, o) for a, b, o in bonds if abs(o) < CIS_CUTOFF]
            geo['omega_max_dev'] = max(devs)
            geo['n_cis'] = len(cis)
            if geo['omega_max_dev'] > omega_tol:
                worst = max(bonds, key=lambda t: _planarity(t[2]))
                reasons.append('peptide bond %s-%s is %.0f deg from planar (limit %.0f)'
                               % (worst[0].num, worst[1].num, _planarity(worst[2]), omega_tol))
            # No cis bonds at all -- not merely no non-proline cis. Both input nanobodies are
            # 100% trans, so this forbids nothing that was native.
            if cis:
                reasons.append('%d cis peptide bond(s): %s' % (
                    len(cis), ', '.join('%s-%s %.0f deg' % (a.num, b.num, o) for a, b, o in cis)))
        geo.update(closure_geometry(des, route))
        if is_ss:
            if geo['ss_dist'] is None:
                reasons.append('no disulfide found for route %s' % route)
            else:
                if not SS_DIST[0] <= geo['ss_dist'] <= SS_DIST[1]:
                    reasons.append('S-S %.2f A outside %.2f-%.2f'
                                   % (geo['ss_dist'], SS_DIST[0], SS_DIST[1]))
                if geo['ss_chi3'] is not None and not (
                        SS_CHI3[0] <= abs(geo['ss_chi3']) <= SS_CHI3[1]):
                    reasons.append('disulfide chi3 %.0f deg outside the %.0f-%.0f gauche window'
                                   % (abs(geo['ss_chi3']), SS_CHI3[0], SS_CHI3[1]))
        elif geo['closure_ON'] is None:
            reasons.append('closure carbonyl O or next N missing')
        elif not CLOSURE_ON[0] <= geo['closure_ON'] <= CLOSURE_ON[1]:
            reasons.append('closure O...N %.2f A outside %.2f-%.2f (ideal 2.25)'
                           % (geo['closure_ON'], CLOSURE_ON[0], CLOSURE_ON[1]))

        intra, intra_polar = min_intra_distance(des)
        geo['intra_min'], geo['intra_min_polar'] = intra, intra_polar
        if intra < INTRA_MIN:
            reasons.append('design contacts itself at %.2f A (limit %.2f)' % (intra, INTRA_MIN))
        if intra_polar < INTRA_MIN_POLAR:
            reasons.append('polar intra-peptide contact %.2f A (limit %.2f)'
                           % (intra_polar, INTRA_MIN_POLAR))

        bad_d = [r for r in d_residues(des) if r.num not in set(expect_d)]
        geo['n_bad_chirality'] = len(bad_d)
        if bad_d:
            reasons.append('unexpected D residue(s): %s'
                           % ', '.join('%s%d' % (r.name, r.num) for r in bad_d))

    out = {'path': design_path, 'motif_rmsd': rms, 'n_paired': npair, 'min_receptor_dist': dist,
           'min_motif_receptor_dist': m_dist,
           'closest_design_res': str(dres) if dres else None,
           'closest_receptor_res': str(rres) if rres else None, 'seed_phi': phi,
           'pass': not reasons, 'reasons': '; '.join(reasons)}
    out.update(geo)
    return out
