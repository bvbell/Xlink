'''Superimpose closed designs back onto the parent complex and filter them.

Two things can go wrong with a macrocycle that closed perfectly:

  * the motif drifted, so the ddG-validated conformation is gone;
  * the ring closed on the receptor side, so the added residues clash with the target.

Both are invisible without the receptor, and the scaffolds are excised motifs with no receptor chain,
so every design is superimposed onto the complex by its motif CA atoms first. Note that Rosetta does
not preserve the input coordinate frame, so RMSD must be measured AFTER superposition -- comparing raw
coordinates reports several angstroms of pure rigid-body displacement and looks like failure.
'''
import math

from . import struct


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


def check(design_path, native_motif, receptor, motif_nums, clash=3.0, motif_rmsd=1.0, seed_num=None,
          offset=0):
    '''-> dict verdict for one closed design.'''
    des = struct.chain_residues(struct.read(design_path), 'A')
    transform, rms, npair = superimpose_on_motif(des, native_motif, offset=offset)
    dist, dres, rres = min_receptor_distance(des, receptor, transform, skip_nums=set(motif_nums))
    reasons = []
    if rms > motif_rmsd:
        reasons.append('motif drifted %.2f A (limit %.2f)' % (rms, motif_rmsd))
    if dist < clash:
        reasons.append('added residue %s within %.2f A of receptor %s (limit %.2f)'
                       % (dres, dist, rres, clash))
    # The seed only counts if the closure algebra did not flip it: a D residue at i+1 must end up at
    # positive phi, or the turn is not the type that was seeded and the register will not continue.
    phi = seed_phi(des, seed_num) if seed_num else None
    if seed_num and phi is not None and phi <= 0:
        reasons.append('seeded residue %d came out at phi %.0f, wrong handedness for a I\'/II\' turn'
                       % (seed_num, phi))
    return {'path': design_path, 'motif_rmsd': rms, 'n_paired': npair, 'min_receptor_dist': dist,
            'closest_design_res': str(dres) if dres else None,
            'closest_receptor_res': str(rres) if rres else None, 'seed_phi': phi,
            'pass': not reasons, 'reasons': '; '.join(reasons)}
