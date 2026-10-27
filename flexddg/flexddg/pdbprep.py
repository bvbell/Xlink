'''PDB cleanup for Rosetta: first altloc only, HETATM -> ATOM, renames for residue types Rosetta lacks.'''
import numpy as np

# Residue types Rosetta has no params for, modeled as their closest standard residue.
# PCA (N-terminal pyroglutamate) -> GLN: OE -> OE1, NE2 is rebuilt by Rosetta.
RENAME = {'PCA': ('GLN', {' OE ': ' OE1'})}


def read_atoms(path):
    '''ATOM/HETATM lines, first altloc only (altloc column blanked), without newlines.'''
    atoms = []
    for l in open(path):
        if not l.startswith(('ATOM', 'HETATM')):
            continue
        if l[16] not in (' ', 'A'):
            continue
        atoms.append((l[:16] + ' ' + l[17:]).rstrip('\n'))
    return atoms


def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def res_key(l):
    '''(chain, resnum, icode, resname)'''
    return (l[21], int(l[22:26]), l[26], l[17:20])


def write_clean(atoms, path):
    with open(path, 'w') as f:
        for l in atoms:
            if l[17:20] in RENAME:
                new, names = RENAME[l[17:20]]
                l = l[:12] + names.get(l[12:16], l[12:16]) + l[16:17] + new + l[20:]
            f.write('ATOM  ' + l[6:] + '\n' if l.startswith('HETATM') else l + '\n')
        f.write('END\n')


def chain_breaks(atoms, cutoff=2.0):
    '''Consecutive residues of one chain whose C(i)-N(i+1) distance exceeds cutoff A.

    Rosetta reads these as bonded neighbors (no terminus variants are added), so the backbone is
    treated as continuous across the gap. Breaks far from the interface are harmless for ddG.
    '''
    C, N, order = {}, {}, []
    for l in atoms:
        k = res_key(l)[:3]
        if not order or order[-1] != k:
            order.append(k)
        if l[12:16] == ' C  ':
            C[k] = xyz(l)
        elif l[12:16] == ' N  ':
            N[k] = xyz(l)
    out = []
    for a, b in zip(order, order[1:]):
        if a[0] == b[0] and a in C and b in N:
            d = float(np.linalg.norm(C[a] - N[b]))
            if d > cutoff:
                out.append((a, b, d))
    return out


def sequence_by_residue(atoms):
    '''{(chain, resnum, icode): resname}'''
    return {res_key(l)[:3]: res_key(l)[3] for l in atoms}
