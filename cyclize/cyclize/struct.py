'''Minimal PDB reading and the geometry helpers the rest of the package needs.

Residues are keyed by (chain, resnum); insertion codes are rejected rather than silently dropped,
because every input here is an excised motif whose numbering must stay comparable to the parent
complex. Only heavy atoms are read.
'''
import math

AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
          'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
          'CYD': 'C'}
# Rosetta's chiral-flip patch renames D residues to these name3 codes, and writes them as HETATM
# rather than ATOM -- a reader that only takes ATOM records silently drops every designed turn residue.
D_AA3to1 = {'DAL': 'a', 'DAR': 'r', 'DAN': 'n', 'DAS': 'd', 'DCS': 'c', 'DGN': 'q', 'DGU': 'e',
            'DHI': 'h', 'DIL': 'i', 'DLE': 'l', 'DLY': 'k', 'DME': 'm', 'DPH': 'f', 'DPR': 'p',
            'DSE': 's', 'DTH': 't', 'DTR': 'w', 'DTY': 'y', 'DVA': 'v'}
AA3to1.update(D_AA3to1)
BACKBONE = ('N', 'CA', 'C', 'O')


class StructError(ValueError):
    pass


class Residue(object):
    __slots__ = ('chain', 'num', 'name', 'atoms')

    def __init__(self, chain, num, name):
        self.chain, self.num, self.name = chain, num, name
        self.atoms = {}

    @property
    def aa1(self):
        return AA3to1.get(self.name, 'X')

    def __repr__(self):
        return '%s%d' % (self.aa1, self.num)


def read(path, chain=None, hetatm=True):
    '''-> {(chain, num): Residue}, optionally a single chain.

    HETATM is read by default because Rosetta writes D-amino acids and other non-canonical residues
    that way; pass hetatm=False to read only ATOM.
    '''
    out = {}
    for line in open(path):
        if not (line.startswith('ATOM') or (hetatm and line.startswith('HETATM'))):
            continue
        ch = line[21]
        if chain is not None and ch != chain:
            continue
        if line[26] != ' ':
            raise StructError('%s: insertion code %r at %s%s; renumber first'
                              % (path, line[26], ch, line[22:26].strip()))
        num = int(line[22:26])
        res = out.get((ch, num))
        if res is None:
            res = out[(ch, num)] = Residue(ch, num, line[17:20].strip())
        name = line[12:16].strip()
        if line[76:78].strip() != 'H' and not name.startswith('H'):
            res.atoms[name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    if not out:
        raise StructError('%s: no residues read%s' % (path, '' if chain is None else ' for chain ' + chain))
    return out


def chain_residues(res, chain):
    '''Residues of one chain, in numbering order.'''
    return [res[k] for k in sorted(res) if k[0] == chain]


def sequence(residues):
    return ''.join(r.aa1 for r in residues)


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def dihedral(p1, p2, p3, p4):
    '''Signed dihedral in degrees.'''
    b0 = [p1[i] - p2[i] for i in range(3)]
    b1 = [p3[i] - p2[i] for i in range(3)]
    b2 = [p4[i] - p3[i] for i in range(3)]
    n = math.sqrt(sum(x * x for x in b1))
    if n == 0:
        raise StructError('coincident atoms in dihedral')
    b1 = [x / n for x in b1]
    d0 = sum(x * y for x, y in zip(b0, b1))
    d2 = sum(x * y for x, y in zip(b2, b1))
    v = [b0[i] - d0 * b1[i] for i in range(3)]
    w = [b2[i] - d2 * b1[i] for i in range(3)]
    cross = [b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0]]
    return math.degrees(math.atan2(sum(x * y for x, y in zip(cross, w)),
                                   sum(x * y for x, y in zip(v, w))))


def phi_psi(residues):
    '''-> {resnum: (phi, psi)} for residues with both neighbors present and bonded.'''
    by = {r.num: r for r in residues}
    out = {}
    for r in residues:
        prev, nxt = by.get(r.num - 1), by.get(r.num + 1)
        try:
            phi = dihedral(prev.atoms['C'], r.atoms['N'], r.atoms['CA'], r.atoms['C']) if prev else None
            psi = dihedral(r.atoms['N'], r.atoms['CA'], r.atoms['C'], nxt.atoms['N']) if nxt else None
        except (KeyError, AttributeError):
            continue
        if phi is not None or psi is not None:
            out[r.num] = (phi, psi)
    return out
