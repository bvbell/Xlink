'''Beta-hairpin register analysis -- the thing every design decision here keys off.

An antiparallel hairpin pairs residue i with residue (K - i) for one constant K, so every backbone
N-H...O=C H-bond in the ladder has i + j == K. Both CDR3s analysed for this project give K = 209 with
no exceptions. Two consequences drive the whole protocol:

  * A motif is "on register" only if first + last == K. Because i and K - i are always an odd number
    of positions apart, an on-register motif always has an EVEN residue count -- parity is the symptom,
    the register is the cause.
  * The two termini of an on-register motif are a reciprocal H-bonded pair, which is why they sit
    4-6 A apart. Closing the ring must therefore add 0 or 2 residues; adding 1 shifts the register by
    one and unzips the ladder.

Disulfide geometry: natural disulfides sit at Cbeta-Cbeta ~3.4-4.6 A (mean ~3.85). Cbeta position is
fixed by the backbone, so a terminal pair outside that window cannot be stapled by mutation to Cys --
no amount of repacking will reach. Those scaffolds need Cys appended on free backbone instead.
'''
from . import struct

HBOND_CUTOFF = 3.4          # N...O heavy-atom distance, generous for a ladder survey
SEQ_SEP = 3                 # ignore |i-j| < 3: those are helical/local, not cross-strand
SS_CB_MIN, SS_CB_MAX = 3.4, 4.6


def hbonds(residues, cutoff=HBOND_CUTOFF):
    '''Cross-strand backbone H-bonds as (donor_num, acceptor_num, distance), donor N -> acceptor O.'''
    out = []
    for d in residues:
        if 'N' not in d.atoms:
            continue
        for a in residues:
            if abs(d.num - a.num) < SEQ_SEP or 'O' not in a.atoms:
                continue
            r = struct.dist(d.atoms['N'], a.atoms['O'])
            if r < cutoff:
                out.append((d.num, a.num, r))
    return sorted(out)


def register(residues, cutoff=HBOND_CUTOFF):
    '''-> (K, n_supporting, n_total): the modal i+j over cross-strand H-bonds.

    n_supporting < n_total means the segment is not a single clean hairpin (a bulge, a kink, or two
    stacked ladders); the caller should say so rather than trusting K.
    '''
    bonds = hbonds(residues, cutoff)
    if not bonds:
        return None, 0, 0
    sums = {}
    for i, j, _ in bonds:
        sums[i + j] = sums.get(i + j, 0) + 1
    k = max(sums, key=lambda s: (sums[s], -abs(s)))
    return k, sums[k], len(bonds)


def pairs_in_register(k, lo, hi):
    '''The H-bond-paired positions of register K that lie within [lo, hi].'''
    return [(i, k - i) for i in range(lo, (k // 2) + 1) if lo <= k - i <= hi and i < k - i]


def terminal_pair(residues, k=None):
    '''Geometry of the two ends that closure has to join.

    Returns a dict with the two termini, their sum (compare to K), Ca-Ca and Cb-Cb distances, and
    whether a disulfide can be made by mutating them (Cb-Cb in the natural window).
    '''
    first, last = residues[0], residues[-1]
    ca = cb = None
    if 'CA' in first.atoms and 'CA' in last.atoms:
        ca = struct.dist(first.atoms['CA'], last.atoms['CA'])
    if 'CB' in first.atoms and 'CB' in last.atoms:
        cb = struct.dist(first.atoms['CB'], last.atoms['CB'])
    return {'first': first, 'last': last, 'sum': first.num + last.num,
            'on_register': None if k is None else (first.num + last.num == k),
            'ca': ca, 'cb': cb,
            'staple_ok': cb is not None and SS_CB_MIN <= cb <= SS_CB_MAX}


# Turn classification. A 2-residue hairpin turn is i, i+1, i+2, i+3 with i and i+3 the H-bonded
# strand pair; only the sign of phi at i+1 and i+2 distinguishes the four types. Beta-sheets have a
# right-handed twist, so hairpins need the mirror types I' and II'; I and II impose the wrong twist.
TURN_TYPES = {(1, 1): "I'", (1, -1): "II'", (-1, -1): 'I', (-1, 1): 'II'}
HAIRPIN_TURNS = ("I'", "II'")


def turn_type(phi1, phi2):
    '''-> "I'", "II'", "I", "II", or None if either phi is missing.'''
    if phi1 is None or phi2 is None:
        return None
    return TURN_TYPES[(1 if phi1 > 0 else -1, 1 if phi2 > 0 else -1)]


def describe(residues, k=None, scan=None):
    '''Human-readable register report for one segment; `scan` maps resnum -> ddG for annotation.'''
    k_obs, n_ok, n_tot = register(residues, )
    k = k or k_obs
    lines = ['%s%d-%s%d  n=%d  %s' % (residues[0].aa1, residues[0].num, residues[-1].aa1,
                                      residues[-1].num, len(residues), struct.sequence(residues))]
    if k_obs is None:
        lines.append('  no cross-strand H-bonds: not a hairpin')
        return '\n'.join(lines)
    lines.append('  register K=%d (%d/%d H-bonds agree)%s'
                 % (k_obs, n_ok, n_tot, '' if n_ok == n_tot else '  <-- not a single clean ladder'))
    for i, j, r in hbonds(residues):
        if i < j:
            lines.append('    %d N...O %d  %.2f A%s' % (i, j, r, '' if i + j == k else '  off-register'))
    t = terminal_pair(residues, k)
    lines.append('  termini %s%d+%s%d = %d (%s)'
                 % (t['first'].aa1, t['first'].num, t['last'].aa1, t['last'].num, t['sum'],
                    'on register' if t['on_register'] else 'OFF REGISTER'))
    lines.append('    Ca-Ca %s  Cb-Cb %s  -> mutate-to-Cys staple %s'
                 % ('%.1f' % t['ca'] if t['ca'] else '?', '%.1f' % t['cb'] if t['cb'] else '?',
                    'viable' if t['staple_ok'] else 'too far (append Cys instead)'))
    if scan:
        lines.append('    terminal ddG: %s=%s  %s=%s'
                     % (t['first'], _fmt(scan.get(t['first'].num)), t['last'], _fmt(scan.get(t['last'].num))))
    return '\n'.join(lines)


def _fmt(v):
    return 'n/a' if v is None else '%.2f' % v
