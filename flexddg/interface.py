'''Interface residues between the moving group and everything else.'''
import numpy as np

from .pdbprep import res_key, xyz


def interface_residues(atoms, moving, partner, cutoff):
    '''{(chain, resnum, icode, resname): min heavy-atom distance} for residues within cutoff A.'''
    heavy = [l for l in atoms if l[76:78].strip() != 'H']
    a = [l for l in heavy if l[21] in partner]
    b = [l for l in heavy if l[21] in moving]
    if not a or not b:
        raise ValueError('no atoms for chains %s or %s' % (''.join(partner), ''.join(moving)))
    XA = np.array([xyz(l) for l in a]); XB = np.array([xyz(l) for l in b])
    d = np.sqrt(((XA[:, None, :] - XB[None, :, :]) ** 2).sum(-1))
    iface = {}
    for i, j in zip(*np.where(d < cutoff)):
        for l in (a[i], b[j]):
            k = res_key(l)
            iface[k] = min(iface.get(k, 99), d[i, j])
    return iface


def disulfide_cys(atoms):
    '''Residue keys of cysteines whose SG is within 2.5 A of another SG.'''
    sg = [l for l in atoms if l[17:20] == 'CYS' and l[12:16] == ' SG ' and l[76:78].strip() != 'H']
    if len(sg) < 2:
        return set()
    S = np.array([xyz(l) for l in sg])
    ds = np.sqrt(((S[:, None] - S[None]) ** 2).sum(-1))
    return {res_key(sg[i]) for i, j in zip(*np.where((ds < 2.5) & (ds > 0)))}


def select(iface, skip, disulf, chain_order):
    '''Sort interface residues by chain order then number; split into (keep, skipped).'''
    keep, skipped = [], []
    for k in sorted(iface, key=lambda k: (chain_order.index(k[0]), k[1], k[2])):
        if k[3] in skip or ('disulfide' in skip and k in disulf):
            skipped.append(k)
        else:
            keep.append(k)
    return keep, skipped
