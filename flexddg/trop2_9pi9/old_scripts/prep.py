"""Clean 9pi9_complex.pdb and pick AB<->E interface residues for alanine scanning."""
import numpy as np
src = '../9pi9_complex.pdb'
atoms = []
for l in open(src):
    if not l.startswith(('ATOM', 'HETATM')):
        continue
    alt = l[16]
    if alt not in (' ', 'A'):
        continue
    l = l[:16] + ' ' + l[17:]
    atoms.append(l.rstrip('\n'))

xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
res = lambda l: (l[21], int(l[22:26]), l[17:20])
heavy = [l for l in atoms if l[76:78].strip() != 'H']

ab = [l for l in heavy if l[21] in 'AB']
e = [l for l in heavy if l[21] == 'E']
XA = np.array([xyz(l) for l in ab]); XE = np.array([xyz(l) for l in e])
d = np.sqrt(((XA[:, None, :] - XE[None, :, :]) ** 2).sum(-1))

CUT = 5.0
iface = {}
for i, j in zip(*np.where(d < CUT)):
    for l, dist in ((ab[i], d[i, j]), (e[j], d[i, j])):
        k = res(l)
        iface[k] = min(iface.get(k, 99), dist)

# PCA E31 distance to partner
pca = [i for i, l in enumerate(e) if l[17:20] == 'PCA']
print('PCA E31 min dist to AB: %.1f' % d[:, pca].min())

# disulfide cysteines (SG-SG < 2.5 A)
sg = [l for l in heavy if l[17:20] == 'CYS' and l[12:16] == ' SG ']
S = np.array([xyz(l) for l in sg])
ds = np.sqrt(((S[:, None] - S[None]) ** 2).sum(-1))
disulf = {res(sg[i]) for i, j in zip(*np.where((ds < 2.5) & (ds > 0)))}

keep, skipped = [], []
for k in sorted(iface, key=lambda k: ('ABE'.index(k[0]), k[1])):
    if k[2] in ('ALA', 'GLY', 'PCA') or k in disulf:
        skipped.append(k)
    else:
        keep.append(k)
print('interface residues (<%.1f A heavy atom):' % CUT, len(iface))
print('skipped (Ala/Gly/PCA/disulfide Cys):', ' '.join('%s%d%s' % (c, n, r) for c, n, r in skipped))
print('to scan:', len(keep))
for c, n, r in keep:
    print('  %s %4d %s  %.1f A' % (c, n, r, iface[(c, n, r)]))

with open('9pi9_clean.pdb', 'w') as f:
    for l in atoms:
        if l[17:20] == 'PCA':
            # N-terminal pyroglutamate contacts the Fab; Rosetta has no PCA type, so model it as GLN
            # (OE -> OE1, NE2 is rebuilt by Rosetta)
            name = ' OE1' if l[12:16] == ' OE ' else l[12:16]
            l = l[:12] + name + l[16:17] + 'GLN' + l[20:]
        f.write('ATOM  ' + l[6:] + '\n' if l.startswith('HETATM') else l + '\n')
    f.write('END\n')
with open('interface_residues.txt', 'w') as f:
    for c, n, r in keep:
        f.write('%s %d %s\n' % (c, n, r))
