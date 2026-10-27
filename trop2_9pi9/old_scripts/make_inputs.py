"""One input directory per interface alanine mutation (moving chain E away from the AB Fab)."""
import os
three = {'ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
         'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
pdb = os.path.abspath('9pi9_clean.pdb')
names = []
for line in open('interface_residues.txt'):
    chain, num, resn = line.split()
    name = '%s_%s%sA' % (chain, three[resn], num)       # e.g. A_Y32A
    d = os.path.join('inputs', name)
    os.makedirs(d, exist_ok=True)
    if not os.path.exists(os.path.join(d, '9pi9_clean.pdb')):
        os.symlink(pdb, os.path.join(d, '9pi9_clean.pdb'))
    open(os.path.join(d, 'chains_to_move.txt'), 'w').write('E\n')
    open(os.path.join(d, 'nataa_mutations.resfile'), 'w').write('NATAA\nstart\n%s %s PIKAA A\n' % (num, chain))
    names.append(name)
open('mutations.txt', 'w').write('\n'.join(names) + '\n')
print(len(names), 'mutations')
