load /scratch/users/bvbell/xlink/9v32_complex.pdb, cplx
create tgtA, cplx and chain A and polymer
load /scratch/users/bvbell/xlink/cldn6.pdb, mob
super mob, tgtA
load /scratch/users/bvbell/xlink/cldn6_9v32fab_complex.pdb, chim
python
from pymol import cmd
d = cmd.rms_cur("mob", "chim and chain A", matchmaker=-1)
print("VERIFY rms_cur(pymol-super-cldn6, my-chainA) = %.6f A" % d)
# Fab coords unchanged?
dh = cmd.rms_cur("cplx and chain H", "chim and chain H", matchmaker=-1)
dl = cmd.rms_cur("cplx and chain L", "chim and chain L", matchmaker=-1)
print("VERIFY chain H rmsd to original = %.6f ; chain L = %.6f" % (dh, dl))
import numpy as np
M = cmd.get_object_matrix("mob")
R = np.array(M).reshape(4,4)[:3,:3]
print("VERIFY det(R) = %.8f" % np.linalg.det(R))
python end
