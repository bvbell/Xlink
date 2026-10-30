#!/bin/bash
# Constrained FastRelax of the CLDN6/9v32-Fab superposition chimera.
#
# Why: the raw chimera has 56 heavy-atom pairs under 2.0 A against the Fab (worst 0.45 A),
# several of them backbone-backbone at the CLDN6 beta1-beta2 loop (39-44) vs CDR-L1/L3.
# flex ddG's own opening MinMover is CA-restrained at 0.5 A and will not clear those.
#
# Two constraint strengths are run side by side so we can pick on evidence rather than taste:
#   tasks 0-9   coord_cst_stdev 0.5  (tight: stays close to the superposition)
#   tasks 10-19 coord_cst_stdev 1.0  (loose: lets the clashing loop actually move)
#
# -in:missing_density_to_jump: chain H has two pre-existing gaps (131->139, 188->193) that
# Rosetta would otherwise read as 13.8 A and 12.5 A covalent bonds. During a bounded ddG
# minimization that is tolerable (and is how the other targets were run), but FastRelax is
# aggressive enough that the phantom bond strain would distort the constant domain. The jump
# applies to this prep step only; the ddG run reads the relaxed PDB with the usual flags.
#
# talaris2014 (-restore_talaris_behavior) to match the scorefunction the ddG protocol uses.
set -uo pipefail

ml purge
ml biology rosetta/3.14

RELAXDIR=/scratch/users/bvbell/xlink/flexddg/cldn6_9v32fab/relax
CHIMERA=/scratch/users/bvbell/xlink/cldn6_9v32fab_complex.pdb

i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 10 ]; then VAR=sd05; SD=0.5; else VAR=sd10; SD=1.0; fi
OUT=$RELAXDIR/decoys/$VAR/$i
mkdir -p "$OUT"

echo "relax task $i variant $VAR coord_cst_stdev $SD on $(hostname)"

exec relax.default.linuxgccrelease \
  -s "$CHIMERA" \
  -nstruct 1 \
  -relax:constrain_relax_to_start_coords \
  -relax:ramp_constraints false \
  -relax:coord_cst_stdev "$SD" \
  -in:missing_density_to_jump \
  -restore_talaris_behavior \
  -in:file:fullatom \
  -ignore_unrecognized_res \
  -ignore_zero_occupancy false \
  -ex1 -ex2 -use_input_sc \
  -constant_seed -jran $((1000 + i)) \
  -out:path:all "$OUT" \
  -out:file:scorefile score.sc \
  -overwrite
