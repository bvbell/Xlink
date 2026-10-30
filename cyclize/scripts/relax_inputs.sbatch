#!/bin/bash
#SBATCH --job-name=relax_inputs
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1
#SBATCH --output=/scratch/users/bvbell/xlink/cyclize/logs/relax_inputs_%A_%a.out
#
# Relax both input complexes so the motifs can be re-excised from a self-consistent pose.
#
# Why (2026-10-30): the whole design set was invalidated by geometry defects, and two of them
# originate upstream of GenKIC:
#   * uOR nanobody residues 102-107 are D-amino acids in uor_complex.pdb, labelled with standard
#     L three-letter codes. Rosetta assigns L types on read and rebuilds CB from ideal L geometry,
#     so every design is built L and then superimposed back onto D coordinates -- CB lands 2.0-2.7 A
#     off and the sidechain points into the receptor. Relaxing the input and re-excising the motif
#     from the SAME relaxed pose makes the motif source and the alignment target consistent, which
#     removes the mismatch by construction.
#   * apj's input is clean on chirality but still benefits from being the same object the designs
#     are later aligned to.
#
# NOTE the uOR caveat: this does not recover the experimental conformation. Rosetta's L rebuild
# moves uOR motif sidechain tips up to 6.7 A (Arg107) and 6.65 A (Phe104). The result is
# self-consistent but its relationship to the J126 density is unverified. apj has no such caveat.
#
# nstruct 10 per target; pick the lowest-scoring decoy that holds the motif in place.
set -uo pipefail

ml purge
ml biology rosetta/3.14

CYC=/scratch/users/bvbell/xlink/cyclize
ROOT=/scratch/users/bvbell/xlink
WORK=$CYC/scripts/relax_inputs
mkdir -p "$WORK" "$CYC/logs"

TARGETS=("apj $ROOT/apj_complex.pdb" "uor $ROOT/uor_complex.pdb")
read -r NAME PDB <<< "${TARGETS[$SLURM_ARRAY_TASK_ID]}"

OUT=$WORK/$NAME
mkdir -p "$OUT"

# CA coordinate restraints at sd 0.5 A keep the fold and the interface where they are; the point is
# to regularise local geometry (and, for uOR, to let Rosetta's L rebuild settle), not to redock.
cat > "$OUT/relax.xml" <<'EOF'
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn" weights="ref2015"/>
    <ScoreFunction name="sfxn_cst" weights="ref2015">
      <Reweight scoretype="coordinate_constraint" weight="1.0"/>
    </ScoreFunction>
  </SCOREFXNS>
  <MOVERS>
    <VirtualRoot name="vroot" removable="true"/>
    <AddConstraints name="addcst">
      <CoordinateConstraintGenerator name="coordgen" sd="0.5" ca_only="true"/>
    </AddConstraints>
    <FastRelax name="relax" scorefxn="sfxn_cst" repeats="3"/>
    <RemoveConstraints name="rmcst" constraint_generators="coordgen"/>
  </MOVERS>
  <PROTOCOLS>
    <Add mover="vroot"/>
    <Add mover="addcst"/>
    <Add mover="relax"/>
    <Add mover="rmcst"/>
  </PROTOCOLS>
  <OUTPUT scorefxn="sfxn"/>
</ROSETTASCRIPTS>
EOF

echo "=== $NAME relax starting on $(hostname) $(date) ==="
cd "$OUT"
rosetta_scripts.default.linuxgccrelease \
  -parser:protocol relax.xml \
  -s "$PDB" \
  -nstruct 10 \
  -out:path:pdb "$OUT" \
  -out:file:scorefile "$OUT/relax.sc" \
  -overwrite \
  -ignore_unrecognized_res \
  -ex1 -ex2
echo "=== $NAME exit $? $(date) ==="
