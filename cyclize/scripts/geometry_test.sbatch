#!/bin/bash
#SBATCH --job-name=geom_test
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-23
#SBATCH --output=/scratch/users/bvbell/xlink/cyclize/logs/geom_test_%A_%a.out
#
# Small test: does relaxing the input, and/or pinning omega trans, produce chemically valid
# macrocycles? Run this BEFORE any production rebuild.
#
# Background (2026-10-30): all 311 designs in the first set were invalid. Three defects:
#   1. closure carbonyl O misplaced in 57% of amides (O...N down to 0.14 A, ideal 2.25)
#   2. 47% of bonds flanking the closure >30 deg from planar, 50% non-proline cis
#   3. uOR motif residues 102-107 are D-amino acids in the input, rebuilt L by Rosetta on read
# Cause of 1 and 2: GenKIC with CA pivots solves phi/psi at the pivots and leaves omega FREE
# (variants.TARGET_LOOP says so), and PeptideCyclizeMover constrains only the closure bond.
#
# Conditions, 2 targets x 2 scaffolds x 3 arms x 2 route groups = 24 tasks:
#   base   relaxed input, omega unconstrained      -- control; should reproduce the failure
#   omega  relaxed input, omega_trans=True         -- the fix under test
#   loop7  omega_trans=True, target_loop=7         -- extends the sampled loop by releasing two
#                                                     more COLD terminal motif residues. Still
#                                                     adds only 2 residues to the motif; released
#                                                     residues keep their identity (no design).
#
# Both closure chemistries are tested. The disulfide routes are independently broken in the
# archived set -- ssapp S-S distance median 2.28 A against an ideal 2.05 and as high as 5.82 A
# (not a bond at all), with |chi3| outside the 60-120 deg gauche window in 43 of 74; the one
# ssmut structure is 1.70 A at chi3 164 deg. omega_trans pins peptide bonds and deliberately
# skips the S-S, so for these routes it fixes the FLANKING bonds while ForceDisulfides +
# FastRelax govern the bond itself.
#
# Success, route-aware:
#   amide  closure O...N 2.15-2.35 A
#   ss     S-S 1.95-2.15 A and |CB-SG-SG-CB| in 60-120 deg
#   both   zero non-proline cis and every bond within 30 deg of planar
set -uo pipefail

ml purge
ml biology rosetta/3.14

CYC=/scratch/users/bvbell/xlink/cyclize
WORK=$CYC/scripts/geom_test
mkdir -p "$WORK" "$CYC/logs"
cd "$CYC"

ARMS=(base omega loop7)
TARGETS=(apj uor)
SCAFFOLDS=(short medium)
ROUTES=(amide ssapp,ssmut)

ai=$(( SLURM_ARRAY_TASK_ID % 3 ))
ti=$(( (SLURM_ARRAY_TASK_ID / 3) % 2 ))
si=$(( (SLURM_ARRAY_TASK_ID / 6) % 2 ))
ri=$(( SLURM_ARRAY_TASK_ID / 12 ))
ARM=${ARMS[$ai]}; TGT=${TARGETS[$ti]}; SCAF=${SCAFFOLDS[$si]}; RTE=${ROUTES[$ri]}

echo "=== $TGT $SCAF arm=$ARM routes=$RTE on $(hostname) $(date) ==="
python3 "$CYC/scripts/run_geometry_test.py" \
    --target "$TGT" --scaffold "$SCAF" --arm "$ARM" --routes "$RTE" \
    --nstruct 5 --out "$WORK"
echo "=== exit $? $(date) ==="
