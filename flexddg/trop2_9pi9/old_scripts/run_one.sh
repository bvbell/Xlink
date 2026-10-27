#!/bin/bash
# Run one flex ddG trajectory: run_one.sh <mutation_name> <nstruct_index>
# Production defaults from Barlow et al. 2018; override via env for testing.
set -euo pipefail
MUT=$1; I=$2
NTRIALS=${NTRIALS:-35000}
STRIDE=${STRIDE:-35000}
MAXMIN=${MAXMIN:-5000}
CONV=${CONV:-1.0}
BASE=$(cd "$(dirname "$0")" && pwd)
IN=$BASE/inputs/$MUT
OUT=$BASE/output/$MUT/$(printf %02d "$I")
[ -f "$OUT/ddG.db3" ] && grep -q "reported success" "$OUT/rosetta.out" 2>/dev/null && { echo "done: $OUT"; exit 0; }
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
rosetta_scripts.default.linuxgccrelease \
  -s "$IN/9pi9_clean.pdb" \
  -parser:protocol "$BASE/ddG-backrub.xml" \
  -parser:script_vars chainstomove=$(cat "$IN/chains_to_move.txt") \
    mutate_resfile_relpath="$IN/nataa_mutations.resfile" \
    number_backrub_trials=$NTRIALS max_minimization_iter=$MAXMIN \
    abs_score_convergence_thresh=$CONV backrub_trajectory_stride=$STRIDE \
  -restore_talaris_behavior -in:file:fullatom -ignore_unrecognized_res \
  -ignore_zero_occupancy false -ex1 -ex2 -constant_seed -jran $((1000 + I)) \
  > rosetta.out 2>&1
