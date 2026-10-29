#!/bin/bash
# Full distillation pass over every target: cluster the passing structures down to distinct
# conformations, then score each representative on rama_prepro and on the interface.
#
# Submitted rather than run interactively because the interface scoring is ~3.6 s/structure and the
# four targets together carry several hundred representatives.
#
#   sbatch -D /scratch/users/bvbell/xlink/cyclize scripts/distill.sbatch
#
# Order matters: rescore and interface both read <target>/representatives/representatives.csv, so
# representatives must run first. compare_prebug needs the reruns finished, which they are.
set -uo pipefail

ml purge
ml biology rosetta/3.14

cd /scratch/users/bvbell/xlink/cyclize
export PYTHONPATH=/scratch/users/bvbell/xlink/cyclize${PYTHONPATH:+:$PYTHONPATH}
PY=/home/groups/possu/miniconda/bin/python3
TARGETS="apj uor apj_i1 uor_i1"

echo "=== distill on $(hostname), $(date) ==="

# Did the rm_lower fix move geometry, or only the reported score? Pre-fix structures are in
# <target>/output_prebug/; this is the measurement the docs currently mark as an open inference.
echo; echo "### compare_prebug"
$PY scripts/compare_prebug.py apj uor

echo; echo "### degeneracy (apj, uor only -- the committed plot covers the original grid)"
$PY scripts/degeneracy.py apj uor --csv-dir=scripts

for t in $TARGETS; do
  echo; echo "### representatives $t"
  $PY scripts/representatives.py "$t" || echo "representatives FAILED for $t"
done

for t in $TARGETS; do
  echo; echo "### rescore $t"
  $PY scripts/rescore.py "$t" || echo "rescore FAILED for $t"
done

for t in $TARGETS; do
  echo; echo "### interface $t"
  $PY scripts/interface.py "$t" || echo "interface FAILED for $t"
done

echo; echo "=== done $(date) ==="
