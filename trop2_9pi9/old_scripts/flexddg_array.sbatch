#!/bin/bash
#SBATCH --job-name=flexddg_ala
#SBATCH --partition=normal,possu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-874%100
#SBATCH --output=logs/%A_%a.out
# Fab-side interface alanine scan: 25 mutations x 35 nstruct, one trajectory per task, 100 at a time.
set -uo pipefail
ml biology rosetta/3.14
cd "$SLURM_SUBMIT_DIR"
NSTRUCT=35
mapfile -t MUTS < mutations_fab.txt
MUT=${MUTS[$((SLURM_ARRAY_TASK_ID / NSTRUCT))]}
I=$(( SLURM_ARRAY_TASK_ID % NSTRUCT + 1 ))
echo "task $SLURM_ARRAY_TASK_ID: $MUT nstruct $I on $(hostname)"
./run_one.sh "$MUT" "$I"
