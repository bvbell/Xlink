#!/bin/bash
#SBATCH --job-name=reexcise
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/users/bvbell/xlink/cyclize/logs/reexcise_%j.out
set -uo pipefail
cd /scratch/users/bvbell/xlink/cyclize
python3 scripts/reexcise_motifs.py --write
