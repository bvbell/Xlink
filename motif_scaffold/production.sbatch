#!/bin/bash
#SBATCH --job-name=rfd_9pi9_prod
#SBATCH --array=0-7
#SBATCH --partition=possu,bioe
#SBATCH --gres=gpu:1
# SE3nv ships torch 1.13.1, which supports only up to sm_86. Restrict to
# A100 (GPU_CC:8.0) / A40 (GPU_CC:8.6); H100 (9.0) dies with SIGILL.
#SBATCH --constraint=GPU_CC:8.0|GPU_CC:8.6
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/bvbell/xlink/motif_scaffold/logs/prod_%A_%a.out
#SBATCH --error=/scratch/users/bvbell/xlink/motif_scaffold/logs/prod_%A_%a.err

# 512 designs at the sweep-optimal contig: register -1 on both segments,
# no flanks, 2-6 linker.  64 designs per array task, offset by design_startnum
# so tasks write disjoint filenames into one directory.
set +u
BASE=/scratch/users/bvbell/xlink/motif_scaffold
OUT=$BASE/runs_prod
NDES=64
START=$(( SLURM_ARRAY_TASK_ID * NDES ))

source /home/groups/possu/miniconda/etc/profile.d/conda.sh
conda activate SE3nv
export PYTHONPATH=/home/groups/possu/RFdiffusion:$PYTHONPATH

mkdir -p $OUT
cd $BASE/logs

echo "task $SLURM_ARRAY_TASK_ID: designs ${START}..$(( START + NDES - 1 )) on $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader

python /home/groups/possu/RFdiffusion/scripts/run_inference.py \
  inference.input_pdb=$BASE/inputs/9pi9_motif_H3L3_pad.pdb \
  inference.output_prefix=$OUT/prod \
  inference.ckpt_override_path=/home/groups/possu/RFdiffusion/models/Complex_base_ckpt.pt \
  'contigmap.contigs=[E175-251/0 A100-106/2-6/B90-94]' \
  'ppi.hotspot_res=[E183,E239,E243,E245]' \
  inference.num_designs=$NDES \
  inference.design_startnum=$START \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0

RC=$?
echo "task $SLURM_ARRAY_TASK_ID rc=$RC  total pdbs now: $(ls $OUT/*.pdb 2>/dev/null | wc -l)"
exit $RC   # propagate failure; a trailing echo would mask it as COMPLETED
