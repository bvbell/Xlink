#!/bin/bash
set +u
BASE=/scratch/users/bvbell/xlink/motif_scaffold
PDB=$BASE/inputs/9pi9_motif_H3L3_pad.pdb
CKPT=/home/groups/possu/RFdiffusion/models/Complex_base_ckpt.pt
NDES=8
source /home/groups/possu/miniconda/etc/profile.d/conda.sh
conda activate SE3nv
export PYTHONPATH=/home/groups/possu/RFdiffusion:$PYTHONPATH
cd $BASE/logs
while IFS=$'\t' read -r NAME CONTIG; do
  [ -z "$NAME" ] && continue
  echo "=== $NAME : [E175-251/0 $CONTIG] ==="
  mkdir -p $BASE/runs/$NAME
  python /home/groups/possu/RFdiffusion/scripts/run_inference.py \
    inference.input_pdb=$PDB \
    inference.output_prefix=$BASE/runs/$NAME/$NAME \
    inference.ckpt_override_path=$CKPT \
    "contigmap.contigs=[E175-251/0 $CONTIG]" \
    'ppi.hotspot_res=[E183,E239,E243,E245]' \
    inference.num_designs=$NDES \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    > $BASE/logs/$NAME.log 2>&1
  RC=$?
  N=$(ls $BASE/runs/$NAME/*.pdb 2>/dev/null | wc -l)
  echo "    rc=$RC  pdbs=$N"
done < $BASE/sweep.tsv
echo "SWEEP COMPLETE"
