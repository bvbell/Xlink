#!/bin/bash
#SBATCH --job-name=relax_iface
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-8
#SBATCH --output=/scratch/users/bvbell/xlink/cyclize/logs/relax_iface_%A_%a.out
#
# Does relaxing the complex rescue the uOR interface score?
#
# Why: 100% of the 68 built uOR design complexes carry a sub-2.0 A interchain pair (mean closest
# approach 1.19 A, worst 0.77 A), against 47% / 1.80 A for apj. The native uor_complex.pdb is clean
# (min 2.27 A, Y102 CD1 - K238 CD at 3.39 A), so the transplant creates the clash: the motif's Tyr
# lands 1.08 A from LYS238 CD. interface.py repacks sidechains but never relieves backbone, so a
# clash that repacking cannot rotamer its way out of shows up as a huge positive dG_separated.
#
# This scores each complex twice through the identical InterfaceAnalyzer protocol:
#   arm A  as-is                 -- reproduces the committed interface.csv number
#   arm B  constrained FastRelax -- CA-restrained (sd 0.5 A), then the same InterfaceAnalyzer
# If B goes strongly negative, the uOR numbers are a clash artifact and are recoverable. If B stays
# positive, the placement itself is wrong and the input model has to be replaced.
#
# apj tasks are controls: they should move much less, since they start nearly clean.
set -uo pipefail

ml purge
ml biology rosetta/3.14

CYC=/scratch/users/bvbell/xlink/cyclize
WORK=$CYC/scripts/relax_iface
mkdir -p "$WORK"

# task -> target, structure. uOR: 5 best by dG plus the ssmut outlier. apj: 3 best as controls.
STRUCTS=(
  "uor long__amide__DPRO-GLY__c00   amide"
  "uor long__amide__GLY-GLY__c00    amide"
  "uor long__amide__ASN-GLY__c00    amide"
  "uor long__ssapp__DCYD-DCYD__c00  ssapp"
  "uor medium__amide__ASN-GLY__c00  amide"
  "uor short__ssmut__CYS-CYS__c00   ssmut"
  "apj short__amide__DPRO-GLY__c00  amide"
  "apj short__ssapp__DCYD-DCYD__c01 ssapp"
  "apj short__ssapp__CYD-CYD__c00   ssapp"
)
read -r TARGET STRUCT ROUTE <<< "${STRUCTS[$SLURM_ARRAY_TASK_ID]}"

PDB=$CYC/$TARGET/interface/complexes/$STRUCT.pdb
if [ ! -f "$PDB" ]; then echo "MISSING $PDB"; exit 1; fi

# Receptor chain is B for both targets (cyclize.json); the design is chain P (interface.py).
RCHAIN=B

# amide closes through a peptide bond and must be declared cyclic before scoring; the disulfide
# routes must NOT get PeptideCyclizeMover -- same rule as rescore.py / interface.py.
if [ "$ROUTE" = "amide" ]; then
  CYCMOVER='    <PeptideCyclizeMover name="cyc" residue_selector="pep"/>'
  CYCSTEP='    <Add mover="cyc"/>'
else
  CYCMOVER=''
  CYCSTEP=''
fi

write_xml () {  # $1 = with_relax (0|1), $2 = outfile
  local relax_movers='' relax_steps=''
  if [ "$1" = "1" ]; then
    relax_movers='    <VirtualRoot name="vroot" removable="true"/>
    <AddConstraints name="addcst">
      <CoordinateConstraintGenerator name="coordgen" sd="0.5" ca_only="true"/>
    </AddConstraints>
    <FastRelax name="relax" scorefxn="sfxn_cst" repeats="1"/>
    <RemoveConstraints name="rmcst" constraint_generators="coordgen"/>'
    relax_steps='    <Add mover="vroot"/>
    <Add mover="addcst"/>
    <Add mover="relax"/>
    <Add mover="rmcst"/>'
  fi
  cat > "$2" <<XML
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn" weights="ref2015"/>
    <ScoreFunction name="sfxn_cst" weights="ref2015">
      <Reweight scoretype="coordinate_constraint" weight="1.0"/>
    </ScoreFunction>
  </SCOREFXNS>
  <RESIDUE_SELECTORS><Chain name="pep" chains="P"/></RESIDUE_SELECTORS>
  <MOVERS>
$CYCMOVER
$relax_movers
    <InterfaceAnalyzerMover name="ia" scorefxn="sfxn" interface="P_$RCHAIN"
                            pack_separated="true" pack_input="true"
                            packstat="false" interface_sc="true"
                            scorefile_reporting_prefix="if"/>
  </MOVERS>
  <FILTERS>
    <Ddg name="ddg" scorefxn="sfxn" jump="1" repeats="1" repack="1" confidence="0"/>
  </FILTERS>
  <PROTOCOLS>
$CYCSTEP
$relax_steps
    <Add mover="ia"/>
    <Add filter="ddg"/>
  </PROTOCOLS>
  <OUTPUT scorefxn="sfxn"/>
</ROSETTASCRIPTS>
XML
}

cd "$WORK"
for ARM in asis relax; do
  [ "$ARM" = "relax" ] && R=1 || R=0
  OUT=$WORK/${TARGET}__${STRUCT}__${ARM}
  mkdir -p "$OUT"
  write_xml "$R" "$OUT/protocol.xml"
  echo "=== $TARGET $STRUCT arm=$ARM (relax=$R) on $(hostname) $(date +%T) ==="
  rosetta_scripts.default.linuxgccrelease \
    -s "$PDB" \
    -parser:protocol "$OUT/protocol.xml" \
    -out:path:all "$OUT" \
    -out:file:scorefile score.sc \
    -nstruct 1 -overwrite -mute all \
    -ignore_unrecognized_res -ignore_zero_occupancy false \
    -constant_seed -jran $((7000 + SLURM_ARRAY_TASK_ID))
  echo "  exit $? ; score:"
  awk 'NR<=2 || /^SCORE:/' "$OUT/score.sc" 2>/dev/null | tail -2 | cut -c1-200
done

echo "=== done $(date) ==="
