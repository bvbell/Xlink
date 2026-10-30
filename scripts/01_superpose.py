"""
Superimpose CLDN6 onto the CLDN18.2 antigen (chain A) of 9v32 and report
superposition quality for three PyMOL methods (super / align / cealign).

Run with the PyMOL 2.5.3 module:
    ml viz pymol/2.5.3
    pymol -cq /scratch/users/bvbell/xlink/scripts/01_superpose.py

Writes:
    /scratch/users/bvbell/xlink/work/transform_super.txt   4x4 transform (row-major)
    /scratch/users/bvbell/xlink/work/align_pairs_super.txt  aligned CA residue pairs
    /scratch/users/bvbell/xlink/work/superpose_report.txt   RMSD table
"""
import os
import json
from pymol import cmd

BASE = "/scratch/users/bvbell/xlink"
WORK = os.path.join(BASE, "work")
os.makedirs(WORK, exist_ok=True)

TARGET_PDB = os.path.join(BASE, "9v32_complex.pdb")
MOBILE_PDB = os.path.join(BASE, "cldn6.pdb")

report = []


def log(s):
    print(s)
    report.append(s)


def fresh_load():
    cmd.delete("all")
    cmd.load(TARGET_PDB, "cplx")
    cmd.create("tgtA", "cplx and chain A and polymer")
    cmd.load(MOBILE_PDB, "mob")
    cmd.delete("cplx")


# ---------------------------------------------------------------- methods
results = {}

# --- super (all-atom, default) ---
fresh_load()
r = cmd.super("mob", "tgtA", object="aln_super")
log("super (all-atom, default):")
log("  RMSD after refinement : %.3f A over %d atoms" % (r[0], r[1]))
log("  refinement cycles     : %d" % r[2])
log("  RMSD before refinement: %.3f A over %d atoms" % (r[3], r[4]))
log("  raw alignment score   : %.2f ; residues aligned: %d" % (r[5], r[6]))
results["super_all"] = r
M_super_all = cmd.get_object_matrix("mob")

# raw alignment pairs (CA level) from the super alignment object
raw = cmd.get_raw_alignment("aln_super")
model_names = {}
pairs = []
for idxs in raw:
    entry = {}
    for obj, idx in idxs:
        if obj not in model_names:
            model_names[obj] = cmd.get_model(obj).atom
        a = model_names[obj][idx - 1]
        entry[obj] = a
    if "mob" in entry and "tgtA" in entry:
        am, at = entry["mob"], entry["tgtA"]
        pairs.append((int(at.resi), at.resn, int(am.resi), am.resn, at.name, am.name))

with open(os.path.join(WORK, "align_pairs_super.txt"), "w") as fh:
    for p in pairs:
        fh.write("%d %s %d %s %s %s\n" % p)

# residue-level identity over aligned region (unique residue pairs)
respairs = sorted(set((p[0], p[1], p[2], p[3]) for p in pairs))
nid = sum(1 for p in respairs if p[1] == p[3])
log("  aligned residue pairs : %d ; identical: %d ; seq id = %.1f%%"
    % (len(respairs), nid, 100.0 * nid / max(1, len(respairs))))

with open(os.path.join(WORK, "respairs_super.txt"), "w") as fh:
    for p in respairs:
        fh.write("%d %s %d %s\n" % p)

with open(os.path.join(WORK, "transform_super.txt"), "w") as fh:
    fh.write(" ".join("%.10f" % v for v in M_super_all))

# --- super, CA only ---
fresh_load()
r = cmd.super("mob and name CA", "tgtA and name CA", object="aln_super_ca")
log("")
log("super (CA only):")
log("  RMSD after refinement : %.3f A over %d CA pairs (%d cycles)" % (r[0], r[1], r[2]))
log("  RMSD before refinement: %.3f A over %d CA pairs" % (r[3], r[4]))
results["super_ca"] = r
M_super_ca = cmd.get_object_matrix("mob")
with open(os.path.join(WORK, "transform_super_ca.txt"), "w") as fh:
    fh.write(" ".join("%.10f" % v for v in M_super_ca))

# --- align (sequence-driven) ---
fresh_load()
r = cmd.align("mob", "tgtA", object="aln_align")
log("")
log("align (all-atom, sequence-driven):")
log("  RMSD after refinement : %.3f A over %d atoms (%d cycles)" % (r[0], r[1], r[2]))
log("  RMSD before refinement: %.3f A over %d atoms" % (r[3], r[4]))
log("  raw alignment score   : %.2f ; residues aligned: %d" % (r[5], r[6]))
results["align_all"] = r
M_align = cmd.get_object_matrix("mob")
with open(os.path.join(WORK, "transform_align.txt"), "w") as fh:
    fh.write(" ".join("%.10f" % v for v in M_align))

# --- cealign ---
fresh_load()
r = cmd.cealign("tgtA", "mob", object="aln_ce")
log("")
log("cealign:")
log("  RMSD : %.3f A over %d aligned positions" % (r["RMSD"], r["alignment_length"]))
results["cealign"] = {"RMSD": r["RMSD"], "n": r["alignment_length"]}
M_ce = cmd.get_object_matrix("mob")
with open(os.path.join(WORK, "transform_cealign.txt"), "w") as fh:
    fh.write(" ".join("%.10f" % v for v in M_ce))

# ------------------------------------------- how different are the methods?
fresh_load()
cmd.super("mob", "tgtA")
cmd.create("mob_super", "mob")
fresh_load_names = None
cmd.delete("mob")
cmd.load(MOBILE_PDB, "mob2")
cmd.align("mob2", "tgtA")
cmd.create("mob_align", "mob2")
cmd.delete("mob2")
cmd.load(MOBILE_PDB, "mob3")
cmd.cealign("tgtA", "mob3")
cmd.create("mob_ce", "mob3")
cmd.delete("mob3")

log("")
log("Pairwise CA RMSD between the CLDN6 placements produced by each method")
log("(same atoms, no further fitting -- how much the answer moves):")
for a, b in [("mob_super", "mob_align"), ("mob_super", "mob_ce"), ("mob_align", "mob_ce")]:
    d = cmd.rms_cur("%s and name CA" % a, "%s and name CA" % b, matchmaker=-1)
    log("  %-10s vs %-10s : %.3f A" % (a.replace("mob_", ""), b.replace("mob_", ""), d))

with open(os.path.join(WORK, "superpose_report.txt"), "w") as fh:
    fh.write("\n".join(report) + "\n")

print("\nDONE")
