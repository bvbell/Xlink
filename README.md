# Xlink

Binder epitope work on the Sherlock cluster, in `$SCRATCH/xlink`: find which residues carry a binding
interface, then rebuild those residues on something smaller than the parent antibody or nanobody.

- **`flexddg/`** — flex ddG alanine scans of protein-protein interfaces (Rosetta). A reusable package
  plus one directory per target: TROP2 (9pi9), CLDN18.2 (9v32), a CLDN6 chimera, and two
  nanobody-GPCR complexes. See `flexddg/README.md`.
- **`motif_scaffold/`** — RFdiffusion scaffolding of the 9pi9 Fab hotspots onto short peptides, using the
  hotspots the TROP2 scan found.
- **`cyclize/`** — closes the nanobody CDR3 hairpins into cyclic peptides by amide or disulfide,
  holding the ddG-validated conformation rigid, for APJ and the mu-opioid receptor. The governing
  constraint is the hairpin H-bond register (`i + j = 209` in both), which is what forces an
  even-length motif and a 0- or 2-residue closure. See `cyclize/README.md`.

`STATUS.md` is the hand-maintained overview of what is running, what is done, and what is still open.

Structures at the top level are the inputs those use:

| file | contents |
| ---- | -------- |
| `9pi9_complex.pdb` | Fab + TROP2 |
| `9v32_complex.pdb` | Fab + CLDN18.2 |
| `cldn6.pdb` | CLDN6 monomer, residues 1-187 |
| `cldn6_9v32fab_complex.pdb` | CLDN6 superimposed into 9v32 in place of CLDN18.2, Fab unmoved (built by `scripts/01_superpose.py` + `02_build_chimera.py`) |
| `uor_complex.pdb` | nanobody + human mu-opioid receptor (OPRM1 68-353), from `Rosetta_9_real_trunc_model_1_flip_fit_J126map-coot-0_flip.pdb` |
| `apj_complex.pdb` | d006 nanobody + APJ/rubredoxin chimera (PDB 5VBL construct), from `d006_lmpnn_fkd2_run7_af3 (1).cif` |

The receptor in each nanobody complex was identified by exact UniProt matching, **not** from the file
names: `Rosetta_9_...` is the mu-opioid receptor and `d006_...` is APJ, which is the reverse of what the
working labels suggested. `prep/` holds the structure-preparation reports and `scripts/` the code that
produced these files.

## What is not in this repo

Run outputs are excluded (see `.gitignore`) because they are gigabytes: `cyclize/*/output/`, `flexddg/*/output/`
(~4.3 GB of TROP2 trajectories, plus CLDN as it runs) and `motif_scaffold/runs*/` (~3.8 GB of
RFdiffusion output). **`$SCRATCH` purges files untouched for 90 days and is not backed up, so it holds
the only copy of all of that.** Anything that must outlive this phase belongs on Oak or `$GROUP_HOME`.

Also excluded: `flex_ddG_tutorial/` (a clone of Kortemme-Lab/flex_ddG_tutorial; the files used here are
vendored under `flexddg/flexddg/`) and `flexddg/pipeline/` (superseded scripts, still read by the running
CLDN array job).
