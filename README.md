# Xlink

Binder epitope work on the Sherlock cluster, in `$SCRATCH/xlink`: find which residues carry a binding
interface, then rebuild those residues on something smaller than the parent antibody or nanobody.

- **`flexddg/`** — flex ddG alanine scans of protein-protein interfaces (Rosetta). A reusable package
  plus one directory per target: TROP2 (9pi9), CLDN18.2 (9v32), a CLDN6 chimera, and two
  nanobody-GPCR complexes. **All five scans are complete** as of 2026-10-30. See `flexddg/README.md`.
- **`motif_scaffold/`** — RFdiffusion scaffolding of the 9pi9 Fab hotspots onto short peptides, using the
  hotspots the TROP2 scan found.
- **`cyclize/`** — closes the nanobody CDR3 hairpins into cyclic peptides by amide or disulfide,
  holding the ddG-validated conformation rigid, for APJ and the mu-opioid receptor. The governing
  constraint is the hairpin H-bond register (`i + j = 209` in both), which is what forces an
  even-length motif and a 0- or 2-residue closure. See `cyclize/README.md`.

`STATUS.md` is the hand-maintained overview of what is running, what is done, and what is still open.
**Read it first** — as of 2026-10-30 the `cyclize` rebuild chain is in flight (relax -> re-excise
-> geometry test, jobs 50318433 / 50319161 / 50321818), and four findings invalidate earlier
numbers:

- **uOR residues 102-107 are D-amino acids** — six consecutive residues at the CDR3 tip of
  `uor_complex.pdb`, i.e. exactly the transplanted motif. No uOR number in this repo is currently
  trustworthy and the uOR work needs rebuilding from a repaired model.
- **The macrocycle amide bonds are not planar** — 89% of all 314 built complexes carry at least one
  peptide bond >30° from planar, and 618 cis bonds sit where no proline follows (vs 25 legitimate
  cis-Pro). With omega free the chains fold through themselves, which is the source of the
  self-intersection and the ~3300 REU of strain. This affects **both** GPCR targets, so uOR carries
  two independent defects and apj one. `validate.py` checks the design against the receptor but
  never against itself or its own geometry.
- The committed `cyclize/*/interface/interface.csv` scores were measuring clash rather than binding,
  and `dG_separated` is blind to intra-peptide energy by construction.
- The CLDN6 chimera turned out to have a one-residue interface.

**Both `cyclize` targets are therefore in a rebuild state, and the `flexddg` CLDN scans are the only
results currently standing unqualified.** The pre-2026-10-30 design set is archived at
`$OAK/bvbell_archive/xlink_cyclize_invalid_geometry_20261030.tar.gz` and documented in
`cyclize/ARCHIVE.md`; the rebuild is described in `STATUS.md`.

Structures at the top level are the inputs those use:

| file | contents |
| ---- | -------- |
| `9pi9_complex.pdb` | Fab + TROP2 |
| `9v32_complex.pdb` | Fab + CLDN18.2 |
| `cldn6.pdb` | CLDN6 monomer, residues 1-187 |
| `cldn6_9v32fab_complex.pdb` | CLDN6 superimposed into 9v32 in place of CLDN18.2, Fab unmoved (built by `scripts/01_superpose.py` + `02_build_chimera.py`). 56 heavy-atom pairs under 2.0 A — superseded as an input by the relaxed model below |
| `cldn6_9v32fab_relaxed_complex.pdb` | the above after constrained FastRelax, picked from 20 decoys by `scripts/05_pick_relaxed_chimera.py`. 0 pairs under 2.0 A. **This is what `flexddg/cldn6_9v32fab/` actually scans** |
| `uor_complex.pdb` | nanobody + human mu-opioid receptor (OPRM1 68-353), from `Rosetta_9_real_trunc_model_1_flip_fit_J126map-coot-0_flip.pdb` |
| `apj_complex.pdb` | d006 nanobody + APJ/rubredoxin chimera (PDB 5VBL construct), from `d006_lmpnn_fkd2_run7_af3 (1).cif` |

The receptor in each nanobody complex was identified by exact UniProt matching, **not** from the file
names: `Rosetta_9_...` is the mu-opioid receptor and `d006_...` is APJ, which is the reverse of what the
working labels suggested. `prep/` holds the structure-preparation reports and `scripts/` the code that
produced these files.

## What is not in this repo

Run outputs are excluded (see `.gitignore`) because they are gigabytes: `cyclize/*/output/`, `flexddg/*/output/`
(~4.3 GB of TROP2 trajectories plus the CLDN scans) and `motif_scaffold/runs*/` (~3.8 GB of
RFdiffusion output). **`$SCRATCH` purges files untouched for 90 days and is not backed up, so it holds
the only copy of all of that.** Anything that must outlive this phase belongs on Oak or `$GROUP_HOME`.
`motif_scaffold/runs_prod/` is the nearest deadline — last written 2026-10-22/23, so it purges around
**2027-01-20**, and it is not regenerable without re-running RFdiffusion.

Also excluded: `flex_ddG_tutorial/` (a clone of Kortemme-Lab/flex_ddG_tutorial; the files used here are
vendored under `flexddg/flexddg/`) and `flexddg/pipeline/` (superseded scripts; the CLDN array job
49679178 that still read them completed long ago, so this directory can now be deleted).
