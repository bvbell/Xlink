# Xlink

Antibody epitope work on the Sherlock cluster, in `$SCRATCH/xlink`. Two projects over the same Fab:

- **`flexddg/`** — flex ddG alanine scans of antibody-antigen interfaces (Rosetta). A reusable package
  plus one directory per target: TROP2 (9pi9) and CLDN18.2 (9v32). See `flexddg/README.md`.
- **`motif_scaffold/`** — RFdiffusion scaffolding of the 9pi9 Fab hotspots onto short peptides, using the
  hotspots the TROP2 scan found.

Structures at the top level are the inputs both use: `9pi9_complex.pdb` (Fab + TROP2),
`9v32_complex.pdb` (Fab + CLDN18.2), `cldn6.pdb` (CLDN6 monomer).

## What is not in this repo

Run outputs are excluded (see `.gitignore`) because they are gigabytes: `flexddg/*/output/`
(~4.3 GB of TROP2 trajectories, plus CLDN as it runs) and `motif_scaffold/runs*/` (~3.8 GB of
RFdiffusion output). **`$SCRATCH` purges files untouched for 90 days and is not backed up, so it holds
the only copy of all of that.** Anything that must outlive this phase belongs on Oak or `$GROUP_HOME`.

Also excluded: `flex_ddG_tutorial/` (a clone of Kortemme-Lab/flex_ddG_tutorial; the files used here are
vendored under `flexddg/flexddg/`) and `flexddg/pipeline/` (superseded scripts, still read by the running
CLDN array job).
