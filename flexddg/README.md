# flex ddG alanine scans of protein-protein interfaces

`flexddg/` is the package; each target is a directory with a `target.json`. Everything lives in `$SCRATCH`,
which purges files not modified for 90 days: push this repo somewhere, and copy final summaries/plots to
`$GROUP_HOME` or Oak.

| dir               | complex                        | groups (chains)                     | moved    | scanned  |
| ----------------- | ------------------------------ | ----------------------------------- | -------- | -------- |
| `trop2_9pi9/`     | `../9pi9_complex.pdb`          | fab (A heavy, B light), antigen (E) | antigen  | fab      |
| `cldn18_9v32/`    | `../9v32_complex.pdb`          | fab (H heavy, L light), antigen (A) | antigen  | fab      |
| `cldn6_9v32fab/`  | `../cldn6_9v32fab_relaxed_complex.pdb` | fab (H heavy, L light), antigen (A) | antigen  | fab      |
| `uor_nb/`         | `../uor_complex.pdb`           | nanobody (A), receptor (B)          | receptor | nanobody |
| `apj_nb/`         | `../apj_complex.pdb`           | nanobody (A), receptor (B)          | receptor | nanobody |

`../flexddg_ala` is a symlink to `trop2_9pi9/` (`motif_scaffold/build_input*.py` use that path).

The last three are described in `../prep/cldn6_chimera_prep.md` and `../prep/nanobody_gpcr_prep.md`, which
record how each input structure was built and what is wrong with it. Read those before trusting a number:

- **`cldn6_9v32fab/`** is a chimera, not a measured complex: CLDN6 superimposed onto CLDN18.2 in 9v32
  (PyMOL `super`, 1.60 A over 763 atom pairs, 38% identity) with the Fab left where it was. The raw
  superposition had 56 heavy-atom pairs under 2.0 A against the Fab, so **the target reads the
  relaxed model, not the raw chimera**. Whether those clashes are a placement artifact or the actual
  reason this Fab is CLDN18.2-specific is an open question the ddG numbers alone will not settle.

  The relax is done (array 49878193, 2026-10-27, `cldn6_9v32fab/relax/relax.sbatch`): 20 constrained
  FastRelax decoys at two coordinate-constraint strengths, `sd05` (0.5 A) and `sd10` (1.0 A).
  `../scripts/05_pick_relaxed_chimera.py` picks among them — lowest `total_score` of the decoys with
  zero sub-2.0 A pairs — and writes `../cldn6_9v32fab_relaxed_complex.pdb` with hydrogens stripped.
  **All 20 cleared every sub-2.0 A pair** and both pre-existing chain-H gaps survived; the winner is
  `decoys/sd10/15` at -766.5 REU, 0.31 A CA RMSD from the superposition (0.67 A over the interface).
  Residual contacts are 0 / 1 / 6 at 2.0 / 2.5 / 3.0 A with the worst at 2.47 A (a Gln-Tyr hydrogen
  bond), against 0 / 1 / 3 for the original 9v32. The `sd10` set beats `sd05` by ~135 REU even after
  subtracting the `coordinate_constraint` term, so it is genuinely better relaxed rather than just
  less restrained.

  Setup on the relaxed model yields **17 Fab positions to scan** (33 interface residues, 4 skipped as
  Ala/Gly), which is not the same 17 as the prep note's "shared with CLDN18.2" set: relaxation drops
  `L_S99`, `L_Y100` and `L_F102` out of the 5 A interface and adds the chimera-only `L_Q27 L_S28
  L_L29 L_L30 L_N37`. For the matched comparison against the CLDN18.2 run, filter at analysis time
  rather than restricting the scan.
- **`uor_nb/`** is a Coot-built model of a nanobody on the human mu-opioid receptor (OPRM1 68-353; chain B
  numbering is OPRM1 minus 65). Clean input: no gaps, breaks, altlocs or missing atoms.
- **`apj_nb/`** is an AlphaFold3 model of the d006 nanobody on APJ. Chain B is the PDB **5VBL** construct,
  so it carries that construct's five mutations (V117A, T177N, **W261K** at the W6.48 toggle switch,
  C325L, C326M) and has *C. pasteurianum* rubredoxin replacing ICL3 at chain-B residues 213-266. The
  rubredoxin is 39 A from the nanobody and has no Fe (four bare thiolates adding score noise, no
  interface effect); it can be deleted if it ever causes trouble.

Both nanobodies bind the **extracellular** face with CDR3 in the orthosteric pocket, so no interface
residue is membrane-buried and scoring the complex in solution is defensible. Neither receptor has a
ligand modelled.

## Use

```bash
bin/flexddg setup  mytarget        # clean PDB, interface scan, inputs/, mutations*.txt, chain_breaks.txt
bin/flexddg smoke  mytarget --sbatch   # one 100-trial trajectory in $SCRATCH/flexddg_smoke/, then parse it
bin/flexddg submit mytarget        # queue only unfinished trajectories + a summarize/plot job after them
bin/flexddg status mytarget        # done / incomplete / missing per mutation, active jobs
bin/flexddg summarize mytarget     # analysis_output/{results,struct_scores,alascan_summary}.csv
bin/flexddg plot   mytarget        # plots/alascan_{ranked,trajectories,hotspot_terms}.{png,pdf}
bin/flexddg compare trop2_9pi9 cldn18_9v32 -o compare   # hotspots side by side + hotspots.csv
```

`submit` refuses while the target has active jobs (pending is not a reason to resubmit), and does not
rerun incomplete trajectories unless you pass `--retry-failed` after `status` has shown you why they failed.
Each submission is recorded in `tasks/submissions.tsv` with its task list; log `logs/<job>_<i>.out` ran
line `i` of that list.

## target.json

```json
{"title": "CLDN18.2 (9v32)", "pdb": "../../9v32_complex.pdb",
 "groups": {"fab": ["H", "L"], "antigen": ["A"]},
 "move": "antigen",
 "chain_labels": {"H": "Heavy chain", "L": "Light chain", "A": "CLDN18.2"},
 "scan": {"side": ["fab"]},
 "mutations": ["H:W33F", "H:W33A+L:Y100A"]}
```

Optional sections and their defaults (see `flexddg/config.py`; unknown keys are rejected):

- `scan`: `mode` interface|none, `side` (groups to mutate, default all), `cutoff` 5.0 A, `to` "A",
  `skip` ["ALA", "GLY", "PCA", "disulfide"]
- `mutations`: explicit mutants run in addition to the scan; `CHAIN:WT<num>[icode]MUT`, `+` joins positions.
  The wild type is checked against the PDB.
- `protocol`: nstruct 35, ntrials 35000, stride 35000, maxmin 5000, conv 1.0, extra_flags [] (Barlow 2018)
- `slurm`: partition normal,possu, time 04:00:00, mem 3G, max_concurrent 100, analysis_time, analysis_mem
- `analysis`: bootstrap 1000, seed 0, hotspot 1.0 kcal/mol, min_nstruct (default: nstruct)
- `rosetta_module`: "biology rosetta/3.14"

Chain breaks are reported in `chain_breaks.txt` but kept as Rosetta reads them (bonded neighbors). Both
current targets have gaps far from the interface and were run that way.

## Layout

- `flexddg/`: `config`, `pdbprep`, `interface`, `inputs` (setup), `runner` (one trajectory), `slurm`,
  `analysis`, `plots`, `smoke`, `__main__` (CLI). `vendor/` and `data/ddG-backrub.xml` are the flex ddG
  tutorial's files, unchanged.
- `tests/test_flexddg.py`: `python3 tests/test_flexddg.py` (~1.5 min). Golden files in `tests/golden/`
  pin the original TROP2 setup and summary and the inputs of the first CLDN run.
- `pipeline/`: the previous scripts, still read by array job 49679178 (CLDN). Delete once it has finished.
- `trop2_9pi9/old_scripts/`: the original one-off TROP2 scripts. TROP2's original outputs keep their names
  (`fab_complete-*.csv`, `fab_alascan_summary.csv`, `plots/fab_alascan_*`); `alascan_summary.csv` etc. are
  the same numbers regenerated by this package.
- Python: `/home/groups/possu/miniconda/bin/python3` (3.8; numpy 1.24, pandas 1.5, matplotlib 3.7), set in
  `bin/flexddg`; override with `FLEXDDG_PYTHON`.
