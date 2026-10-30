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

  **The scan finished 2026-10-30** (array 50209855, 595 trajectories; analysis 50209859). It is a
  much weaker interface than CLDN18.2's: **one hotspot above 1 REU** (`L_Y38A`, 2.94 REU GAM), 5.3
  REU over all 17 positions against 9.7 over CLDN18.2's 16. On the 12 shared positions CLDN6 sums
  3.9 against 7.4, and the loss is specific -- `H_N103A` 1.62 -> 0.31, `L_K36A` 1.30 -> 0.00,
  `H_W100A` 1.06 -> -0.10, while `L_Y38A` alone *gains* (1.90 -> 2.94). The five chimera-only
  contacts recover only ~1.4 REU between them.

  That is a partial answer to the open question above, in the direction of "real": the chimera does
  not merely clash, it loses four of CLDN18.2's five hotspots while keeping one. Full numbers and
  the matched table are in `../STATUS.md`.

  **The mechanism, worked out 2026-10-30.** This is a clean same-antibody comparison: the Fab
  coordinates are *identical* between the two complexes (1588 H + 1661 L atoms; the relax moved the
  Fab only 0.32 A CA RMSD), so every ddG difference is attributable to the antigen.

  **The 9V32 epitope is the least conserved part of the protein.** Structurally aligning the two
  antigens (nearest-CA, 2.5 A) gives 34% identity globally but **1 of 9 aligned epitope residues
  conserved -- 11%**. The lone survivor is V54 -> V55. Four further CLDN18.2 epitope residues (V40,
  T41, A42, S57) have no partner within 2.5 A at all: the backbone itself diverges there.

  Every lost hotspot maps onto a specific substitution, and each loses its partner twice over -- by
  chemistry *and* by distance:

  | Fab position | ddG 18.2 -> 6   | partner in CLDN18.2      | in CLDN6     | min distance      |
  | ------------ | --------------- | ------------------------ | ------------ | ----------------- |
  | `L_K36`      | 1.30 -> **0.00**  | E62, NZ-OE1 salt bridge  | E -> **Q**63 | 3.70 -> **7.41 A** |
  | `H_N103`     | 1.62 -> **0.31**  | E56, N-OE1 H-bond        | E -> **Q**57 | 3.45 -> **6.44 A** |
  | `H_W100`     | 1.06 -> **-0.10** | F60, aromatic/backbone   | F -> **Q**61 | 4.12 -> **6.22 A** |
  | `L_Y100`     | 1.51 -> *dropped* | S57, V40, A42, T41       | unalignable  | --                |

  Three acidic/aromatic -> glutamine swaps (E56Q, F60Q, E62Q) wipe out three of CLDN18.2's five
  hotspots, and the partners do not merely change identity -- they recede to 6-7 A, well past
  contact. `L_Y38` survives and *gains* (1.90 -> 2.94) because it picks up four contacts in CLDN6
  against two in CLDN18.2.

  **The dead heavy chain is an artifact of the construction, not a finding.** Contacts under 4.5 A
  go heavy 58 -> 48 but light **72 -> 172**: the chimera levers the heavy chain off the antigen while
  jamming the light chain into it. The light chain's new contacts buy almost nothing (`L_L29` 0.86,
  `L_Q27` 0.57, `L_S28` 0.04, `L_L30` -0.06, `L_N37` 0.01) -- nonspecific packing, not recognition.
  This resolves the question flagged above; no further check of heavy-chain placement is needed.

  **Interpretive limit, which bounds everything in this subsection.** The chimera holds the Fab in
  its CLDN18.2 pose, so the scan measures what 9V32 would get *if it approached CLDN6 exactly as it
  approaches CLDN18.2*. A real CLDN6 complex would re-dock. This is therefore strong evidence for
  **why** the Fab is CLDN18.2-specific, and a sound lower bound on cross-reactivity, but it is **not
  a prediction of a CLDN6 binding mode** and the 5.3 REU total must not be read as "CLDN6 affinity".
- **`uor_nb/`** is a Coot-built model of a nanobody on the human mu-opioid receptor (OPRM1 68-353; chain B
  numbering is OPRM1 minus 65). No gaps, breaks, altlocs or missing atoms.

  > **This scan's results are not currently trustworthy (found 2026-10-30).** Nanobody residues
  > **102-107 are D-amino acids** in `../uor_complex.pdb` -- six consecutive residues at the CDR3
  > tip, verified by the improper dihedral N-CA-C-CB (-50 to -92 deg against +57.6 +- 1.3 for the
  > rest of the chain; 0 inverted in the receptor, in the apj nanobody, and in experimental 9v32).
  > Bond lengths are normal, so it is a mirror-imaged backbone, almost certainly a local Coot
  > artifact -- the source filename contains "flip" twice and the prep report checked occupancies,
  > altlocs, missing atoms and bond lengths but never chirality.
  >
  > **Y102 is one of the six.** Since this scan is what identified the uOR hotspots and those
  > hotspots are what selected the cyclize motif, the ddG values here were measured on D residues
  > scored against L Ramachandran and rotamer terms. **Re-run this scan on the relaxed model**
  > (`../uor_complex_relaxed.pdb`, see `../STATUS.md` "The 2026-10-30 rebuild") before relying on any uOR hotspot, and treat the "8 hotspots > 1 REU, 13.6 REU total"
  > line in `../STATUS.md` as provisional. See `../cyclize/README.md`, "Actual root cause".
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

Chain breaks are reported in `chain_breaks.txt` but kept as Rosetta reads them (bonded neighbors). The
two nanobody-GPCR targets have gaps far from the interface and were run that way.

## Layout

- `flexddg/`: `config`, `pdbprep`, `interface`, `inputs` (setup), `runner` (one trajectory), `slurm`,
  `analysis`, `plots`, `smoke`, `__main__` (CLI). `vendor/` and `data/ddG-backrub.xml` are the flex ddG
  tutorial's files, unchanged.
- `tests/test_flexddg.py`: `python3 tests/test_flexddg.py` (~1.5 min). Golden files in `tests/golden/`
  pin the original TROP2 setup and summary and the inputs of the first CLDN run.
- `pipeline/`: the previous scripts. Array job 49679178 (CLDN), the last reader, completed long ago, so this directory can now be deleted.
- `trop2_9pi9/old_scripts/`: the original one-off TROP2 scripts. TROP2's original outputs keep their names
  (`fab_complete-*.csv`, `fab_alascan_summary.csv`, `plots/fab_alascan_*`); `alascan_summary.csv` etc. are
  the same numbers regenerated by this package.
- Python: `/home/groups/possu/miniconda/bin/python3` (3.8; numpy 1.24, pandas 1.5, matplotlib 3.7), set in
  `bin/flexddg`; override with `FLEXDDG_PYTHON`.
