# Project status

Last updated 2026-10-29. Hand-maintained; `cyclize status <target>` and `flexddg status <target>` are
the live sources of truth for job state.

## Running right now

| job                     | what                                     | tasks | state   |
| ----------------------- | ---------------------------------------- | ----- | ------- |
| `50154360` cyclize_apj  | APJ CDR3 macrocycle closure, 27 variants  | 27    | PENDING |
| `50154361`              | APJ validation, `afterany:50154360`       | 1     | PENDING |
| `50154369` cyclize_uor  | uOR CDR3 macrocycle closure, 28 variants  | 28    | PENDING |
| `50154371`              | uOR validation, `afterany:50154369`       | 1     | PENDING |

All four were queued behind `reloop_screen2` (464 tasks), so expect a wait before the first results.
Each cyclize task runs 50 structures in one Rosetta process, ~40 min, 2 h limit.

Check with `cyclize status apj` / `cyclize status uor`. Results land in
`cyclize/<target>/analysis_output/{results,summary}.csv` once the validation job runs; it is chained
`afterany`, so partial results still get analysed.

## Where each piece stands

**Done — flex ddG alanine scans** (`flexddg/`, see its README):

| target            | scan | notes                                                        |
| ----------------- | ---- | ------------------------------------------------------------ |
| `trop2_9pi9/`     | done | original TROP2 run, outputs keep their original names          |
| `cldn18_9v32/`    | done | committed in `655f87b`                                         |
| `apj_nb/`         | done | 10 hotspots > 1 REU, 26.2 REU total                            |
| `uor_nb/`         | done | 8 hotspots > 1 REU, 13.6 REU total; noisier (Coot-built model)  |
| `cldn6_9v32fab/`  | **not run** | set up but `analysis_output/` is empty; chimera input, see README caveats |

**Submitted — CDR3 macrocyclisation** (`cyclize/`, see its README). Six motifs, three closure routes,
55 variants. The design constraint that drives everything is the hairpin register **i + j = 209**;
all six scaffolds are on-register and therefore even-length.

**Not started:**
- Flex ddG on the closed macrocycles. The current validation only filters on motif RMSD, seed
  handedness and receptor clash -- it says nothing about binding energy. Re-scoring survivors through
  `flexddg/` would give numbers directly comparable to the nanobody scans.
- `cldn6_9v32fab/` scan.
- Whether APJ's R27 + Y33 (5.6 REU of CDR1) contact the same receptor surface as the CDR3 tip. This
  decides how much affinity a CDR3-only macrocycle can possibly retain, and is still an open question.
- `design_i2`: designing i+2 of the turn instead of taking it from the seed pair.

## Open questions and known problems

- **Hotspot coverage is 60-78%.** No CDR3-only macrocycle recovers the CDR1/CDR2 contribution --
  APJ loses R27 + Y33 (5.6 REU), uOR loses E54 + L32 (3.0 REU). APJ is the worse case. Expect the
  macrocycles to be substantially weaker binders than the parent nanobodies regardless of how
  faithfully the motif is transplanted.
- **`ssmut` gives a strained disulfide** (1.70 A measured vs 2.05 ideal) because the motif is pinned.
  Only uOR/short qualifies for that route at all.
- **~30% of closed amide designs will fail on seed handedness.** D-Pro cannot be a GenKIC pivot and
  keep its chirality, and moving the pivot off it stops closure entirely, so wrong-handed solutions
  are filtered after the fact rather than prevented.
- **uOR inputs are a Coot-built model**, APJ is an AlphaFold3 model carrying the 5VBL construct's five
  mutations including W261K at the W6.48 toggle. Neither receptor has a ligand modelled. Details in
  `prep/nanobody_gpcr_prep.md`.

## Resuming from a fresh session

Everything needed is committed: the package, both `cyclize.json` files, the scaffolds, the complexes,
the alascan summaries, `variants.txt` and `setup.log`. `inputs/` and `output/` are gitignored but
`inputs/` is fully regenerable with `cyclize setup <target>` (deterministic). Start with
`cyclize/README.md`, then `cyclize status apj`.

```bash
cd cyclize && python3 tests/test_cyclize.py      # ~2 s, no Rosetta, confirms the package is sane
./bin/cyclize status apj                          # job state and per-variant closure counts
```

## Repo hygiene

Everything is in `$SCRATCH`, purged at 90 days of no modification. Nothing here is a sole copy of
anything that matters -- but the working tree has a lot of uncommitted content: `cyclize/`,
`apj_scaffolds/`, `uOR_scaffolds/`, `flexddg/{apj_nb,uor_nb,cldn6_9v32fab}/`, `prep/`, `scripts/`, and
several loose PDB/CIF files at the root. Commit and push before relying on any of it. Copy final
designs and summaries to `$GROUP_HOME` or Oak.

Rosetta: use `ml rosetta/3.14` for anything you run; `/home/groups/possu/Rosetta` is an older checkout
kept for reading source, offline docs and integration tests.
