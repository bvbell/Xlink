# Project status

Last updated 2026-10-29. Hand-maintained; `cyclize status <target>` and `flexddg status <target>` are
the live sources of truth for job state.

## Running right now

Arrays `50154360` (apj, 27 variants) and `50154369` (uor, 28 variants) are part way through, each with
a validation job chained `afterany` (`50154361`, `50154371`). 50 structures per variant, one Rosetta
process per array task.

Closure results so far, by route -- the pattern is clean and worth knowing before reading any numbers:

| route   | outcome                                                                      |
| ------- | ---------------------------------------------------------------------------- |
| `amide` | **closes reliably, 50/50 on every variant tried**, both targets                 |
| `ssmut` | **50/50** (uor/short, the only scaffold eligible)                               |
| `ssapp` | **closed 0/50 on first submission** -- fixed in `b612855`, 8/8 on retest        |

### Outstanding: ssapp variants need a retry

The first submission of `ssapp` sampled nothing on the two terminal Cys, so the sulfurs could never
meet (details in `cyclize/README.md`). Fixed, but the already-run tasks must be redone. `apj` inputs
were regenerated before its array started, so some of its tasks picked the fix up and some did not,
depending on when they launched.

**Once both arrays finish** (`cyclize status <target>` shows no active job), rerun the failures:

```bash
cyclize status apj                 # confirm the empty ones are ssapp and nothing else
cyclize submit apj --retry         # reruns empty/partial variants only
cyclize submit uor --retry
```

`submit` refuses while an array is active, so this cannot be done early. Expect roughly 10-12 empty
`ssapp` variants per target.

## Where each piece stands

**Done — flex ddG alanine scans** (`flexddg/`, see its README):

| target            | scan | notes                                                        |
| ----------------- | ---- | ------------------------------------------------------------ |
| `trop2_9pi9/`     | done | original TROP2 run, outputs keep their original names          |
| `cldn18_9v32/`    | done | committed in `655f87b`                                         |
| `apj_nb/`         | done | 10 hotspots > 1 REU, 26.2 REU total                            |
| `uor_nb/`         | done | 8 hotspots > 1 REU, 13.6 REU total; noisier (Coot-built model)  |
| `cldn6_9v32fab/`  | **not run** | set up but `analysis_output/` is empty; chimera input, see README caveats |

**Running — CDR3 macrocyclisation** (`cyclize/`, see its README). Six motifs, three closure routes,
55 variants. The design constraint that drives everything is the hairpin register **i + j = 209**;
all six scaffolds are on-register and therefore even-length. Amide and ssmut closure both work; ssapp
needs the retry described above.

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
  Only uOR/short qualifies for that route at all. `ssapp`, which builds its Cys on free backbone,
  reaches 2.04-2.43 A instead -- so the strain really is the pinning, not the chemistry.
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

Everything is in `$SCRATCH`, purged at 90 days of no modification, so the git remote
(`git@github.com:bvbell/Xlink.git`) is the durable copy.

Committed: the `cyclize` package and both target configs, the scaffolds, both complex PDBs, and the
`apj_nb`/`uor_nb` alascan summaries and plots. The ~6 GB of trajectories under `*/output/` is
gitignored and regenerable.

Still uncommitted and worth a decision: `prep/`, `scripts/`, `flexddg/cldn6_9v32fab/`, and several
loose PDB/CIF files at the repo root (`apj_complex.pdb` and `uor_complex.pdb` are now tracked, but
`cldn6_9v32fab_complex.pdb`, `cldn6_9v32fab_relaxed_complex.pdb`, the Coot `Rosetta_9_real_trunc_*`
model and `d006_lmpnn_fkd2_run7_af3 (1).cif` are not). `prep/` matters most: the two prep notes it
holds are what the flexddg README points at for input provenance.

Copy final designs and summaries to `$GROUP_HOME` or Oak.

Rosetta: use `ml rosetta/3.14` for anything you run; `/home/groups/possu/Rosetta` is an older checkout
kept for reading source, offline docs and integration tests.
