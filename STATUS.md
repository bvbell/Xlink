# Project status

Last updated 2026-10-29. Hand-maintained; `cyclize status <target>` and `flexddg status <target>` are
the live sources of truth for job state.

## Results — all 55 variants complete

Both arrays and both retries finished; every variant closed 50/50. **2750 structures, 1200 passed
validation.** Numbers from `cyclize/<target>/analysis_output/summary.csv`.

| target | route   | variants | closed | passed | rate |
| ------ | ------- | -------- | ------ | ------ | ---- |
| apj    | `amide` | 15       | 750    | 403    | 54%  |
| apj    | `ssapp` | 12       | 600    | 232    | 39%  |
| uor    | `amide` | 15       | 750    | 318    | 42%  |
| uor    | `ssapp` | 12       | 600    | 197    | 33%  |
| uor    | `ssmut` | 1        | 50     | 50     | 100% |

**Closure is not the discriminator — every variant closed 50/50.** The scaffold is:

| scaffold | apj          | uor          |
| -------- | ------------ | ------------ |
| short    | **89%**      | **72%**      |
| medium   | 24%          | 27%          |
| long     | 28%          | 19%          |

The shortest motif wins on both targets by a factor of 3, which is the **opposite** of what hotspot
coverage predicts: short keeps only 60-67% of the alanine-scan signal but is far and away the most
buildable. Fewer released motif residues means less for the sampler to distort. Whether that trade is
worth it cannot be settled from these numbers -- it needs ddG on the survivors.

Best per target:

- **apj**: `short__amide__ASN-GLY`, `short__amide__DALA-GLY`, `short__amide__DPRO-GLY` -- all 50/50,
  best motif RMSD **0.09 A**. Four of the five turn seeds pass at 100% on short.
- **uor**: `short__ssmut__CYS-CYS` 50/50 at **0.00 A** (it reproduces the native backbone exactly, by
  construction), then `short__amide__DALA-GLY` 46/50 and `short__ssapp__DCYD-DCYD` 44/50.

`ASN-GLY` is the sharpest scaffold split on apj -- 50/50 on short, 4/50 on long, 2/50 on medium -- so
turn choice is strongly scaffold-dependent and the five-seed sweep earned its keep.

### What to do next

1. **Flex ddG the survivors.** Validation only checks geometry; nothing here is a binding number.
2. **Fix the uOR input.** It carries a 1.09 A Y102-K238 clash (see `cyclize/README.md`), which makes
   fine distance filtering on that target untrustworthy.
3. Consider pushing `short` harder -- more nstruct, and the unimplemented `design_i2`.

## Where each piece stands

**Done — flex ddG alanine scans** (`flexddg/`, see its README):

| target            | scan | notes                                                        |
| ----------------- | ---- | ------------------------------------------------------------ |
| `trop2_9pi9/`     | done | original TROP2 run, outputs keep their original names          |
| `cldn18_9v32/`    | done | committed in `655f87b`                                         |
| `apj_nb/`         | done | 10 hotspots > 1 REU, 26.2 REU total                            |
| `uor_nb/`         | done | 8 hotspots > 1 REU, 13.6 REU total; noisier (Coot-built model)  |
| `cldn6_9v32fab/`  | **not run** | set up but `analysis_output/` is empty; chimera input, see README caveats |

**Done — CDR3 macrocyclisation** (`cyclize/`, see its README). Six motifs, three closure routes,
55 variants, all complete; results above. The design constraint that drives everything is the hairpin
register **i + j = 209**; all six scaffolds are on-register and therefore even-length.

**Not started:**
- Flex ddG on the closed macrocycles -- the single biggest gap. Validation only filters on motif RMSD,
  seed handedness and receptor clash; **no number here is a binding energy.** Re-scoring the 1200
  survivors through `flexddg/` would give values directly comparable to the nanobody scans, and is the
  only way to test whether `short` winning on buildability survives contact with affinity.
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
- **Seed handedness costs real yield.** D-Pro cannot be a GenKIC pivot and keep its chirality, and
  moving the pivot off it stops closure entirely, so wrong-handed solutions are filtered after the fact
  rather than prevented. Measured across both targets, amide passes 42-54% of closed structures.
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
