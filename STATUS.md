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

### The 1200 passers are 112 distinct conformations

Clustering every passing structure on N/CA/C backbone RMSD (`cyclize/scripts/degeneracy.py`,
plotted in `cyclize/scripts/degeneracy.png`) gives **112 distinct conformations at 1.0 A** -- apj 51,
uOR 61. Representatives are extracted to `cyclize/<target>/representatives/` by
`scripts/representatives.py`.

**The best-passing variants are the most collapsed, not the most productive.** apj
`short__amide__DALA-GLY` and `short__amide__DPRO-PRO` are 50/50 and still a single cluster at 0.25 A:
fifty trajectories, one backbone. Pooling across turn seeds is worse -- apj's 402 short structures
from nine seed/route combinations are **5** backbones, uOR's 122 medium are **2**. The pinned motif,
not the turn seed, is dictating the closure.

This reframes the "short wins 3x" result above: part of that margin is that short has essentially one
accessible basin, so every trajectory finds it and it scores 100%. More nstruct buys nothing, and the
effective sample size behind the medium/long numbers is far smaller than 50 -- `ASN-GLY`'s "sharpest
scaffold split" rests on 2-4 structures.

### What to do next

1. **The i+1 turn sweep is complete** (2026-10-29; arrays 50164053 apj / 50164055 uor, in
   `cyclize/apj_i1/` and `cyclize/uor_i1/`). Sequence was the only axis left once sampling saturated,
   and the one position the closure chemistry frees is i+1 of the amide turn. 21 legal identities x 3
   scaffolds x 2 targets = 126 variants, all 63/63 done on both targets, **2483 passing structures**
   (apj 1403, uOR 1080). See `cyclize/README.md` for why this is enumeration rather than Rosetta
   design, and why `PHI_POS` now admits every D residue.

   **D-Ala is not special.** On apj/short, `ASN`, `ASP`, `DALA`, `DASN`, `DASP` and `DGLN` all close
   50/50 and all pass 50/50 — D-Ala's apparent win in the original five-seed sweep was an artifact of
   only five seeds being tried, not a property of D-Ala. That is the result the sweep was run to get,
   and it reopens i+1 as a position that can be chosen on chemistry rather than on buildability.
   Scaffold ordering is unchanged: apj short 850/1050 passing against long 340 and medium 213; uOR
   short 549 against medium 299 and long 232.

   **These are raw pass counts and have not been de-duplicated yet.** The 55-variant grid collapsed
   5:1 under clustering and this one will too, so the next step is `representatives.py` then
   `rescore.py` and `interface.py` on `apj_i1`/`uor_i1`. No conclusion about *which* identity is best
   should be drawn from the counts above.
2. **Score the sweep** on total `rama_prepro`, loop-only `rama_prepro`, and InterfaceAnalyzer
   `dG_separated`. The mechanics are now verified against rosetta/3.14 and written up under "Scoring
   a closed design" in `cyclize/README.md`. The headline: **reading a cyclic pose naively gives the
   wrong `rama_prepro`, silently**, and wrong at exactly the turn residues -- on
   `short__amide__ASN-GLY` the naive total is 22.8 against a true 29.4, with the error landing on two
   of the three GenKIC loop residues. `<PeptideCyclizeMover/>` fixes it for amide designs; for the
   disulfide routes the naive read is already correct and applying that mover instead *creates* a
   bogus number. Per-residue sums need `residue_energy_breakdown`, not `PerResidueEnergyMetric`
   (which splits each pair 50/50 across its two residues). Interface scoring needs a two-chain
   complex built via `validate.superimpose_on_motif`, and must repack identically across all designs
   or it ranks clash severity rather than binding.

   **The `rama_prepro` half is done** -- `cyclize/scripts/rescore.py`, route-aware, writes
   `<target>/rescore/rescore.csv` (112 representatives in ~25 s/target). It reproduces the reference
   measurement exactly, and the run-time-versus-clean comparison separates the routes the way the
   `ssapp` bug predicts: `ssapp` +10.61 REU mean over a tight +9.4 to +12.2 band (the phantom),
   `amide` +4.06 over a loose +0.8 to +7.7 (unexplained, probably PDB coordinate rounding through a
   steep potential, and harmless to ranking since every structure takes the identical path).
   Rescored numbers are therefore **not** comparable to `score.sc` numbers.

   **The interface half is done too** -- `cyclize/scripts/interface.py` builds the two-chain complex
   (the designs are peptide-only; chain B is never in the pose) and scores `dG_separated` and `ddg`
   with an identical sidechain repack on every design, ~3.6 s/structure. Two results:

   - **apj looks sane and reorders the ranking.** Best is `short__amide__DPRO-GLY` at **-25.2 REU**
     over 1429 A^2 buried, then `short__amide__ASN-GLY` at -21.7. Note that
     `short__amide__DALA-GLY` -- the 50/50 pass-rate winner with the 0.09 A motif RMSD -- comes in
     at only -15.4. **Buildability and interface score disagree**, which is the first evidence that
     the pass-rate ranking in the table above was measuring the wrong thing.
   - **uOR numbers are unusable**, at +370 to +413 REU. That is the known 1.09 A Y102-K238 clash:
     Y102 is a pinned motif residue in all three uOR scaffolds and the contact is interface-side, so
     it does not cancel in `dG_separated`. The spread among the best few is the same size as the
     differences being ranked. **Fix the input before reading anything into uOR binding.**
3. **Fix the uOR input.** It carries a 1.09 A Y102-K238 clash (see `cyclize/README.md`), which makes
   fine distance filtering on that target untrustworthy. It also shows up in the energies: uOR's
   `rama_prepro` runs 33-37 against apj's 10-14, so **only within-target ranking is meaningful**.
4. Flex ddG remains undone and is no longer treated as a prerequisite -- the motifs come from real
   binders, so the sweep is ranked on torsion quality plus a fast interface score instead.

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
- **Every `ssapp` run carries a ~10 REU phantom energy on residue 1.** `build.xml` strips
  `LOWER_TERMINUS_VARIANT` from residue 1 for all routes, but `ssapp` closes through S-S and never
  bonds anything to that N, so Rosetta scores its rama at phi = 0 -- a forbidden region. The artifact
  was present during `FastRelax` in all 600 apj + 600 uOR `ssapp` structures, so it biased the
  sampling, not just the reported score. **Fixed 2026-10-29** -- the two movers are now guarded on
  the closure chemistry (`build.py`, emitted only when the ring closes through a peptide bond). The
  completed `ssapp` and `ssmut` runs were produced *before* the fix and still carry it; rerun them
  if their geometry matters. `scripts/rescore.py` reads the PDBs as-is and does not reproduce the
  artifact, so the rescored numbers are clean. Details in `cyclize/README.md`.
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
