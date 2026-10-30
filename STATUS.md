# Project status

Last updated 2026-10-30. Hand-maintained; `cyclize status <target>` and `flexddg status <target>` are
the live sources of truth for job state.

**In flight as of 2026-10-30:** the rebuild chain — relax inputs (array **50318433**) ->
re-excise motifs (**50319161**) -> geometry test (array **50321818**), each chained `afterok` on
the previous. Read them out with `sacct -j <id>` and see "The 2026-10-30 rebuild" below. The two
jobs outstanding on 10-29 both finished and are read out below: job 50216944 (uOR clash relief)
and the `cldn6_9v32fab` scan (array 50209855 + analysis 50209859).

> ### Read this first: uOR residues 102-107 are D-amino acids
>
> Found 2026-10-30. Six consecutive residues at the CDR3 tip of `uor_complex.pdb` chain A have
> inverted backbone chirality — **exactly the transplanted motif**, Y102 included. It is a local
> Coot artifact (bond lengths are normal; the source filename contains "flip" twice), and it
> propagated into all three uOR scaffolds, the `uor_nb` alanine scan that selected the motif, and
> every uOR design.
>
> **It is the single cause of every uOR anomaly recorded in this file** — the `rama_prepro` of
> 33-37, the native-pair control scoring as a non-binder, and the 68-of-68 interface clash. Earlier
> explanations of those three ("global strain", "a strained input model", "a sidechain rotamer
> problem in the designs") are **all superseded**; they are left in place below with corrections
> attached, because the reasoning trail matters, but do not act on them.
>
> **No uOR number anywhere in this file is currently trustworthy, and uOR work cannot be rescored
> into validity — it has to be rebuilt from a repaired model.** Evidence and the repair plan:
> `cyclize/README.md`, "Actual root cause: uOR residues 102-107 are D-amino acids".
>
> ### And read this second: the macrocycles have broken amide bonds
>
> Also found 2026-10-30, independent of the above, and it affects **both** GPCR targets. The rings
> close with correct bond *lengths* (1.31-1.34 Å) but **omega is not planar**. Across all 314
> built complexes in all four targets (314 PDBs on disk; 311 of them carry a scored
> `interface.csv` row, which is the denominator used for validity below):
>
> - **281 / 314 (89%)** have at least one peptide bond more than 30° from planar
> - **476 / 3716 bonds (13%)** are non-planar; worst deviation median 65°, **max 90°** — a 90°
>   amide has no partial double-bond character at all, i.e. it is not a peptide bond
> - **618 cis bonds that are not before proline**, against 25 legitimate cis-Pro. Non-proline cis
>   is ~0.03% of bonds in real proteins; here it is **17%**
>
> With omega free, the chain folds back in ways real peptides cannot — which is what produces the
> self-intersection (`CA8–CA10` at **1.87 Å**, worst heavy-atom pair **0.64 Å**) and the ~3300 REU
> of intra-peptide `fa_rep`. Self-intersection is the *symptom*; non-planar omega is the cause.
>
> **The non-planarity is introduced by the build, not inherited** (checked 2026-10-30). No bond in
> any input chain exceeds 30° from planar: apj nanobody tops out at 25.2°, uOR at 10.3°,
> experimental 9v32 at 21.1°. Those input maxima are mildly non-ideal pre-proline bonds and are
> normal. The ≥30° population is entirely the pipeline's doing, and it lands specifically **at the
> turn and the ring-closure seam** — in `apj_i1 long__amide__DALA-GLY__c01` the two broken bonds
> are `1 ALA→2 LYS` at 90.0° (adjacent to the seam) and `13 DAL→14 GLY` at 70.8° (in the turn),
> while every motif-interior bond between them matches the input exactly.
>
> **~~Correspondingly, the good news: some structures introduce none at all.~~ Retracted the same
> day — see below. That claim rested on `apj medium__amide__DPRO-GLY__c00`, whose closure carbonyl
> oxygen turns out to sit 0.38 Å from the next nitrogen. There is no existence proof of a clean
> closure.**
>
> ### And read this third: the amide closure itself is wrong
>
> Found by visual inspection 2026-10-30, after every automated check above had passed it. **The
> closure carbonyl oxygen is misplaced.** At the ring-closure bond of
> `apj medium__amide__DPRO-GLY__c00` the `O···N` distance is **0.38 Å** against an ideal
> **2.25 Å**, with an `O=C–N` angle of **16.3°** against 123°. The oxygen is sitting essentially
> on top of the next residue's nitrogen.
>
> *(Attribution corrected by audit: that structure was the geometry-ranked **inspection set's**
> top pick, not this file's dG-best design. The dG-best, `apj short__amide__DPRO-GLY`, has a
> correct closure at 2.23 Å — it fails on flanking omega instead, at −75° and +45°. 75 of 239
> amides do have a sound closure oxygen, which is why 0-of-311 rests on the combination of
> criteria rather than on this one.)*
>
> It is exclusively a closure defect — **every one of 2627 internal peptide bonds measures
> 2.23–2.29 Å** — and it affects **136 / 239 amide complexes (57%)**, median 1.85 Å, minimum
> 0.14 Å. It is present in the raw build output, not introduced by scoring.
>
> **And the closure is paid for by wrecking its neighbours.** Bucketing every bond by ring-distance
> from the closure:
>
> | distance | n | non-planar >30° | cis (non-Pro) |
> | --- | --- | --- | --- |
> | 0 (the closure bond) | 239 | 0% | 0% |
> | **1 (flanking)** | 478 | **47%** | **50%** |
> | 2 | 478 | 24% | 43% |
> | 3 | 478 | 1% | 0% |
> | 4+ | 1193 | 0% | 0% |
>
> `PeptideCyclizeMover` constrains the closure bond itself — hence the pristine row 0 — but nothing
> constrains the flanking bonds or the carbonyl oxygen, so GenKIC buys the closure by flipping the
> neighbours cis and twisting them out of plane.
>
> **Why no filter caught it:** `O` and `N` are a **1-3 (geminal) pair through `C`**, and Rosetta
> excludes 1-3 pairs from `fa_rep` — so a 0.14 Å oxygen–nitrogen overlap costs *literally nothing*
> in score (`corr(closure O···N, total_score) = +0.10`). The intra-peptide clash scan skipped the
> closure pair by construction. `validate.py` checks motif RMSD, seed phi and receptor distance,
> none of which involve the oxygen. The defect was invisible to every metric and visible by eye.
>
> ### The bottom line: 0 of 311 complexes are geometrically clean
>
> | amide (239) | surviving |
> | --- | --- |
> | closure `O···N` 2.15–2.35 Å | 75 |
> | + no non-proline cis near the closure | **1** |
> | + flanking bonds planar (<30°) | **0** |
>
> | disulfide (72) | surviving |
> | --- | --- |
> | no non-proline cis near the closure | 1 |
> | + flanking planar | 1 |
> | + no self-intersection (≥3.0 Å) | **0** |
>
> Every amide design has a broken closure oxygen, or a non-proline cis bond beside the closure, or
> a non-planar flanking bond — usually more than one. **The GenKIC closure protocol as configured
> does not produce chemically valid macrocycles**, and this cannot be fixed by filtering, because
> there is nothing to filter down to.
>
> **One correction to the self-intersection claim above, in the designs' favour.** The `intra_min`
> criterion of ≥3.4 Å did not distinguish atom types, and a 2.8–3.0 Å N···O contact is a hydrogen
> bond, not a clash. Of 143 apj structures the closest intra-peptide contact is: **80 genuinely
> clashing** below 2.5 Å (e.g. 1.72 Å O···C), 56 polar at 2.5–3.0 Å (hydrogen bonds, fine), and 10
> nonpolar at 2.5–3.0 Å. So roughly **56%** self-intersect, not 100% — "the best in the set is
> 2.99 Å against a ≥3.4 Å norm" overstated it.
>
> **The 0-of-311 conclusion does not depend on that threshold.** Re-checked on criteria that are
> unambiguous — correct closure oxygen, no non-proline cis beside the closure, planar flanking
> bonds, and (for the disulfides) no sub-2.5 Å contact — the counts are still **0 / 239 amide** and
> **0 / 72 disulfide**.
>
> **The scorefunction flagged it and nothing was listening**: build-time `score.sc` reads `omega`
> 34.7 and `rama_prepro` 40.9, but `fa_rep` 3278.8 swamps them in `total_score` and no filter ever
> looked. GenKIC runs `loop_bump_check` and `rama_prepro_check` but **has no omega filter**, and
> `FastRelax` cannot repair it because the MoveMap sets `bb=0` on motif residues, freezing their
> omega at whatever the closure produced.
>
> `dG_separated` then **cancels** intra-peptide energy entirely, so the defect carries no score
> penalty while the collapsed mass buries extra surface — the ranking was actively selecting for
> it. All seven top apj designs by dG are among the worst offenders.
>
> **So both targets need rebuilding. uOR carries two independent defects (D-residues *and* broken
> closure); apj carries one.** The CLDN flex ddG targets are unaffected by both and are the only
> results currently standing unqualified. See "Recommendations" immediately below.

## The 2026-10-30 rebuild

**The entire pre-2026-10-30 design set is archived and superseded.** 0 of 311 scored complexes
were chemically valid; the full diagnosis is in `cyclize/ARCHIVE.md`.

```
$OAK/bvbell_archive/xlink_cyclize_invalid_geometry_20261030.tar.gz   72 MB, 10274 files
```

One tarball, not loose files — OAK was at 93% of its inode quota and `$GROUP_HOME` at 99% of
volume. The originals stay under `cyclize/<target>/` until `$SCRATCH` purges them.

### The plan

1. **Relax both input complexes** (`cyclize/scripts/relax_inputs.sbatch`, array **50318433**).
   CA-restrained FastRelax, sd 0.5 Å, 3 repeats, nstruct 10 per target.
2. **Re-excise the motifs from the relaxed poses** (`cyclize/scripts/reexcise_motifs.py`, job
   **50319161**, chained `afterok`). Residue ranges are read from the existing scaffolds rather
   than hard-coded, so the new scaffolds cover exactly the same positions and the `i + j = 209`
   register check applies unchanged — apj short 101-108 / medium 100-109 / long 99-110, uOR short
   101-108 / medium 100-109 / long 98-111. Writes `*_complex_relaxed.pdb` and
   `{apj_scaffolds,uOR_scaffolds}_relaxed/`.
3. **Realign designs to those same relaxed models.** This is the point of the exercise: taking the
   motif source *and* the alignment target from one pose removes the D-vs-L mismatch by
   construction, because Rosetta's L rebuild is then baked into both sides.
4. **Geometry test before any production run** (`cyclize/scripts/geometry_test.sbatch`, array
   **50321818**, chained). 2 targets x 2 scaffolds x 3 arms x **2 route groups** = 24 tasks,
   nstruct 5. **Both closure chemistries are tested**, because both are independently broken:
   `ssapp` S-S distance runs a median 2.28 Å against an ideal 2.05 and as high as **5.82 Å**
   (not a bond at all), with |chi3| outside the 60-120° gauche window in 43 of 74; the one `ssmut`
   structure is 1.70 Å at chi3 164°. `omega_trans` pins peptide bonds and skips the S-S by design,
   so for the disulfide routes it fixes the *flanking* bonds while `ForceDisulfides` + `FastRelax`
   govern the bond itself.

### The three test arms

| arm | what it changes |
| --- | --- |
| `base` | relaxed input only, omega unconstrained — control, should reproduce the failure |
| `omega` | `protocol.omega_trans=True` — pin every sampled peptide bond trans |
| `loop7` | `omega_trans` + `target_loop=7` — release two more cold terminal motif residues |

Validity is **route-aware**: amide needs closure `O···N` 2.15-2.35 Å; disulfide needs S-S
1.95-2.15 Å with |CB-SG-SG-CB| in 60-120°. Both additionally need **zero cis peptide bonds at all**
— not merely no non-proline cis — and every bond within 30° of planar. Both input nanobodies are
100% trans, so the stricter rule forbids nothing that was native; it does mean a cis-Pro turn,
which `base` sometimes produced, is now a rejection.

**The `omega` arm is the actual fix.** GenKIC with CA pivots solves phi/psi at the pivots and
leaves omega free — `variants.TARGET_LOOP`'s own comment says so — and `PeptideCyclizeMover`
constrains only the closure bond. Because omega is *not* a solved DOF, a `set_dihedral` perturber
setting it to 180° holds through the closure solve. Implemented in `build.py` behind
`protocol.omega_trans` (default on), with test coverage in `tests/test_cyclize.py`.

**`loop7` is the fallback and only runs because it is cheap to test alongside.** It extends the
sampled loop by releasing two further *cold* terminal motif residues — `pick_release` never
releases a hotspot. Note what this does and does not do: released residues **keep their identity**
(nothing is designed at those positions, only their backbone is sampled), and the motif still gains
**only 2 residues**. Use it only if `omega` alone cannot close.

### Decision: the relaxed structures are the design basis (2026-10-30)

**Standing decision — the relaxed models are what this application designs against, and the
chirality repair is not a prerequisite.** Rationale: designing L-peptide macrocycles requires a
self-consistent L target, and relaxing the input produces exactly that. Rosetta was already
L-ifying the motif silently on every read; relaxing makes that explicit and puts the motif source
and the alignment target in the same pose, which is what removes the D-vs-L mismatch.

What this accepts, stated plainly so it is not rediscovered later: the relaxed uOR motif is
**Rosetta's conformation, not the experimental one**. The L rebuild moves motif sidechain tips up
to **6.70 Å** (Arg107) and **6.65 Å** (Phe104) while the backbone holds within 0.72 Å. So uOR
designs are internally comparable and chemically sound, but the uOR motif conformation is a model
prediction rather than a density-derived observation. **apj is unaffected** — it is clean on
chirality and its relaxed model is a regularised version of the same conformation.

**What this does *not* resolve, and what still needs doing:** the `flexddg/uor_nb/` alanine scan
that *selected* the uOR motif ran on the D-residue model. Choosing to design against relaxed
structures does not make that scan's hotspots correct — they were measured on inverted
stereocentres. Re-running `uor_nb` on the relaxed model is still required before any uOR result is
interpreted as hotspot-driven, and it may change which residues the motif should contain. This is
independent of the design-basis decision and is tracked separately below.

Settling the hand against the original J126 density remains worthwhile for the structural record,
but it is **no longer on the critical path** for design.

## Recommendations (2026-10-30)

**Headline: do not rebuild either target yet. Both rebuilds depend on an unanswered geometry
question, and running them first would reproduce the same defects at the same scale.**

### Replace the closure protocol. Filtering cannot save this one.

**0 of 311 complexes are geometrically clean** (table in the banner above), so there is no
surviving subset to filter down to and no existence proof that the current protocol can produce a
valid macrocycle. An earlier version of this section said the opposite, on the strength of one
structure whose closure oxygen had not been measured; that has been retracted.

The defect is specific and consistent: `PeptideCyclizeMover` constrains the closure `C–N` bond and
nothing else, so the carbonyl oxygen and the two flanking peptide bonds are unconstrained, and
GenKIC satisfies the ring by distorting exactly those. Three things are therefore unconstrained
that must not be — O placement, flanking omega, and flanking cis/trans.

**Recommended: move to `simple_cycpep_predict`.** It is Rosetta's supported application for
cyclic peptides, it handles omega, cis-proline and ring closure natively rather than leaving them
to a generic KIC call, and `rosetta/3.14` on Sherlock has it (`ml rosetta/3.14`). The existing
`build.py` machinery — motif pinning, route selection, seed enumeration — is still useful; what
should be replaced is the GenKIC block and the cyclisation, not the surrounding pipeline.

If the GenKIC route is kept instead, all of the following are required, not optional:

1. **Constrain the closure carbonyl oxygen**, not just `C–N`. This is the defect no score term can
   see (O and N are a 1-3 pair, excluded from `fa_rep`).
2. **Constrain omega on the flanking bonds** and filter any non-proline cis. GenKIC has no native
   omega filter — that is the hole all of this came through.
3. **Let `FastRelax` minimize omega near the closure** — the MoveMap currently freezes motif
   backbone (`bb="0"`), so even a correct constraint cannot act there.
4. **Re-measure yield honestly afterwards**, against the full criteria list, before any production
   run. Expect it to be low; the current protocol's effective yield of valid structures is zero.

### APJ — the better target, and where to spend effort

Recommended: **make apj the lead target.** It carries one defect rather than two, its input is
structurally sound (AlphaFold3 model, 0 inverted residues, clean scaffolds), and its hotspot signal
is roughly twice uOR's (10 hotspots / 26.2 REU against 8 / 13.6). Fix the closure pipeline here,
where a clean result is interpretable, and only then apply the validated pipeline to uOR.

Order: gate test above -> fix `validate.py` (omega gate, intra-peptide clash gate, stop skipping
motif residues in the receptor-clash test) -> re-run validation over existing output to size the
damage -> rebuild -> rescore with a native reference.

Standing caveats that do not block this: the apj input carries the 5VBL construct's five mutations
including **W261K at the W6.48 toggle**, has no ligand modelled, and the CDR1 contribution
(R27 + Y33, 5.6 REU) is unrecoverable by any CDR3-only macrocycle. apj is the better target, not a
clean one.

### uOR — recommend deprioritising

Recommended: **park uOR until the apj pipeline is proven.** The case against spending on it now:

- it carries **both** defects — the D-residue motif *and* the same broken closure;
- its input is the least reliable object in the project: a Coot-built model, no ligand, and an
  undetected chirality error in precisely the six residues the whole design depends on;
- the chirality repair is not a quick fix, and **a relax will not do it**. Tested 2026-10-30: a
  CA-restrained `FastRelax` on chain A drops `rama_prepro` 31.96 → 9.50 and leaves only 1 of 6
  D residues — but that is not minimization inverting stereocentres (impossible); it is **Rosetta
  assigning L residue types from the standard three-letter codes and rebuilding CB from ideal L
  geometry on read**. The price is sidechain tips moving up to **6.70 Å** (Arg107) and **6.65 Å**
  (Phe104) while the backbone holds within 0.72 Å. A relax silently picks L and cannot adjudicate
  the hand; **only the original J126 density can**, and whether it is available and resolves the
  question is itself unverified;
- **everything that selected uOR as a target is downstream of the error.** The `uor_nb` alanine
  scan ran on the D-residue model, so the 8 hotspots and the 13.6 REU that justified the work are
  not independent evidence. That has to be re-established before, not after, the designs.

If uOR is resumed, the order is: repair 102-107 -> re-run the `uor_nb` scan -> **re-decide whether
the target is worth it on the corrected hotspots** -> re-excise scaffolds -> rebuild. Note the
decision point in the middle; it is a real one.

### CLDN6 — recommend no-go on scaffolding

Separately settled: the epitope is 11% conserved, the interface is a one-residue pharmacophore, and
the apparent light-chain contacts are an artifact of holding the Fab in its CLDN18.2 pose. Good
evidence for *why* 9V32 is CLDN18.2-specific; a poor foundation for a binder. Details below.

## Results — the original 55-variant grid

> **Read the banner above before this section.** "Closed" and "passed validation" both turned out
> to mean less than they appear: the apj rings close with correct bond geometry while passing
> *through themselves*, and validation never checked for it. Every pass rate below is a pass rate
> against filters now known to be incomplete.

Every variant closed 50/50. **2750 structures, 1205 passed validation.** Numbers from
`cyclize/<target>/analysis_output/summary.csv`, regenerated 2026-10-29 after the disulfide routes
were rerun under the `rm_lower` fix — so the `ssapp`/`ssmut` rows differ slightly from earlier
versions of this file (apj `ssapp` was 232, uOR `ssapp` was 197; see Open questions).

A separate 126-variant i+1 sweep lives in `apj_i1/` and `uor_i1/` and is reported below.

| target | route   | variants | closed | passed | rate |
| ------ | ------- | -------- | ------ | ------ | ---- |
| apj    | `amide` | 15       | 750    | 403    | 54%  |
| apj    | `ssapp` | 12       | 600    | 240    | 40%  |
| uor    | `amide` | 15       | 750    | 318    | 42%  |
| uor    | `ssapp` | 12       | 600    | 194    | 32%  |
| uor    | `ssmut` | 1        | 50     | 50     | 100% |

**Closure is not the discriminator — every variant closed 50/50.** The scaffold is:

| scaffold | apj      | uor      |
| -------- | -------- | -------- |
| short    | **91%**  | **71%**  |
| medium   | 22%      | 27%      |
| long     | 30%      | 19%      |

The shortest motif wins on both targets by a factor of 3, the **opposite** of what hotspot coverage
predicts: short keeps only 60-67% of the alanine-scan signal but is far and away the most buildable.

**Two later results cut against reading much into this table, and both are below.** Pass rate
largely measures how collapsed a variant's conformational ensemble is, not how good it is — the 1205
passers are only **120** distinct backbones. And on apj the interface score reorders the winners
outright: `short__amide__DALA-GLY` tops this table at 50/50 but scores only **-13.7** REU, against
**-25.2** for `short__amide__DPRO-GLY`.

Best per target **by pass count**, which is the weakest of the three rankings now available:

- **apj**: `short__amide__ASN-GLY`, `short__amide__DALA-GLY`, `short__amide__DPRO-GLY` -- all 50/50,
  best motif RMSD **0.09 A**. Four of the five turn seeds pass at 100% on short.
- **uor**: `short__ssmut__CYS-CYS` 50/50 at **0.00 A** (it reproduces the native backbone exactly, by
  construction), then `short__amide__DALA-GLY` 46/50 and `short__ssapp__DCYD-DCYD` 44/50.

`ASN-GLY` looked like the sharpest scaffold split on apj -- 50/50 on short, 4/50 on long, 2/50 on
medium. Treat that with caution: the medium and long arms rest on 2-4 structures, which clustering
showed collapse to 1-2 distinct backbones, and the i+1 sweep has since found six identities tying
`ASN-GLY` on short. The five-seed grid was too small to support the conclusion drawn from it.

### The passers are far fewer distinct conformations

Clustering every passing structure on N/CA/C backbone RMSD (`cyclize/scripts/degeneracy.py`,
plotted in `cyclize/scripts/degeneracy.png`) gives **120 distinct conformations at 1.0 A** -- apj 52
of 643 passing, uOR 68 of 562. Representatives are extracted to `cyclize/<target>/representatives/`
by `scripts/representatives.py`. (Final numbers, regenerated by job 50173650 against the post-rerun
set of 1205; the earlier 112 / apj 51 / uOR 61 came from the pre-rerun 1200.)

**The best-passing variants are the most collapsed, not the most productive.** apj
`short__amide__DALA-GLY` and `short__amide__DPRO-PRO` are 50/50 and still a single cluster at 0.25 A:
fifty trajectories, one backbone. Pooling across turn seeds is worse -- apj's 402 short structures
from nine seed/route combinations are **5** backbones, uOR's 122 medium are **2**. The pinned motif,
not the turn seed, is dictating the closure.

This reframes the "short wins 3x" result above: part of that margin is that short has essentially one
accessible basin, so every trajectory finds it and it scores 100%. More nstruct buys nothing, and the
effective sample size behind the medium/long numbers is far smaller than 50 -- `ASN-GLY`'s "sharpest
scaffold split" rests on 2-4 structures.

> **Reinterpreted 2026-10-30, and it is a sharper reading than the one above.** Given that every
> apj structure is self-intersecting, "one accessible basin" may not be a property of the target at
> all -- it may be that **the only closures GenKIC can find are self-intersecting ones, and they
> all live in one basin.** The clustering measured real degeneracy, but the explanation on offer
> here (the pinned motif dictates the closure) assumed the closures were valid structures. If the
> ring geometry is infeasible for a 2-residue turn, the collapse is a symptom of that, not of
> conformational rigidity. Settle the feasibility question before citing this section's
> interpretation again; the cluster *counts* stand either way.

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

   **The de-duplication has now landed, and it was worse than the 5:1 the grid showed: i+1 identity
   does not change the backbone at all.** Clustering across variants of the same scaffold (job
   50173650 for `representatives`, then `degeneracy.py` over all four targets on 2026-10-29),
   apj/short is **20 i+1 identities, 850 passing structures, 2 backbones** at 1.0 A with 656 of them
   in one cluster. apj/long is 21 identities and 340 structures in **3** clusters (330 in one);
   uOR/medium is 20 identities and 299 structures in **1**; uOR/long 21 identities and 232 in **2**
   (229 in one). Only uOR/short retains any spread, at 8 clusters over 549 structures.

   So the "six identities pass 50/50" result **is** one backbone wearing six side chains, exactly as
   feared. The stronger and cleaner statement that replaces "D-Ala is not special" is that **i+1 is
   conformationally inert**: pass counts across i+1 identities carry no structural information, and
   the position should be chosen on chemistry and interface score alone. Note the i+1 sweep is
   amide-only, so part of the collapse against the grid's 6 clusters per scaffold is route
   homogeneity -- the `ssapp` route is what supplied most of the grid's diversity.

   Per-target within-variant clusters: apj_i1 1403 passing -> **91** distinct (28 singletons),
   uOR_i1 1080 -> **100** (22 singletons). Read `cyclize/logs/distill_50173650.out` and
   `cyclize/logs/degeneracy_all4.out`; `scripts/degeneracy_{variants,cross}.csv` now cover all four
   targets rather than just the grid.
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
   `<target>/rescore/rescore.csv` (~25 s/target over the representatives). It reproduces the reference
   measurement exactly, and the run-time-versus-clean comparison separates the routes the way the
   `ssapp` bug predicts: `ssapp` +10.61 REU mean over a tight +9.4 to +12.2 band (the phantom),
   `amide` +4.06 over a loose +0.8 to +7.7 (unexplained, probably PDB coordinate rounding through a
   steep potential, and harmless to ranking since every structure takes the identical path).
   Rescored numbers are therefore **not** comparable to `score.sc` numbers.

   **The interface half is done too** -- `cyclize/scripts/interface.py` builds the two-chain complex
   (the designs are peptide-only; chain B is never in the pose) and scores `dG_separated` and `ddg`
   with an identical sidechain repack on every design, ~3.6 s/structure. Two results:

   All four targets were regenerated by job 50173650 on 2026-10-29, so the `ssapp`/`ssmut` rows are
   no longer stale. Numbers below supersede every earlier version of this file.

   - **apj looks sane and reorders the ranking.** On the grid, best is `short__amide__DPRO-GLY` at
     **-25.2 REU** over 1429 A^2 buried, and the post-fix disulfides now take 2nd and 3rd --
     `short__ssapp__DCYD-DCYD` **-21.1** and `short__ssapp__CYD-CYD` **-19.1**, both at much better
     motif RMSD (0.21 and 0.11 A) than the amide winner's 0.65. `short__amide__ASN-GLY` is -17.9
     (was reported -21.7 pre-rerun) and `short__amide__DALA-GLY` -- the 50/50 pass-rate winner with
     the 0.09 A motif RMSD -- comes in at only **-13.7** (was -15.4). **Buildability and interface
     score disagree**, which is the first evidence that the pass-rate ranking in the table above was
     measuring the wrong thing.
   - **Across the i+1 sweep, the best interface score in the project is `apj_i1`
     `short__amide__DGLU-GLY` at -27.0 REU**, then `DASN-GLY` -25.7, `DASP-GLY` -25.1, `DGLN-GLY`
     -24.3. Since i+1 is conformationally inert (item 1), these four are the *same backbone* scored
     with four side chains -- which is exactly the regime where an interface score, rather than a
     pass count, is the right discriminator. Read them as a side-chain ranking on one scaffold.
   - **uOR numbers are unusable**, at +370 to +413 REU. **The clash explanation is back, and the
     retraction in `083402e` was measured against the wrong object** (corrected 2026-10-29). Both
     things are true at once:
     - the **native** `uor_complex.pdb` is clean -- 0 interchain pairs under 2.0 A, closest approach
       2.27 A, and Y102 CD1 - K238 CD really is 3.39 A apart, which is what `083402e` measured;
     - but the **built design complexes** are not. **All 68 of 68** carry a sub-2.0 A interchain
       pair, mean closest approach **1.19 A**, worst **0.77 A**, and the recurring offender is the
       transplanted motif Tyr at 1.08 A from LYS238 CD. apj, for contrast, is 26 of 55 (47%) with a
       mean of 1.80 A.

     So the **transplant creates the clash** -- it is not in the input pose and it is not global.
     `interface.py` repacks sidechains but never relieves backbone, so a clash no rotamer can escape
     lands as a large positive `dG_separated`. Measure it with the snippet under Open questions.

     **Resolved 2026-10-30 by job 50216944: the clash is fully recoverable, and the uOR designs are
     not non-binders.** A CA-restrained FastRelax (sd 0.5 A) before the identical InterfaceAnalyzer
     flips every uOR structure from strongly positive to normally negative, and clears the clash
     outright:

     | structure                        | as-is dG | relaxed dG | closest approach |
     | -------------------------------- | -------- | ---------- | ---------------- |
     | uor `long__amide__GLY-GLY`       | +369.2   | **-33.0**  | 0.77 -> 2.88 A   |
     | uor `long__ssapp__DCYD-DCYD`     | +401.1   | **-29.9**  | 1.45 -> 2.85 A   |
     | uor `medium__amide__ASN-GLY`     | +526.0   | **-31.2**  | 1.14 -> 2.45 A   |
     | uor `long__amide__ASN-GLY`       | +367.3   | **-24.1**  | 1.29 -> 2.73 A   |
     | uor `short__ssmut__CYS-CYS`      | +913.0   | **-20.8**  | 1.09 -> 2.64 A   |
     | uor `long__amide__DPRO-GLY`      | +367.6   | **-14.2**  | 1.15 -> 2.76 A   |
     | apj `short__ssapp__DCYD-DCYD`    | -19.4    | -44.5      | 2.49 -> 2.52 A   |
     | apj `short__amide__DPRO-GLY`     | -20.7    | -37.4      | 2.44 -> 2.57 A   |
     | apj `short__ssapp__CYD-CYD`      | -22.4    | -36.8      | 2.36 -> 2.51 A   |

     The as-is arm is validated: it reproduces `long__amide__DPRO-GLY__c00` at 367.635 exactly. **So
     the ~400 REU was a scoring-pipeline artifact, not the model** -- `interface.py` never relieved
     backbone, and one un-escapable clash dominated the whole number.

     Three caveats, none of which undo that conclusion but all of which bound it:

     - **The relax arm is not a neutral control.** The apj structures improve by 15-25 REU too, so
       relaxed and as-is numbers are not comparable and the whole set has to be rescored before
       anything is ranked. The relaxed ordering already differs: `short__ssapp__DCYD-DCYD` (-44.5)
       now beats `short__amide__DPRO-GLY` (-37.4), reversing the as-is result.
     - **The relief is bought by moving the peptide, and motif fidelity was not re-measured.** The
       receptor barely moves (CA RMSD 0.13-0.17 A) but the design does: 0.43-0.72 A on uOR and up to
       **1.26 A (max 2.88 A) on apj `short__amide__DPRO-GLY`**. Against original motif RMSDs of
       0.09-0.65 A, that is a shift comparable to or larger than the entire motif-fidelity budget.
       **Post-relax motif RMSD against the parent nanobody is the missing measurement** -- until it
       exists, we cannot say the relaxed poses still present the pharmacophore they were built for.
     - `dG_cross` reads 0.000 in every relax arm (the VirtualRoot changes the fold tree). Use
       `dG_separated` and `ddg`, which are consistent across both arms.

     **This also undercuts the native-pair "floor".** Relaxed uOR designs (-14 to -33) score *better*
     than the +34 to +45 floor the native-pair control reported -- so that control, run as-is, was
     itself measuring clash rather than strain. A residual uOR-vs-apj gap of roughly 10-25 REU
     survives under relax and is plausibly the real model-strain penalty, but **the floor is now
     unquantified and has to be re-measured under the relax protocol.** ~~The independent strain
     evidence is untouched: uOR `rama_prepro` still runs 33-37 against apj's 10-14.~~
     **Superseded later the same day:** that `rama_prepro` gap is the D residues being scored
     against L Ramachandran tables, not model strain — relaxing chain A takes it from 31.96 to
     9.50 without any refitting. It is not independent evidence of anything.

     *Superseded 2026-10-30, kept for the reasoning trail:* **The two findings compose rather than
     compete.** The native-pair control (uOR +34.4 / +41.7 / +45.2 by scaffold against apj's -5.7 /
     -12.9 / -19.1) stands and still says the Coot model is strained -- it sets a **floor**. What it
     explicitly did not explain was the designs landing ~330 REU below that floor, and the universal
     transplant clash is the leading candidate for exactly that residual. If it holds, uOR is
     **two** problems: a strained model setting the floor, and a transplant artifact adding the rest
     -- and only the first needs a new input structure.

     *What actually happened:* the clash half was confirmed, the floor half was not. The "two
     problems" framing survives only in weakened form -- the transplant clash is real and accounts
     for the ~330, but the floor it was supposedly stacked on top of was an artifact of the same
     as-is protocol, so its size is unknown. Re-measure before quoting +34 / +41 / +45 again.

     Caveat on the control: it is a real per-scaffold measurement but **no committed script
     reproduces it**, so it cannot be re-run after an input change. Rebuilding it as a `--native`
     flag on `interface.py` would make it usable as the regression check it was meant to be -- now
     doubly needed, since the control has to be redone under relax regardless.
3. **~~Replace the uOR input model~~ — deprioritised 2026-10-30; fix `interface.py` instead.**
   50216944 answered the sizing question and the answer was "most of it is the pipeline": the
   ~330 REU residual clears entirely under a restrained relax, so the model buys back only the
   residual 10-25 REU gap, not the 400. **The cheap fix comes first** — add a restrained-relax step
   to `interface.py` and rescore all four targets. Replacing the model stays on the list (the
   `rama_prepro` strain was then believed real and independent) but is no longer blocking any
   number. **That `rama_prepro` premise did not survive either** — see item 8: the 33-37 is the
   D residues scored against L tables, which the relax resolves (31.96 -> 9.50).

   The paragraph below is the pre-50216944 reasoning, kept because the parts about *interpretability*
   still hold: uOR absolute numbers from the committed `interface.csv` remain unusable, and the
   `flexddg/uor_nb/` alanine scan that chose the uOR motifs was run on this same model.

   The native uOR
   complex scores as a non-binder and uOR's `rama_prepro` runs 33-37 against apj's 10-14, so the
   Coot-built model is genuinely strained throughout and a local repair will not fix *that* part.
   Until it is replaced, **no uOR number is interpretable in absolute terms and only within-target
   ranking has any meaning** — and that includes the `flexddg/uor_nb/` alanine scan that selected the
   uOR motifs in the first place, which was run on this same model.

   What changed on 2026-10-29 is the **size** of what a new model would buy. The strained model sets
   a floor at +34 to +45 REU; the designs sit ~330 REU below it because of a transplant clash that is
   present in 68 of 68 built complexes and absent from the input pose (see item 2 and Open questions).
   A new input structure fixes the floor, not the residual. If 50216944 shows the residual clears
   under a restrained relax, then **most of the uOR gap is a scoring-pipeline bug in `interface.py`,
   not the model**, and the cheaper fix comes first.

   *(That conditional resolved on 2026-10-30: the residual did clear, so the "most of the uOR gap is
   a scoring-pipeline bug" branch is the one that happened. Note the "+34 to +45 floor" in this
   paragraph did not survive — it was itself an as-is measurement, i.e. clash, and the relaxed
   designs score below it. See item 2.)*
4. Flex ddG remains undone and is no longer treated as a prerequisite -- the motifs come from real
   binders, so the sweep is ranked on torsion quality plus a fast interface score instead.

## Where each piece stands

**Done — flex ddG alanine scans** (`flexddg/`, see its README):

| target            | scan | notes                                                        |
| ----------------- | ---- | ------------------------------------------------------------ |
| `trop2_9pi9/`     | done | original TROP2 run, outputs keep their original names          |
| `cldn18_9v32/`    | done | committed in `655f87b`                                         |
| `apj_nb/`         | done | 10 hotspots > 1 REU, 26.2 REU total                            |
| `uor_nb/`         | **redo** | 8 hotspots > 1 REU, 13.6 REU total — **provisional: run on a model whose residues 102-107 are D-amino acids, Y102 among them** |
| `cldn6_9v32fab/`  | done | 17 positions x 35 traj; finished 2026-10-30, **1 hotspot > 1 REU, 5.3 REU total** |

**`cldn6_9v32fab/` was unblocked on 2026-10-29.** The 20-decoy constrained relax (array 49878193) had
finished on 10-27 and the winning decoy had been picked on 10-28, but nothing consumed it: no scan was
ever submitted and `target.json` still pointed at the **unrelaxed** chimera. The target now reads
`../../cldn6_9v32fab_relaxed_complex.pdb` (= `relax/decoys/sd10/15`, the lowest-scoring of 20 decoys
that all cleared every sub-2.0 A pair; 56 -> 0 clashes under 2.0 A, worst contact now 2.47 A). `setup`
on the relaxed model gives **17 Fab positions** over 33 interface residues. The smoke test passed on
the relaxed input (job 50208459, `H_S99A`, one 100-trial trajectory, parsed clean), and the scan
ran as **array 50209855, 595 trajectories, 100 at a time, with analysis job 50209859 chained
behind it** (`--dependency=afterok`) — both completed 2026-10-30. Full provenance and the
shifted interface list are in `flexddg/README.md`. Note the position list is *not* the prep note's
"17 shared with CLDN18.2" set -- relaxation drops `L_S99`/`L_Y100`/`L_F102` and adds five
chimera-only light-chain contacts; do the matched comparison by filtering at analysis time.

**It finished 2026-10-30, and the chimera interface is much weaker and much more concentrated than
CLDN18.2's.** `analysis_output/alascan_summary.csv`, plots in `plots/`:

- **One hotspot above 1 REU: `L_Y38A` at 2.94 REU GAM** (3.34 talaris), `frac_pos` 1.0. Then a long
  tail -- `L_W56A` 0.88, `L_L29A` 0.86, `L_Q27A` 0.57 -- and nothing else above 0.32. Total GAM over
  all 17 positions is **5.3 REU**, against **9.7** for `cldn18_9v32/` over its 16.
- **The matched comparison is the real result.** On the 12 positions the two scans share, CLDN6 sums
  **3.9 REU** against CLDN18's **7.4**, and the loss is specific rather than uniform:

  | position  | cldn6 | cldn18 | delta     |
  | --------- | ----- | ------ | --------- |
  | `L_Y38A`  | 2.94  | 1.90   | **+1.04** |
  | `L_W56A`  | 0.88  | 0.41   | +0.47     |
  | `H_N103A` | 0.31  | 1.62   | **-1.31** |
  | `L_K36A`  | 0.00  | 1.30   | **-1.30** |
  | `H_W100A` | -0.10 | 1.06   | **-1.16** |
  | `L_Y98A`  | -0.35 | 0.33   | -0.68     |
  | `L_N34A`  | 0.09  | 0.48   | -0.39     |
  | `H_S99A`  | -0.01 | 0.37   | -0.38     |

  CLDN18's four other >1 REU hotspots all collapse on CLDN6 while `L_Y38` alone gets stronger. The
  five chimera-only contacts that relaxation added (`L_Q27`, `L_S28`, `L_L29`, `L_L30`, `L_N37`)
  recover only ~1.4 REU between them, so they do not make up the difference.
- **The epitope is the least conserved part of the protein, and that is the whole explanation**
  (worked out 2026-10-30). The Fab coordinates are *identical* between the two complexes (1588 H +
  1661 L atoms; relax moved the Fab 0.32 A CA RMSD), so every ddG difference is the antigen's doing.
  Structurally aligning the antigens gives 34% identity globally but **1 of 9 aligned epitope
  residues conserved -- 11%** (only V54 -> V55); four more CLDN18.2 epitope residues have no partner
  within 2.5 A at all. Each lost hotspot maps to a specific substitution and loses its partner
  twice, by chemistry and by distance:

  | Fab position | ddG 18.2 -> 6   | partner in CLDN18.2 | in CLDN6     | min distance       |
  | ------------ | --------------- | ------------------- | ------------ | ------------------ |
  | `L_K36`      | 1.30 -> **0.00**  | E62 salt bridge     | E -> **Q**63 | 3.70 -> **7.41 A** |
  | `H_N103`     | 1.62 -> **0.31**  | E56 H-bond          | E -> **Q**57 | 3.45 -> **6.44 A** |
  | `H_W100`     | 1.06 -> **-0.10** | F60 aromatic        | F -> **Q**61 | 4.12 -> **6.22 A** |
  | `L_Y100`     | 1.51 -> *dropped* | S57/V40/A42/T41     | unalignable  | --                 |

  Three acidic/aromatic -> Gln swaps wipe out three of five hotspots. `L_Y38` survives and *gains*
  because it picks up four contacts in CLDN6 against two in CLDN18.2.
- **The dead heavy chain is an artifact, not a finding** -- resolved, no further check needed.
  Contacts under 4.5 A go heavy 58 -> 48 but light **72 -> 172**: the chimera levers the heavy chain
  off the antigen and jams the light chain into it. The light chain's new contacts buy almost
  nothing (`L_L29` 0.86, `L_Q27` 0.57, rest ~0) -- nonspecific packing, not recognition.
- **Interpretive limit.** The chimera holds the Fab in its CLDN18.2 pose, so this measures what 9V32
  would get *if it approached CLDN6 exactly as it approaches CLDN18.2*; a real complex would re-dock.
  Strong evidence for **why** the Fab is CLDN18.2-specific and a sound lower bound on
  cross-reactivity, but **not** a prediction of a CLDN6 binding mode -- 5.3 REU is not "CLDN6
  affinity".
- Consequence for motif design: a CLDN6 binder built from this interface is **a one-residue
  pharmacophore** -- one tyrosine's worth of real signal plus a light chain pressed into a surface it
  does not read. That is a far worse starting point than APJ (10 hotspots, 26.2 REU) or even uOR
  (8, 13.6), and it should be weighed before any scaffolding effort is spent on it.

**Done — CDR3 macrocyclisation** (`cyclize/`, see its README). Six motifs, three closure routes,
55 variants, all complete; results above. The design constraint that drives everything is the hairpin
register **i + j = 209**; all six scaffolds are on-register and therefore even-length.

**Not started:**
- Flex ddG on the closed macrocycles. **Deliberately dropped as a gate** -- the motifs come from real
  binders, so the designs are ranked on torsion quality plus a fast interface score instead
  (`scripts/rescore.py`, `scripts/interface.py`). Flex ddG at ~19 CPU-hours per point remains the
  validation tool for a final shortlist, not a ranking tool.
- Whether APJ's R27 + Y33 (5.6 REU of CDR1) contact the same receptor surface as the CDR3 tip. This
  decides how much affinity a CDR3-only macrocycle can possibly retain, and is still an open question.
- `design_i2`: designing i+2 of the turn instead of taking it from the seed pair.

## Open questions and known problems

- **Hotspot coverage is 60-78%.** No CDR3-only macrocycle recovers the CDR1/CDR2 contribution --
  APJ loses R27 + Y33 (5.6 REU), uOR loses E54 + L32 (3.0 REU). APJ is the worse case. Expect the
  macrocycles to be substantially weaker binders than the parent nanobodies regardless of how
  faithfully the motif is transplanted.
- **Every built uOR design complex clashes with the receptor; the input pose does not.** All 68 of 68
  uOR complexes in `uor/interface/complexes/` have an interchain heavy-atom pair under 2.0 A (mean
  closest approach 1.19 A, worst 0.77 A) against 26 of 55 and 1.80 A for apj, while native
  `uor_complex.pdb` is clean at 2.27 A. The transplant, not the model, puts the motif Tyr 1.08 A into
  LYS238 CD. **Answered 2026-10-30 by job 50216944: recoverable.** A CA-restrained FastRelax clears
  the clash on all 6 uOR test structures (0.77-1.45 A -> 2.45-2.88 A) and flips `dG_separated` from
  +367/+913 to -14/-33, while the 3 apj controls barely move on geometry (2.36-2.49 -> 2.51-2.57 A).
  Full table and the three caveats -- the relax arm is not neutral, the peptide moves up to 1.26 A
  CA RMSD so motif fidelity needs re-measuring, and `dG_cross` is void in the relax arm -- are under
  item 2 of "What to do next". Results in `scripts/relax_iface/`, logs
  `logs/relax_iface_50216944_*.out`, protocol `scripts/relax_interface.sbatch`.

  **What remains open is not whether it is recoverable but what the relaxed poses are worth**: the
  clash relief comes from moving the design, and no post-relax motif RMSD has been computed. A pose
  that scores -33 by abandoning the transplanted motif is not a win.

  **ROOT CAUSE FOUND 2026-10-30: residues 102-107 of the uOR nanobody are D-amino acids.** Six
  consecutive residues, exactly the CDR3 tip, exactly the transplanted motif. Verified by the
  improper dihedral N-CA-C-CB: 102-107 run -50 to -92 deg against +57.6 +- 1.3 for every other
  residue in the chain, with controls at 0 inverted for the uOR receptor (277), the apj nanobody
  (107) and experimental 9v32 (545, which fixes the sign convention). Bond lengths are normal, so
  it is a mirror-imaged backbone -- a local Coot artifact; the source filename contains "flip"
  twice and `prep/nanobody_gpcr_prep.md` never checked chirality. All three uOR scaffolds inherit
  all six; apj scaffolds have none.

  Rosetta builds the designs correctly as **L** and CA-superimposes them onto **D** targets, so the
  backbone matches (dCA 0.14-0.57 A) while CB lands on the opposite side (**dCB 2.02-2.68 A**),
  driving the Tyr to 1.01 A of LYS238:NZ. The clincher is residue 106 -- the one motif residue the
  build left as D -- whose CB moves **0.24 A** where the five rebuilt as L move 2.0-2.7. Nothing but
  a chirality mismatch produces that pattern.

  **This supersedes two earlier diagnoses on this page, both now struck:** the "global strain"
  reading of `rama_prepro` 33-37 (it is D residues scored against L tables) and the "sidechain
  rotamer / never optimised against the receptor" reading (the designs are right; the target is
  wrong). No packer setting and no relax can fix mirror-imaged target coordinates. Full evidence in
  `cyclize/README.md` under "Actual root cause: uOR residues 102-107 are D-amino acids".

- **`validate`'s 2.2 A clash filter cannot see motif clashes, by construction.**
  `min_receptor_distance` is called with `skip_nums={motif residues}` (`cyclize/validate.py:96`), so
  the filter only ever measures residues the design *added*. The skip is deliberate and its stated
  rationale is sound -- without it, native contacts like Y102 are flagged as clashes and all 50
  structures fail -- but it assumes the motif keeps its native sidechain geometry, which is exactly
  what breaks above. **That is why all 68 uOR complexes passed validation carrying a 1.08 A motif
  clash.** Compounding it, `superimpose_on_motif` fits on **CA only** (`validate.py:39`), so nothing
  constrains the sidechains of a transplant whose whole purpose is reproducing sidechain contacts.
  This filter is not catching what it looks like it is catching, and it will not catch the same
  problem in any future target either.

  ```bash
  cd /scratch/users/bvbell/xlink/cyclize   # count complexes with a sub-2.0 A interchain pair
  python3 - <<'PY'
  import glob; import numpy as np
  for t in ('uor','apj'):
      bad=tot=0
      for f in glob.glob(t+'/interface/complexes/*.pdb'):
          P=[];R=[]
          for l in open(f):
              if l.startswith(('ATOM','HETATM')) and l[76:78].strip()!='H':
                  (P if l[21]=='P' else R).append([float(l[30+8*i:38+8*i]) for i in range(3)])
          if not P or not R: continue
          d=np.linalg.norm(np.array(P)[:,None]-np.array(R)[None],axis=-1); tot+=1; bad+=d.min()<2.0
      print(t, bad, '/', tot)
  PY
  ```

- **`total_score` varies by ~5000 REU across the interface set, which undermines the ranking.** In
  `apj_i1/interface/interface.csv` the scored complexes run from -175 REU (`medium__amide__DTHR-GLY`)
  to +4904 (`short__amide__HIS-GLY`), and several designs in the dG top-6 sit above +3200 --
  including the headline `short__amide__DGLU-GLY` at +3364. This is exactly the failure mode
  `interface.py` was written to avoid ("must repack identically across all designs or it ranks clash
  severity rather than binding"), so the condition needs checking before the dG ordering is trusted
  as more than a coarse sort. A sanity filter on `total_score` before ranking is the obvious guard.

  **Checked 2026-10-30, and the worry is half wrong: `corr(dG_separated, total_score) = -0.04`
  across all 142 apj structures.** The strain is *not* driving the dG ordering, so the ranking is
  not corrupted in the way this entry assumed and a `total_score` filter is not the right guard.
  But a 16-residue peptide carrying +3300 REU is still not a physical pose, and **the strain is
  unsourced**: it is not interchain clash (the top designs sit at a 2.44 A gap) and it is not
  scaffold-specific (median `total_score` by scaffold/route runs 997-1533 across all six
  combinations while individual structures span -175 to +5375). The one systematic signal is that
  `ssapp` is consistently lower-strain than `amide` at every scaffold length. **This needs its own
  diagnosis before any amide design is preferred over a disulfide one.**

  The guard that *does* work is filtering on geometry rather than energy -- see the apj triage
  entry below.

- **apj triage (2026-10-30): filter on geometry, and the winner changes.** apj is free of the uOR
  chirality problem (0 inverted in the complex, in all three scaffolds, and in the built designs),
  so these numbers measure something real -- but half the set is junk. Of 142 scored structures
  (excluding the known-bad `long__amide__DALA-GLY__c01`), only **65 score `dG_separated` < 0** and
  **77 carry a sub-2.0 A interchain pair**. Binned by closest approach:

  | closest approach | n  | mean dG | dG < 0    |
  | ---------------- | -- | ------- | --------- |
  | 0 - 1.5 A        | 54 | +72.2   | 21/54     |
  | **1.5 - 2.0 A**  | 23 | +135.2  | **0/23**  |
  | 2.0 - 2.5 A      | 65 | +14.0   | 44/65     |

  **Nothing in the 1.5-2.0 A band scores as a binder -- not one of 23.**

  The raw dG top-6 is one backbone wearing six side chains: `short__amide__` DGLU/DASN/DPRO/DASP/
  DGLN/HIS-GLY at dG -27.0 to -24.0, all at **identical gap 2.44 A and motif RMSD 0.64-0.66**, each
  carrying ~3300 REU of unsourced strain. That is exactly what "i+1 is conformationally inert"
  predicts; read it as a side-chain ranking on one pose, not six results.

  Filtering on **gap >= 2.0 A and motif_rmsd <= 0.3** leaves 36 structures and removes both failure
  modes at once. Within that set the ranking is:

  | design                            | dG    | total_score | motif RMSD | gap    |
  | --------------------------------- | ----- | ----------- | ---------- | ------ |
  | apj `short__ssapp__DCYD-DCYD`     | -21.1 | 780         | **0.21**   | 2.49 A |
  | apj `short__ssapp__CYD-CYD`       | -19.1 | 1098        | **0.11**   | 2.36 A |
  | apj_i1 `medium__amide__DTYR-GLY`  | -18.6 | 140         | **0.08**   | 2.33 A |
  | apj `short__amide__GLY-GLY`       | -17.7 | 197         | 0.10       | 2.42 A |
  | apj_i1 `medium__amide__DGLN-GLY`  | -15.0 | **-152**    | 0.05       | 2.42 A |

  ~~`short__ssapp__DCYD-DCYD` is the most defensible design in the project.~~ **Withdrawn hours
  later the same day — see the next entry. No design in the apj set is defensible.**

- **The apj macrocycles self-intersect — ~56% of them.** Found and verified 2026-10-30; this
  supersedes the triage entry above and withdraws its recommendation. *(Originally written as
  "all 143 of them" on a >=3.4 A cutoff that ignored atom type. Corrected: 80 of 143 genuinely
  clash below 2.5 A; 56 more have their closest contact at 2.5-3.0 A between N/O, which is a
  hydrogen bond. The closure defects below are what make it 0 of 311, not this.)*

  The ring *closes* correctly -- every peptide bond in `short__amide__DPRO-GLY` is 1.31-1.34 A --
  but the chain passes through itself:

  ```
  CA8  - CA10         1.87 A     (two non-adjacent alpha carbons)
  CA7  - CA10         2.40 A
  CA1  - CA9          2.84 A
  7:C  - 10:N         0.64 A     worst non-adjacent heavy-atom pair
  ```

  **GenKIC found a valid closure that is not a valid structure.** Across all 143 apj structures the
  worst non-adjacent intra-peptide contact has median **2.43 A**, p25 **1.39 A**, minimum
  **0.38 A**, and the largest is **2.99 A**. Read that last figure with the atom-type caveat
  above: a 2.5-3.0 A N...O pair is a hydrogen bond, so the honest count of genuine
  self-intersection is **80 of 143** below 2.5 A, not all of them.

  **It is in the builds, not the scoring.** The original
  `apj/output/short__amide__DPRO-GLY/short_0001.pdb` carries the identical defect (0.64 A, CA8-CA10
  1.87 A) and its build-time `score.sc` reads `total_score` 3358.7 with **`fa_rep` 3278.8**.
  `FastRelax` never relieved it because the motif is pinned and the ring cannot open.

  **This also sources the ~3300 REU strain flagged as unexplained in the entry above: it is
  intra-peptide `fa_rep`.** That question is now closed.

  **Why it inverted the ranking.** `dG_separated` is complex minus separated chains, so
  intra-peptide energy appears in both terms and **cancels exactly**. A self-intersecting peptide
  therefore pays *no* dG penalty, while its collapsed mass buries more surface and gains `fa_atr`.
  The metric actively rewards the defect:

  | worst intra contact | n  | median `total_score` | median dG |
  | ------------------- | -- | -------------------- | --------- |
  | < 1.5 A (severe)    | 39 | 2016                 | **-8.3**  |
  | 1.5 - 2.5 A         | 40 | 1191                 | +9.1      |
  | 2.5 - 3.0 A         | 64 | 1017                 | +9.4      |

  `corr(intra, total_score) = -0.61`. Every one of the top 7 by dG is among the worst offenders --
  `long__amide__DALA-GLY` at 0.38 A, then the six `short__amide__*-GLY` at 0.60-0.69 A.

  **The disulfides are better but not clean, and the earlier recommendation is withdrawn.**
  `short__ssapp__DCYD-DCYD` runs 1.72-1.93 A across its three copies. One check worth recording,
  since a naive scan would get this wrong: for `ssapp` the Cys pair is non-adjacent in sequence, so
  the S-S bond itself can look like a clash -- it does not here. SG-SG is a legitimate 2.37-2.52 A,
  and the offending contacts are backbone overlaps (`2:O - 10:CB` 1.93 A, `2:O - 10:C` 1.72 A,
  `8:O - 10:O` 1.98 A). Severely strained rather than topologically threaded like the amide
  winners, but not acceptable either.

  **Root cause: `validate.py` never checks the peptide against itself.** It gates on motif RMSD,
  seed phi, and design-versus-*receptor* distance (`min_receptor_distance`). Self-intersecting
  closures pass validation, enter the representatives, and win the ranking. Together with the
  motif-skipping clash filter noted above, that is **two** blind spots in the same function, and
  neither failure was visible to any existing filter.

- **One interface row is certainly junk: `apj_i1 long__amide__DALA-GLY__c01`.** It reports
  `dG_separated` **-100.4** with `ddg` **+69.8** -- inconsistent signs -- at `total_score` +4489 and
  `dG_dSASA_norm` -5.2, against -1.8 or shallower for every other design. Its sibling `__c00` scores
  a normal -2.8. It is the only design in the set below -30 REU and should be excluded, not treated
  as the best binder in the project.

- **Every `ssapp` run carries a ~10 REU phantom energy on residue 1.** `build.xml` strips
  `LOWER_TERMINUS_VARIANT` from residue 1 for all routes, but `ssapp` closes through S-S and never
  bonds anything to that N, so Rosetta scores its rama at phi = 0 -- a forbidden region. The artifact
  was present during `FastRelax` in all 600 apj + 600 uOR `ssapp` structures. **Fixed 2026-10-29** --
  the two movers are now guarded on the closure chemistry (`build.py`, emitted only when the ring
  closes through a peptide bond), and all 25 disulfide variants were **resubmitted** the same day
  (arrays 50171634 apj / 50171636 uor). The pre-fix runs are kept in `<target>/output_prebug/` so
  the effect can be measured rather than assumed: `scripts/compare_prebug.py` compares psi(1), the
  S-S distance and the yields across the fix. `scripts/rescore.py` reads PDBs as-is and never
  reproduced the artifact, so the rescored numbers were already clean.

  **The reruns barely moved the yields**: apj `ssapp` 232 -> **240** passed of 600 closed, uOR
  `ssapp` 197 -> **194**, uOR `ssmut` 50 -> **50**. So the validation filters (motif RMSD, seed phi,
  clash) are insensitive to the artifact. That does not settle the geometry question -- those
  filters are coarse, and a turn can shift well within them. Seeding makes this **partly** paired and
  partly not (`run_task` seeds from position in the task file): the eight uOR `long`/`medium` ssapp
  variants happen to reuse their seeds and so are an exact paired comparison, while the five
  remaining uOR variants and all twelve apj variants were reseeded and are two independent samples,
  where a +-8/600 difference is within sampling noise. Which is which, and why, is the seed-collision
  entry below. Distribution-level comparison via `compare_prebug.py` is valid throughout.

  **Measured 2026-10-29 (job 50173650, `compare_prebug.py`): the bug moved scores, not geometry.**
  psi(1) shifts by **+3.18 deg** on apj (-8.64 -> -5.47, against a within-variant spread of ~100 deg)
  and **-3.48 deg** on uOR; the S-S distance shifts by **+0.01 A** on apj and **-0.00 A** on uOR.
  Both are far inside the sampling noise, so the earlier inference that relax had optimised psi(1)
  against a bogus phi is **not supported**. The rescored numbers were already clean and the
  structures are usable as built.

- **The pre-versus-post disulfide comparison is seed-confounded, and it is confounded the *opposite*
  way from what was first recorded.** Explained in full 2026-10-29. Eight uOR variants (all four
  `long__ssapp__*`, all four `medium__ssapp__*`) are **50/50 byte-identical** to their pre-fix
  structures; `short__ssmut__CYS-CYS` is 28/50; the four `short__ssapp__*` and all twelve apj
  variants are fully new.

  They are **not** copies -- the files have distinct inodes, `links=1`, and fresh 13:23-13:26 mtimes,
  so Rosetta re-ran and wrote them. The cause is the seeding rule in `cyclize/slurm.py:143`
  (`run_task`): the seed is the variant's **position in the task file**
  (`seed = first + i + 1` -> `-constant_seed -jran $((1000 + seed))` in `runner.py`), deliberately so
  that each variant keeps a distinct RNG stream however the file is chunked. The consequence not
  anticipated is that **the same position in a later submission reuses the same seed.**

  - uOR's retry file `tasks/20261029-113412.txt` and its fix-rerun file `tasks/20261029-131905.txt`
    open with the identical eight entries (`long__ssapp__*` x4 then `medium__ssapp__*` x4), so those
    eight drew identical seeds. From position 9 the files diverge (`short__ssapp__CYD-DCYD` vs
    `short__ssapp__CYD-CYD`), so every `short` variant was reseeded.
  - apj never collided: its retry file `tasks/20261029-113409.txt` **omits `long__ssapp__CYD-CYD`**,
    which shifts every later variant by one position against `tasks/20261029-131902.txt`.

  **So the eight identical variants are the only clean test of the fix, not the broken ones.** With
  the seed held constant and only the protocol differing, output was bit-identical --
  `ModifyVariantType` is deterministic and consumes no RNG, so the phantom rama on residue 1
  demonstrably never changed the sampling. That independently confirms `compare_prebug`'s "scores,
  not geometry" result, on a controlled comparison rather than two independent samples.

  Conversely, **apj's pre/post comparison is the confounded one**: every apj variant was reseeded, so
  its differences are reseeding noise and carry no information about the fix. The earlier note that
  "apj's comparison stands" had it backwards. `short__ssmut__CYS-CYS` at 28/50 was also reseeded, but
  `ssmut` is near-deterministic (psi(1) 136.79 +/- 0.00), so most structures converge regardless.

  Practical consequence: **nothing here invalidates the structures or the yields.** The only wrong
  claim was calling "uOR ssapp 197 -> 194" an independent post-fix validation; it is 8 bit-identical
  variants plus 5 reseeded ones. If a genuine paired comparison is ever wanted, pin the seed
  explicitly instead of relying on task-file position.

  ```bash
  cd /scratch/users/bvbell/xlink/cyclize   # the collision, in one line
  diff <(tr ' ' '\n' < uor/tasks/20261029-113412.txt | head -8) \
       <(tr ' ' '\n' < uor/tasks/20261029-131905.txt | head -8) && echo "same 8 positions -> same seeds"
  ```
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
cd ../flexddg && ./bin/flexddg status cldn6_9v32fab   # all five flexddg scans are done
```

**Picking up where 2026-10-30 left off**, in priority order. The 10-29 jobs are read out; the
rebuild chain (50318433 -> 50319161 -> 50321818) is in flight and is item 1.

1. **Read out the rebuild chain** — relax (50318433) -> re-excise (50319161) -> geometry test
   (50321818). The decisive number is how many structures in each arm and route are valid
   (route-aware: amide closure `O...N` 2.15-2.35 A, or disulfide S-S 1.95-2.15 A with |chi3|
   60-120 deg; both plus zero non-proline cis and every bond within 30 deg of planar). Expect
   `base` ~0,
   `omega` to be the fix; `loop7` only matters if `omega` cannot close. Do **not** start a
   production run before this reads out. See "The 2026-10-30 rebuild".
2. **Re-run the `flexddg/uor_nb/` alanine scan on the relaxed uOR model.** This is the uOR item
   that survives the design-basis decision: the scan that *selected* the uOR motif ran on the
   D-residue model, so its 8 hotspots and 13.6 REU are measured on inverted stereocentres and
   Y102's hotspot status is not trustworthy. Designing against relaxed structures does not fix
   that, and the corrected scan may change which residues the motif should contain. Do this before
   any uOR result is interpreted as hotspot-driven.

   *(The former item 1, "repair uOR 102-107 to L chirality", is superseded by the 2026-10-30
   decision to design against the relaxed structures — the relax performs the L conversion and
   makes the motif source and alignment target consistent, which is what the design needs.
   Settling the hand against the J126 density is still worthwhile for the structural record but is
   off the critical path.)*
3. **Fix the two blind spots in `cyclize/validate.py`, then re-audit everything already built.**
   Both failures above slipped through the same function:
   a. an **intra-peptide clash gate** -- it currently checks the design against the receptor but
      never against itself. Add it, re-run over the existing output, and size the damage properly
      (the 143-structure sample says 0 survive at 3.0 A, but that is apj's interface set, not the
      full 2750);
   b. a **chirality gate** (improper dihedral N-CA-C-CB, three lines), and the same check in the
      prep pipeline -- `prep/nanobody_gpcr_prep.md` verified occupancies, altlocs, missing atoms
      and bond lengths but never chirality;
   c. note the existing receptor-clash filter skips motif residues by construction, which is the
      third blind spot (see Open questions).
   This is the highest-leverage item on the list: without it, a rebuild reproduces the same defects.
4. **Work out why GenKIC cannot close these rings physically, before rebuilding apj.** A 2-residue
   turn may simply be too short to span the hairpin gap, which would tie directly to the `i+j=209`
   register constraint that forces even-length motifs and 0- or 2-residue closures. If the geometry
   is infeasible, no amount of sampling fixes it and the register itself has to be revisited. This
   is a design question, not a scripting one, and it gates the apj rebuild.
5. **Only then rescore apj.** The scoring fixes below are real but pointless against bad builds:
   a. test `-ex1 -ex2` (the emitted XML uses the default rotamer library with no TaskOperations),
      then a motif-sidechain-restrained relax if needed, keeping the route-aware
      `PeptideCyclizeMover` guard;
   b. a `--native` flag that scores the native peptide through the identical pipeline and reports
      **design - native**. `dG_separated` is explicitly not a ddG -- and, as the self-intersection
      finding shows, it is blind to intra-peptide energy entirely, so it must never be the only
      ranking criterion again.
6. **`cldn6_9v32fab` needs a go/no-go decision, and the evidence now points at no-go.** Resolved
   2026-10-30: the interface is a one-residue pharmacophore because the epitope is only 11% conserved
   (see the scan section above), and the dead heavy chain is an artifact of the chimera holding the
   Fab in its CLDN18.2 pose. Nothing further to check -- the question is whether to spend scaffolding
   effort on one tyrosine, and the honest answer is that APJ and uOR are far better starting points.
7. **Make `run_task` seed explicitly instead of by task-file position** (`cyclize/slurm.py:143`).
   The 8 bit-identical uOR variants are explained -- a seed collision between two submissions whose
   task files share their first eight entries, see Open questions -- so this is no longer a mystery,
   just a latent trap: any future resubmission that happens to preserve ordering silently reruns the
   identical trajectory and looks like a confirmation. Recording the seed per variant in
   `tasks/submissions.tsv` would make pre/post comparisons paired on purpose rather than by accident.
8. **~~Replace the uOR input model.~~ Closed 2026-10-30 by the decision to design against the
   relaxed structures.** The relaxed model *is* the working input now. For the record, the two
   earlier reasons for replacing it have both dissolved: 50216944 showed the transplant clash, not
   the model, accounted for ~330 of the ~400 REU; and the `rama_prepro` strain of 33-37 was the
   D residues scored against L tables, which the relax resolves (31.96 -> 9.50 on chain A). What
   remains is not a modelling task but a measurement one — re-run the `uor_nb` scan (item 2) —
   plus the optional density work noted there.
9. `design_i2`, and the APJ CDR1/CDR3 surface question, both still open and both unstarted.

Do **not** re-run the i+1 sweep or add nstruct to it: item 1 under "What to do next" shows the
sampling is saturated and i+1 is conformationally inert, so more of either buys nothing.

## Repo hygiene

Everything is in `$SCRATCH`, purged at 90 days of no modification, so the git remote
(`git@github.com:bvbell/Xlink.git`) is the durable copy.

Committed: the `cyclize` package and both target configs, the scaffolds, both complex PDBs, and the
`apj_nb`/`uor_nb` alascan summaries and plots. The ~6 GB of trajectories under `*/output/` is
gitignored and regenerable.

All of that is now committed (2026-10-29): `prep/`, `scripts/`, `flexddg/cldn6_9v32fab/`, the chimera
and relaxed-chimera PDBs, and the raw Coot `Rosetta_9_real_trunc_*` and `d006_lmpnn_fkd2_*.cif`
models, which `prep/nanobody_gpcr_prep.md` and `scripts/prep_structures.py` reference as the input
provenance for `uor_complex.pdb` and `apj_complex.pdb`. The 17 MB of relax decoys under
`flexddg/cldn6_9v32fab/relax/decoys/` stays gitignored and is regenerable from `relax.sbatch`.

**Pushed through `7955955` (verified against the remote with `git ls-remote`).** The claim that
nothing had ever been pushed was wrong. What is *not* on the remote is **all of 2026-10-30**: the
invalidation findings, `cyclize/ARCHIVE.md`, the `omega_trans` fix, the rebuild scripts and the
inspection set — 22 uncommitted paths. The remote is the durable copy and `$SCRATCH` purges at
90 days, so committing today's work is the outstanding housekeeping step.

Uncommitted as of 2026-10-30: this file, `flexddg/README.md`, and the new
`flexddg/cldn6_9v32fab/analysis_output/` + `plots/` from the finished scan. The latter two match
what is tracked for `cldn18_9v32/` and should be committed the same way. Job 50216944's outputs in
`cyclize/scripts/relax_iface/` are untracked and not covered by `.gitignore`; they are cheap to
regenerate (~3 min for all 9) but the `score.sc` files are the only record of the as-is/relax
comparison, so commit those if the protocol changes under them.

`motif_scaffold/` checked 2026-10-29: the 512 rows of `scores_prod.csv` correspond 1:1 to 512 PDBs in
`motif_scaffold/runs_prod/` (plus 1024 trajectory frames under `runs_prod/traj/`), so that scored set
is intact and nothing was lost. `outputs/` holds only a hydra log directory, which is why it looks
empty -- the designs are in `runs_prod/`, not `outputs/`. Those files were last written 2026-10-22/23,
which puts them in the purge window around **2027-01-20**; they are untracked and not regenerable
without re-running RFdiffusion, so they are the clearest candidate for a copy to Oak.

Copy final designs and summaries to `$GROUP_HOME` or Oak.

Rosetta: use `ml rosetta/3.14` for anything you run; `/home/groups/possu/Rosetta` is an older checkout
kept for reading source, offline docs and integration tests.
