# Cyclising CDR3 hairpin motifs into macrocyclic peptides

`cyclize/` is the package; each target is a directory with a `cyclize.json`. It takes the CDR3 motifs
selected from the flex ddG alanine scans (`../flexddg/`) and closes them into cyclic peptides by amide
or disulfide, holding the ddG-validated binding conformation rigid.

| dir     | complex                | scaffolds                   | receptor | scan source          |
| ------- | ---------------------- | --------------------------- | -------- | -------------------- |
| `apj/`  | `../apj_complex.pdb`   | `../apj_scaffolds/*.pdb`    | chain B  | `../flexddg/apj_nb/` |
| `uor/`  | `../uor_complex.pdb`   | `../uOR_scaffolds/*.pdb`    | chain B  | `../flexddg/uor_nb/` |
| `apj_i1/`, `uor_i1/` | as above | as above | chain B | as above |

`apj_i1/` and `uor_i1/` are the i+1 turn sweep (`redesign_i1`, amide route only); same inputs as
`apj/` and `uor/`, separate directories so `setup` cannot overwrite the finished grid.

Everything lives in `$SCRATCH`, which purges files not modified for 90 days: push this repo, and copy
final designs and summaries to `$GROUP_HOME` or Oak.

## The register is the whole design constraint

Both CDR3s are antiparallel beta-hairpins, and every backbone H-bond in both obeys **i + j = 209**:

```
APJ:  95<->114  97<->112  99<->110  101<->108  103<->106
uOR:  95<->114  97<->112  99<->110  101<->108
```

Three consequences, and they drive every decision in this package:

1. **A motif is usable only if `first + last == 209`.** All six scaffolds satisfy it. Because i and
   209 - i are always an odd number of positions apart, an on-register motif always has an **even**
   residue count -- parity is the symptom, the register is the cause. `cyclize register <target>`
   checks this; `setup` refuses an off-register scaffold outright.
2. **Closure must add 0 or 2 residues, never 1.** The two termini of an on-register motif are a
   reciprocal H-bonded pair (101 N...O 108 at 2.81 A and back), which is why they sit only 4.6-6.3 A
   apart. An odd insertion shifts the ladder by one and unzips it. This is why asymmetric disulfides
   (append one end, mutate the other, = +1) are deliberately not enumerated.
3. **Direct head-to-tail amide is geometrically dead.** A trans peptide bond needs Ca-Ca ~3.8 A; the
   termini are 4.6-6.3 A apart. Amide closure starts at 2 added residues.

## Routes

| route   | added | ring       | what it does                                                     |
| ------- | ----- | ---------- | ---------------------------------------------------------------- |
| `amide` | +2    | motif + 2  | a seeded beta-turn bridging the terminal pair                    |
| `ssapp` | +2    | motif + 2  | one Cys beyond each terminus (prepend + append), S-S between them |
| `ssmut` | 0     | motif      | both terminal residues mutated to Cys, stapling the native pair   |

`ssmut` is only offered where the terminal pair's **Cb-Cb already sits in the natural disulfide window
(3.4-4.6 A)**. Cb position is fixed by the backbone, so a pair outside that window can never be
stapled by mutation, no matter how it repacks. Measured: only **uOR/short (4.4 A)** qualifies; the
other five are 5.6-6.0 A and need Cys appended on free backbone instead.

## Turn seeding

A hairpin-capping turn must be type I' or II' (the mirror types; I and II impose the wrong sheet
twist), and both need **positive phi at i+1**. Residue identity is what makes that reachable, so
seeding the identity is seeding the geometry:

| residue | accessible phi      | why                                                       |
| ------- | ------------------- | --------------------------------------------------------- |
| Gly     | either sign         | no Cbeta                                                   |
| D-Pro   | **+63 +/- 15 only** | pyrrolidine ring                                           |
| L-Pro   | **-63 +/- 15 only** | ring, mirrored -- so **never legal at i+1**                |
| Asn     | - preferred, + reachable (~+55) | side-chain-to-backbone H-bond stabilises alphaL |

`variants.check_seed` enforces this, so a seed pair that cannot make the turn type it claims is a
hard error rather than a wasted run. The five seeds tried are in `variants.TURN_SEEDS`.

### Sweeping i+1 (`redesign_i1`)

The five hand-picked seeds leave almost nothing to chemistry: four of the five put **Gly at i+2**, and
four of the five put Gly or D-Pro at **i+1**. The one position the closure chemistry allows us to
redesign without disturbing anything else is i+1 of the amide route -- i+2 stays Gly, the disulfide
routes add only cysteines, and the pinned motif is never edited. `variants.sweep_seeds()` enumerates
it, gated by `FIXED_IDENTITIES` (Gly, Pro, D-Pro and the cysteines are there for geometry or for the
S-S, so redesigning them would destroy what makes the variant work).

**This is enumeration, not Rosetta design, and that is deliberate.** Identity at i+1 *is* the
geometry: GenKIC's `randomize_backbone_by_rama_prepro` perturber draws that residue's torsions from
its own Ramachandran table, so the backbone must be re-sampled per identity regardless. Mutating i+1
on an already-closed backbone would force the new residue into the torsions chosen for the old one,
and any `rama_prepro` measured afterwards would be scoring it against the wrong backbone. Since the
space is a single position, enumeration is exhaustive where packing would only sample it -- and it
needs no PackRotamers, no resfile and no D-rotamer handling anywhere in the protocol.

`PHI_POS` now admits **every** D residue. `rama_prepro` has no D tables of its own: it mirrors the L
table (`core/scoring/RamaPrePro.cc:100`, *"Score tables are loaded for L-amino acids; D-versions are
mirrored"*), so a D residue's favoured basin is the mirror of its L counterpart's, which is the
positive-phi region. D-Pro remains the special case -- its ring narrows that to about +63 rather than
merely allowing it. D residues are spelled `D` + name3, generated by the chiral-flip patch
(`database/.../patches/D_AA.txt`, `NAME D`, selector `PROPERTY L_AA`); there is no D-Gly.

That gives **21 legal identities** at i+1: ASN ASP HIS SER and 17 D residues. Probed end to end
against rosetta/3.14 with `cyclize smoke` -- `DTRP`, `DARG`, `DGLU` and `DILE` each closed 2/2, placed
the right residue (`DTR`, `DAR`), and came out with positive phi at i+1 (+101.5, +140.9). All 19 D
names were separately confirmed valid for `PeptideStubMover`. Two are **not**, for the same reason
each time -- the patch selects on `PROPERTY L_AA` and names the result `"D" + name`, so anything that
is not an L_AA base type has no D form: **`DGLY` is invalid** (glycine is achiral -- use `GLY`), and
`DCYD` is invalid (`CYD` is a variant, not a base type -- use `DCYS:disulfide`). Resname validity is
worth probing for anything new here; a bad one does not error, it silently closes 0/50.

**Result (2026-10-29): i+1 identity does not change the backbone.** The sweep ran to completion --
126 variants, 2483 passing structures -- and then clustered to almost nothing. On apj/short, 20 i+1
identities and 850 passing structures give **2 backbones** at 1.0 A, 656 of them in a single cluster;
apj/long gives 3 over 340 structures, uOR/medium gives **1** over 299, uOR/long gives 2 over 232.
Only uOR/short keeps any spread (8 clusters over 549).

So the sweep answered its question in the negative: i+1 is **conformationally inert**, and a pass
count per identity carries no structural information. The earlier reading ("D-Ala is not special --
six identities all close 50/50") is true but trivial, because those six are one backbone wearing six
side chains. What this does buy is the thing the sweep was for: i+1 can now be chosen on chemistry
and on interface score, with the confidence that doing so will not cost closure. Caveat: the sweep is
amide-only, so some of the collapse against the 55-variant grid's 5-6 clusters per scaffold is route
homogeneity -- `ssapp` supplied most of the grid's conformational diversity.

Reproduce with `python3 scripts/degeneracy.py apj uor apj_i1 uor_i1 --csv-dir=scripts`; the
cross-variant block at the end of each target is the relevant one, and the CSVs it writes now cover
all four targets.

## What the sampling actually does, and what it cannot

Measured during development on apj/long; these are the non-obvious results worth keeping:

- **GenKIC's 3-residue loop minimum is a hard floor, not a working configuration.** With CA pivots on
  three consecutive residues all six pivot torsions are solved by the closure algebra and only the
  omegas stay free, so the system is effectively determined: 2000 attempts gave zero solutions.
  `variants.TARGET_LOOP` = 5 releases enough cold motif residues to leave real sampling freedom.
  Scaffolds that cannot reach it are flagged `underdetermined` in `setup.log`.
- **`rama_prepro_check` is the dominant throttle.** 4 structures x 5000 attempts: filter off 4/4
  closed, cutoff 5.0 gave 1/4, cutoff 1.0 on the seeded residue gave 0/4. It is therefore **off by
  default** (`protocol.rama_cutoff: null`); backbone quality is judged after the fact instead.
- **A seeded residue cannot be a pivot and keep its chirality.** D-Pro as pivot 1 came out at phi -55
  and +145, which its ring forbids. But moving the pivot off it starves the closure (0/4 at any
  cutoff). So pivots stay spread and `validate` checks the seed's phi sign afterwards -- expect to
  discard roughly 30% of closed designs on handedness.
- **The motif must be pinned explicitly.** An unrestrained `FastRelax` moves it ~5 A and discards the
  entire point of grafting a validated conformation. Relax is confined to the GenKIC loop by MoveMap
  and the motif is pinned by CA coordinate constraints; measured motif RMSD is then 0.56-0.99 A.
- **Rosetta does not preserve the input coordinate frame.** RMSD must be measured after superposition;
  comparing raw coordinates reports several angstroms of pure rigid-body displacement.
- **Rosetta writes D-amino acids as `HETATM`**, with name3 `DPR`, `DAL`, `DCS` etc. A reader that only
  takes `ATOM` silently drops every designed turn residue.
- **There is no `DCYD`.** The disulfide-bonded cysteine is `CYD` for L, but the D form must be spelled
  **`DCYS:disulfide`** -- `DCYD` fails with "The residue DCYD could not be generated", which is why 9
  of 12 `ssapp` variants per target closed 0/50 on the first submission while all three `CYD-CYD` ones
  reached 50/50. Probed against rosetta/3.14: `CYD` ok, `DCYS` ok, `DCYS:disulfide` ok, `DCYD` invalid,
  `CYD:D` invalid. `variants.CYS_FORMS` therefore carries `(label, resname)` pairs, since the colon
  cannot appear in a variant id.
- **A disulfide closure needs explicit `randomize_dihedral` perturbers.** `randomize_backbone_by_rama_prepro`
  aborts on a residue missing a polymer connection, and in the ssapp topology the ring closes through
  the S-S, so both Cys are polymer-terminal. Simply excluding them from that perturber is not enough:
  nothing then moves them off the extended geometry PeptideStubMover built, and the sulfurs never meet
  -- ssapp closed **0/50** that way, with no error, 100,000 attempts. `build.ss_torsions` mirrors
  simple_cycpep_predict's terminal_disulfide setup (the S-S torsion, both side chains CA-CB-SG, and
  what backbone each terminal Cys has), which takes it to 8/8.
- **Sampling saturates long before nstruct 50, so pass counts badly overstate how many designs there
  are.** Clustering the passing structures on N/CA/C RMSD (`scripts/degeneracy.py`) gave **120
  distinct conformations at 1.0 A** -- apj 52, uOR 68. (The earlier 112 / 51 / 61 were measured on
  the pre-rerun set of 1200; see ../STATUS.md.) The variants with the best pass rates are the
  most collapsed, not the most productive: apj `short__amide__DALA-GLY` and `DPRO-PRO` are 50/50 and
  still a single cluster at **0.25 A**. Pooling across turn seeds is worse -- apj's 402 short
  structures from nine seed/route combinations are **5** backbones, uOR's 122 medium are **2**. The
  pinned motif, not the seed, is dictating the closure. Two consequences: more nstruct buys nothing,
  and anything downstream should run on `scripts/representatives.py` output rather than on all 50
  members of a cluster.

The motif is located **by position in the design, not by residue number**. Routes that add residues
make Rosetta renumber the pose 1..N, but `ssmut` adds nothing so the input PDB numbering survives
(101-108). A numeric range matches one and not the other, which made uOR `ssmut` report its own motif
hotspot Y102 as an "added residue" clashing with the receptor and fail all 50 structures. For `ssmut`
the reported receptor distance is `inf`, correctly: it adds no residues, so there is nothing to check.

## Clash threshold

`validate.clash` defaults to **2.2 A, not 3.0**. The native APJ interface packs motif residues to
2.44-2.80 A of the receptor (F104 2.44, E108 2.65, K110 2.70), so a 3.0 A cutoff rejects the native
structure itself. 2.2 A is below anything observed natively, so it flags only real overlap.

## Scoring a closed design

Nothing in the package reads energies yet -- `validate` filters on geometry only, and the `score.sc`
Rosetta drops in each `output/<vid>/` (21 columns, `rama_prepro` among them) is written and never
read. The mechanics below were verified against rosetta/3.14 before any of it gets wired in, because
the obvious way to do it is wrong.

### Reading a cyclic pose naively gives the wrong `rama_prepro`, silently

`rama_prepro` is a **two-body** term on consecutive pairs: the value "for residue i" lives on the pair
`(i, i+1)`, and the cyclic wrap term is the row `(1, N)`, belonging to residue N. Reading an
amide-closed design with `score_jd2`, `per_residue_energies`, `residue_energy_breakdown` or plain
`rosetta_scripts` reports a phi(1) taken from the atom-tree root stub rather than from the geometry
(measured: -46.2 deg against a true 176.9 deg), and so gets the seam residues wrong.

The decisive test is **cyclic-permutation invariance** -- renumber the same ring so residue 8 becomes
residue 1 and rescore. A correctly handled cyclic pose is bit-identical; the naive read is not, and it
moved the total from 31.78 to 36.63 purely by changing where the numbering starts. The error is
therefore silent *and* variable, not a constant offset you could ignore.

On `apj short__amide__ASN-GLY/short_0001.pdb` the naive total is 22.755 against a true 29.356, and the
understated residues are 10 (by 3.6 REU) and 1 (by 3.0 REU) -- which are, exactly, two thirds of that
variant's GenKIC loop `[9, 10, 1]`. The term is wrong precisely where the turn is.

Rosetta writes a `LINK` record for the cyclization bond and reads it back, but the LINK only
half-works: it suppresses the charged-terminus patches while leaving the torsion bookkeeping wrong.
Stripping the LINK is worse -- residue 1 picks up `NtermProteinFull`, `rama_prepro(1)` collapses to
0.000 and `fa_rep` jumps from 199 to 986. **No command-line flag fixes this.** The fix is one mover:

```xml
<MOVERS><PeptideCyclizeMover name="cyc" residue_selector="pep"/></MOVERS>
```

It replaces `ModifyVariantType x2 + DeclareBond` (verified identical), and works even on a PDB with
the LINK stripped. It also adds cyclization constraints to the pose -- harmless under `ref2015`, which
carries no constraint terms, but not under a constraint-weighted function.

Two further rules:

- **For a per-residue breakdown, use `residue_energy_breakdown`, not `PerResidueEnergyMetric`.** The
  metric splits every pair 50/50 between its two residues (verified: res1 = 5.487 = half of
  `(1,2)`=5.017 plus half of `(1,14)`=0.4695), which leaks half of each boundary residue's rama into
  its neighbour -- fatal when the whole point is to sum over a subset.
- **Pass the cyclized pose on as a binary silent file.** It round-trips the declared bond and variant
  types; a PDB does not, and writing one puts you back in the naive state.

### The disulfide routes take the opposite rule

For `ssapp`/`ssmut` the backbone genuinely is linear, closed through S-S (measured on
`long__ssapp__CYD-CYD`: d(C14,N1) = 4.09 A, d(SG1,SG14) = 2.42 A). Rosetta picks the disulfide up from
`SSBOND`/geometry with no flag, so **the naive read is the correct one here**.

**Never apply `PeptideCyclizeMover` to an ssapp or ssmut design.** It will force an N-to-C amide
across that 4.09 A gap and report 27.861 instead of 17.828, with no error raised.

Note also that residue 1 has no preceding residue and residue N no following one, so `rama_prepro`
does not exist for either appended Cys -- the two residues `ssapp` actually adds cannot be judged by
this term at all.

### `scripts/rescore.py`

All of the above is implemented. It splits the structures by route, runs each group under the right
protocol, hands the pose on as a binary silent, and sums the breakdown's pair rows over the variant's
`rosetta_loop` from `input.json`. Output is `<target>/rescore/rescore.csv`; the representatives
score in about 25 s per target.

Two details in the arithmetic worth knowing before trusting a column:

- **The wrap pair `(1, N)` belongs to residue N, not residue 1.** Attributing it to 1 moves the
  turn's largest contribution onto a pinned motif residue and quietly flatters the loop score.
- **`rama_loop` is blank, never 0.0, when nothing scored.** `ssmut` has an empty GenKIC loop, and a
  zero there sorts to the top of a ranking as if it were the best turn in the set. `n_loop_scored`
  shows how many of the loop's residues actually carried a value, which also makes visible that on a
  disulfide route the two chain-terminal residues have none.

Validated two ways. Rescoring `short__amide__ASN-GLY/short_0001.pdb` reproduces 29.356 with every
pair matching to the digit. And comparing each representative's run-time `score.sc` against its clean
rescore separates the routes exactly as the bug predicts:

| route | run-time minus clean rescore | n |
| ----- | ---------------------------- | - |
| `ssapp` | **+10.61 mean, +9.42 to +12.22** | 28 |
| `amide` | +4.06 mean, +0.84 to +7.69 | 23 |

The `ssapp` offset is the phantom: nearly constant, the predicted size, and absent from the clean
read. The `amide` offset is smaller and variable; it is most likely PDB coordinate rounding pushed
through a steep potential, but that is **not confirmed**. It does not affect ranking, since every
structure is rescored by the identical path -- but it does mean a rescored number is not comparable
to a number lifted from `score.sc`.

For the interface half, see `scripts/interface.py` below.

### Interface score

The outputs are **peptide-only, chain A**; the receptor is never present during a run. Any interface
score therefore needs a two-chain complex built first, by applying `validate.superimpose_on_motif`'s
transform and concatenating the receptor chain.

Use the RosettaScripts form (`PeptideCyclizeMover` + `InterfaceAnalyzerMover` + `Ddg` filter,
`pack_separated="true"`) rather than the standalone `InterfaceAnalyzer` app: ~2.4 s/structure against
~2.0 s, and the premium buys cyclization in the same process, so the peptide is scored as the
macrocycle it is rather than as a linear chain with charged termini sitting in the binding site. Feed
it a `-l` list so the ~3 s startup is paid once per batch. For comparison, the flex ddG protocol costs
**~19 CPU-hours per point** (1944 s x 35 trajectories, measured in `flexddg/apj_nb/`) -- a validation
tool for a shortlist, never a ranking tool.

**The caveat that decides whether any of this means anything:** on a raw superposed complex both
numbers are dominated by `fa_rep`. An unrepacked probe complex gave `dG_separated` +226 and `ddg`
+93, because even a 0.71 A motif superposition leaves the turn clashing into the receptor. Without a
repack before scoring, **you are ranking clash severity, not binding** -- and whatever that step is,
it must be identical across every design or the ranking is meaningless. `dG_separated` is a binding
score in REU with no reference-state correction, so it is only ever a relative ranking.

`scripts/interface.py` implements this: it builds the complex, then scores with
`pack_separated="true" pack_input="true"` so every design gets the same sidechain repack. That is
enough to put apj into sensible territory -- the 52 representatives run **-25.2 to positive**, with
`short__amide__DPRO-GLY` best at -25.2 REU over 1429 A^2 of buried surface. It is a **sidechain**
repack only; backbone clash from the superposition is not relieved, which is intended, since
`validate` has already rejected anything closer than 2.2 A. Cost is ~3.6 s/structure, so the 51
representatives take about three minutes.

> **That "which is intended" was wrong, and it is the bug this whole section failed to catch**
> (established 2026-10-30, job 50216944). `validate`'s 2.2 A filter is measured on the **motif**, not
> on every transplanted sidechain, so it lets through complexes with sub-1.0 A interchain pairs
> elsewhere -- and a sidechain repack cannot escape them. Adding a CA-restrained relax before scoring
> moves uOR from +367/+913 to -14/-33 and apj by a further 15-25 REU. **Every number in this section
> and in the committed `interface.csv` files is an as-is number and is therefore suspect.** See
> "Resolved 2026-10-30" below before using any of them.

**uOR interface numbers are not usable.** Its 68 representatives come out at **+370 to +413 REU**
against apj's -25. In the *input pose* this is not a single bad contact: measured in the current
`uor_complex.pdb`, Y102:CD1 to K238:CD is **3.39 A**, an ordinary van der Waals contact, and the
whole chain A/B interface has **nothing under 2.0 A**, closest approach 2.27 A. So the 1.08 A figure
under Known caveats is not a property of the input being used.

**But it is a property of the built complexes** (measured 2026-10-29). The transplant puts it back:
**68 of 68** complexes in `uor/interface/complexes/` have an interchain heavy-atom pair under 2.0 A,
mean closest approach **1.19 A**, worst **0.77 A**, and the recurring offender is the transplanted
motif Tyr at **1.08 A** from LYS238:CD. apj, for contrast, is 26 of 55 (47%) at a mean of 1.80 A --
which tracks its usable scores. Backbone-backbone pairs under 2.8 A are zero in every complex
checked, so these are sidechain overlaps that a *restrained* relax may be able to clear even though
the sidechain-only repack in this protocol cannot.

This matters because it is a candidate explanation for the residual the native control below leaves
open. `validate` rejects anything closer than 2.2 A, so these complexes passed the clash filter --
worth knowing that the filter is measured on the motif, not on every transplanted sidechain.

The real finding is worse, and comes from scoring the **native** motif:receptor pair -- the actual
nanobody CDR3 in its actual position, a real binder -- through the identical protocol:

| control | `dG_separated` | `ddg` | `total_score` |
| ------- | -------------- | ----- | ------------- |
| apj short / medium / long  | -5.7 / -12.9 / -19.1 | -4.8 / -15.5 / -13.3 | ~-240 |
| uOR short / medium / long  | **+34.4 / +41.7 / +45.2** | +97.3 / +91.6 / +93.8 | **+180 to +209** |

> *Superseded 2026-10-30 — kept for the reasoning trail. The control below was an as-is
> measurement dominated by the transplant clash, and the `rama_prepro` evidence behind it is the
> D-vs-L table mismatch. Do not quote the +34 to +45 floor. See "Resolved 2026-10-30" and "Actual
> root cause" below.*

**The native uOR complex scores as a non-binder.** apj's natives sit at -5.7 to -19.1 REU with a
total score near -240; uOR's are positive on every measure. That is a global property of the
Coot-built model, not a single repairable contact, and it sets a floor under every uOR number: no
design can score well against a receptor whose own native complex does not.

**The two findings compose; they do not compete.** The native control puts a floor at +34 to +45, and
the designs land ~330 REU below *that*, which this table explicitly did not explain. The universal
transplant clash above is the leading candidate for that residual: the native pair is clean at 2.27 A
and the designs average 1.19 A, which is exactly the kind of gap `fa_rep` turns into hundreds of REU.
If so, the uOR problem is **two** problems -- a strained model setting the floor, and a transplant
artifact adding the rest -- and only the first needs a new input structure.

### Resolved 2026-10-30: the clash half was right, the floor half was not

`scripts/relax_interface.sbatch` (job 50216944) scored 6 uOR and 3 apj complexes as-is and after a
CA-restrained FastRelax (sd 0.5 A) through the identical protocol; results in `scripts/relax_iface/`,
logs `logs/relax_iface_50216944_*.out`. The as-is arm is validated -- it reproduces
`long__amide__DPRO-GLY__c00` at `dG_separated` 367.635, matching the committed CSV exactly.

The relax clears the clash on every uOR structure and takes the score negative with it:

| structure                     | as-is `dG_separated` | relaxed | closest approach |
| ----------------------------- | -------------------- | ------- | ---------------- |
| uor `long__amide__GLY-GLY`    | +369.2 | **-33.0** | 0.77 -> 2.88 A |
| uor `long__ssapp__DCYD-DCYD`  | +401.1 | **-29.9** | 1.45 -> 2.85 A |
| uor `medium__amide__ASN-GLY`  | +526.0 | **-31.2** | 1.14 -> 2.45 A |
| uor `long__amide__ASN-GLY`    | +367.3 | **-24.1** | 1.29 -> 2.73 A |
| uor `short__ssmut__CYS-CYS`   | +913.0 | **-20.8** | 1.09 -> 2.64 A |
| uor `long__amide__DPRO-GLY`   | +367.6 | **-14.2** | 1.15 -> 2.76 A |
| apj `short__ssapp__DCYD-DCYD` | -19.4  | -44.5     | 2.49 -> 2.52 A |
| apj `short__amide__DPRO-GLY`  | -20.7  | -37.4     | 2.44 -> 2.57 A |
| apj `short__ssapp__CYD-CYD`   | -22.4  | -36.8     | 2.36 -> 2.51 A |

So the ~400 REU was **this protocol's refusal to relieve backbone**, not the input model: one
sidechain overlap no rotamer can escape dominated the entire number. The prediction two paragraphs
up -- "sidechain overlaps that a *restrained* relax may be able to clear" -- held exactly.

Three things bound that conclusion:

- **The relax arm is not a neutral control.** apj improves 15-25 REU too, so relaxed and as-is
  numbers are not comparable and the whole set must be rescored before anything is ranked. The
  relaxed ordering already differs from the as-is one: `short__ssapp__DCYD-DCYD` (-44.5) overtakes
  `short__amide__DPRO-GLY` (-37.4).
- **The relief comes from moving the design, and motif fidelity was not re-measured.** The receptor
  barely moves (CA RMSD 0.13-0.17 A) but the peptide does: 0.43-0.72 A on uOR, up to **1.26 A (max
  2.88 A)** on apj `short__amide__DPRO-GLY`. Against original motif RMSDs of 0.09-0.65 A that is a
  shift comparable to the whole motif-fidelity budget. **Post-relax motif RMSD against the parent
  nanobody is the missing measurement**; until it exists, a relaxed pose scoring -33 cannot be
  claimed to still present the pharmacophore it was built for. If it turns out large, the protocol
  needs motif-atom restraints rather than CA-only.
- `dG_cross` reads 0.000 in every relax arm (the `VirtualRoot` changes the fold tree). Use
  `dG_separated` and `ddg`, which are consistent across both arms.

**This also undercuts the native-pair floor in the table above.** Relaxed uOR designs (-14 to -33)
score *better* than the +34 to +45 that control reported -- so the control, also run as-is, was
itself measuring clash rather than strain. A residual uOR-vs-apj gap of roughly 10-25 REU survives
under relax and is plausibly the real model-strain penalty, but **the floor is now unquantified**
and has to be re-measured under the relax protocol before +34 / +41 / +45 is quoted again. The
~~independent strain evidence is untouched: uOR `rama_prepro` still runs 33-37 against apj's
10-14.~~ **Superseded 2026-10-30:** that gap is D residues scored against L Ramachandran tables,
not strain — a CA-restrained relax of chain A takes `rama_prepro` from 31.96 to 9.50.

Caveat on the native control itself: it is a real per-scaffold measurement but **no committed script
reproduces it**, so it cannot currently be re-run after an input change. Rebuilding it as a
`--native` flag on `interface.py` is the obvious fix and would make the control usable as intended --
now doubly needed, since the control has to be redone under relax regardless.

### The uOR clash is a chirality error in the input model (2026-10-30)

> **Correction.** This section first diagnosed the clash as a sidechain *rotamer* problem in the
> designs and recommended fixing the repack. **That was wrong.** The measurements below are correct
> and still worth reading, but the cause is one level deeper and is given under "Actual root cause":
> **residues 102-107 of the uOR nanobody are D-amino acids.** Do not act on the fix order that
> follows the table without reading to the end of this section.

Comparing the transplanted motif Tyr against native Y102 in `../uor_complex.pdb`:

| complex                       | Tyr->K238 | backbone vs native Y102 | **sidechain vs native** |
| ----------------------------- | --------- | ----------------------- | ----------------------- |
| `long__amide__ASN-GLY__c00`   | 1.69 A    | 0.37 A                  | **2.59 A** (max 3.25)   |
| `long__amide__DALA-GLY__c00`  | 1.26 A    | 0.15 A                  | **2.93 A** (max 3.64)   |
| `long__amide__DALA-GLY__c01`  | 1.33 A    | 0.50 A                  | **3.93 A** (max 5.01)   |
| `long__amide__DPRO-GLY__c00`  | 1.48 A    | 0.22 A                  | **2.83 A** (max 3.49)   |

Native contact is **CD1**-CD at 3.39 A; every design contacts through **CE2**. The backbone lands
within 0.15-0.50 A of native -- the CA superposition is working exactly as intended -- while the ring
has rotated into a different chi. The clash originates in the design's own rotamer.

**Root cause: the designs are built and relaxed as peptide-only poses with the receptor absent**
(see this script's own docstring, "chain B is never in the pose"). The motif sidechains therefore
optimise against solvent in a free macrocycle, settle into rotamers that are perfectly reasonable for
an isolated peptide, and only then does the receptor get pasted in around them. The transplant does
not so much *create* the clash as **reveal that the motif sidechains were never optimised in the
presence of the thing they are meant to bind.**

Two details in `cyclize/validate.py` let this through unnoticed:

- `superimpose_on_motif` fits on **CA only** (`atom='CA'`, `validate.py:39`). For a transplant whose
  entire purpose is reproducing specific sidechain contacts, nothing constrains the sidechains.
- `min_receptor_distance` is called with `skip_nums={motif residues}` (`validate.py:96`), so the
  2.2 A clash filter is measured **only on residues the design added**. The skip is deliberate and
  its stated rationale is sound -- without it, native contacts like Y102 get flagged as clashes and
  all 50 structures fail. But the rationale assumes the motif keeps its native sidechain geometry,
  which is precisely what breaks here. **A motif Tyr at 1.08 A is invisible to validation by
  construction**, which is why 68 of 68 passed.

So the fix order is cheapest-first, and a restrained relax is *not* the first thing to try:

1. **Check why `pack_input=true` is not already fixing a Tyr chi flip.** The emitted XML uses the
   default rotamer library -- no `-ex1 -ex2`, no TaskOperations. A repack able to reach the native
   rotamer should resolve most of this for free, with no relax and no peptide movement. Test on one
   structure; if extra rotamers fix it, this is a packer-resolution bug and everything downstream
   gets simpler.
2. **If that is not enough, restrain motif *sidechain* atoms to native during the relax**, not just
   CA. That is the targeted version of 50216944 and it addresses that job's open caveat directly:
   it relieves the clash without letting the peptide walk off the pharmacophore.
3. **For a real ddG, score the native peptide through the identical pipeline and report
   design - native.** `dG_separated` is explicitly not a ddG. A native reference is what cancels the
   uOR model's intrinsic strain and makes its numbers comparable to apj -- worth considerably more
   than any change to the alignment target, and it is the `--native` flag already on the backlog.

Items 1-2 patch this at **scoring** time, which is the right place to start because it is cheap and
changes nothing upstream. The **durable** fix is further back: relax the final designs with the
receptor present, so the motif sidechains are optimised against the surface they have to bind rather
than against solvent. That is a change to the build, not to `interface.py`, and it should wait until
item 1 has established how much of the problem is merely packer resolution.

### Actual root cause: uOR residues 102-107 are D-amino acids

**The fix order above does not apply to uOR.** No packer setting and no relax can fix this, because
the coordinates being matched are mirror-imaged. Verified 2026-10-30 by the improper dihedral
N-CA-C-CB, whose sign is the backbone's chirality:

```
uor_complex.pdb chain A        other chains (controls)
  100 VAL  +57.7                uOR receptor chain B   0 / 277 inverted
  101 SER  +58.1                apj nanobody chain A   0 / 107 inverted
  102 TYR  -50.3  <-- D         9v32 Fab H / L / A     0 / 545 inverted
  103 ASN  -53.5  <-- D
  104 PHE  -53.7  <-- D        L-residues in this chain: +57.6 +- 1.3 deg
  105 PHE  -64.2  <-- D        (sign convention fixed against experimental
  106 PHE  -91.7  <-- D         9v32, where every residue is L by definition)
  107 ARG  -53.7  <-- D
  108 SER  +57.8
```

Six consecutive residues, **exactly the CDR3 tip, exactly the transplanted motif**, with a sharp
boundary at 101/102 and 107/108. Bond lengths through the region are normal (N-CA 1.48-1.51, CA-C
1.49-1.55, CA-CB 1.52-1.53 A), so this is not broken connectivity -- it is a mirror-imaged backbone,
almost certainly a local Coot refit artifact. The source filename contains "flip" twice:
`../Rosetta_9_real_trunc_model_1_flip_fit_J126map-coot-0_flip.pdb`. `../prep/nanobody_gpcr_prep.md`
checked occupancies, altlocs, missing atoms and bond lengths, but never chirality.

**It propagated into every uOR scaffold**: `../uOR_scaffolds/{short,medium,long}.pdb` each carry all
six (short is 8 residues of which **6 are D**). apj scaffolds: 0.

**How that produces the clash.** Rosetta builds the designs correctly as **L** -- in
`uor/interface/complexes/long__amide__ASN-GLY__c00.pdb`, 14 of 15 residues are +57. The design is
then CA-superimposed onto D target coordinates. Backbone matches; CB lands on the opposite side:

| native | dCA | dCB | |
| ------ | --- | --- | --- |
| 100, 101 | 0.52, 0.49 | 0.88, 0.46 | L in native -- normal |
| 102 TYR | 0.42 | **2.02** | D in native, rebuilt L |
| 103 ASN | 0.25 | **2.33** | D in native, rebuilt L |
| 104 PHE | 0.28 | **2.22** | D in native, rebuilt L |
| 105 PHE | 0.34 | **2.59** | D in native, rebuilt L |
| 106 PHE | 0.14 | **0.24** | **D in native, and the design kept it D (-91.6 vs -91.7)** |
| 107 ARG | 0.28 | **2.68** | D in native, rebuilt L |
| 108-111 | 0.18-0.57 | 0.07-0.61 | L in native -- normal |

**Residue 106 is the control that settles it.** It is the one motif residue the build left as D, and
it is the one motif residue whose CB does *not* move -- 0.24 A against 2.0-2.7 A for the five that
were rebuilt as L. Nothing but a chirality mismatch produces that pattern. The consequence at the
interface: the design's Tyr reaches 1.01 A of LYS238:NZ, where native Y102's closest receptor
contact is an ordinary 2.48 A.

The per-atom gradient on `long__amide__DALA-GLY__c00` vs native Y102 makes the same point in one
line, and rules out the rotamer explanation outright:

```
N 0.23   CA 0.17   C 0.09   O 0.09   |   CB 2.21   CG 2.56   CZ 3.26   OH 3.64
```

**CB is not a rotamer degree of freedom** -- its position is fixed by the backbone. A 2.21 A CB
displacement sitting behind a 0.17 A CA is the signature of inverted chirality and of nothing else;
no chi angle can move CB at all. The growing displacement out along the sidechain (2.21 -> 3.64) is
just that error being levered outward.

**So every uOR anomaly recorded in these docs has one cause.** The `rama_prepro` of 33-37 against
apj's 10-14 is D residues scored against L Ramachandran tables, not "global strain". The native-pair
control at +34 to +45 "scoring as a non-binder" is the same six residues, not a bad receptor. The
68-of-68 clash is structurally unavoidable rather than a packing failure. And the earlier reading on
this page -- that the designs' sidechains were simply never optimised against the receptor -- is
**wrong**: the designs are right, the target is wrong.

**What this means for the work.** uOR designs cannot be rescored into validity; they were pinned to
mirror-imaged targets and need rebuilding. The repair is cheaper than the "replace the uOR input
model" item it supersedes -- it is a six-residue local geometry fix, not a new structure -- but it
sits upstream of everything: repair the model, re-run the `../flexddg/uor_nb/` alanine scan (the
hotspots that **selected this motif** were measured on D residues, so Y102's hotspot status is not
yet trustworthy), re-excise the scaffolds from the repaired pose, then rebuild. **apj is untouched
by all of this** and the scoring fixes above still apply to it.

One caution on the repair itself: reflecting CB across the N-CA-C plane is the minimal correction,
but the sidechains then need rebuilding and the loop needs refinement against the original J126
density, which is what the hand should ultimately be decided by. Treat "idealize and refit" as the
shape of the fix, not as a one-command operation.

#### Rosetta has been silently re-imposing L geometry all along (tested 2026-10-30)

A CA-restrained `FastRelax` (3 repeats, sd 0.5 A) on chain A alone, scored before and after
(`scripts/relax_uor_test/`):

| | `total_score` | `rama_prepro` | `omega` | D-residues |
| --- | --- | --- | --- | --- |
| input | +45.3 | **31.96** | 3.01 | **6** |
| relaxed | -339.1 | **9.50** | 9.29 | **1** |

At face value that looks like a repair -- `rama_prepro` falls by 22 REU and five of the six D
residues read out as L afterwards. **It is not a repair, and minimization did not do it.** Energy
minimization cannot invert a stereocentre; that needs bonds broken and reformed. What actually
happens is that the PDB labels these residues with standard L three-letter codes (`TYR`, `ASN`,
`PHE`, `ARG`), so **Rosetta assigns L residue types on read and rebuilds CB and the sidechain from
ideal L internal coordinates**, silently discarding the input's hand.

The cost shows in what moved -- backbone stays put, sidechains swing:

| res | N | CA | C | CB | sidechain tip |
| --- | --- | --- | --- | --- | --- |
| 102 TYR | 0.36 | 0.28 | 0.69 | **1.94** | 2.14 |
| 103 ASN | 0.49 | 0.45 | 0.42 | **1.24** | 3.20 |
| 104 PHE | 0.42 | 0.38 | 0.33 | **2.39** | **6.65** |
| 105 PHE | 0.72 | 0.68 | 0.68 | **2.42** | 5.43 |
| 106 PHE | 0.65 | 0.67 | 0.72 | 0.59 | 2.61 |
| 107 ARG | 0.42 | 0.42 | 0.45 | **2.44** | **6.70** |

**This is the mechanism behind every uOR anomaly, now observed directly and in isolation.** Rosetta
L-ifies the motif on read; the designs are therefore built L; they are then CA-superimposed back
onto the original **D** coordinates; CB lands 2.0-2.7 A away and the sidechain points into the
receptor. That is the 68-of-68 clash.

Residue **106 PHE resisted** -- improper unchanged at -88.3, CB moving only 0.59 A, with no missing
atoms to explain it. -88.3 is neither clean L nor clean D, so it is distorted in some further way
the rebuild did not resolve. It is also the one motif residue whose CB does not move in the built
designs (0.24 A), which is the control described above.

**Practical consequences.**

1. **A relax cannot decide the hand, and will always silently choose L** because of the residue
   naming. Do not use it as the repair, and do not read the `rama_prepro` drop as evidence of one.
2. **The hand must be settled against the original J126 density.** The two cases need different
   fixes: if the Coot backbone is mis-built and the protein is really L, Rosetta's L-ification
   points the right way, but the sidechain positions it invents are packer guesses sitting 2-7 A
   from the deposited ones and must be re-refined against density rather than accepted. If the
   density instead supports the deposited sidechains, the rebuild is wrong and the backbone has to
   be rebuilt to match.
3. **If a D residue is ever intended, it must carry a D residue type** (`DTY`, `DPH`, `DAR`, ...).
   Under standard L names Rosetta converts it without warning -- which is exactly how this survived
   an alanine scan and two full design sweeps unnoticed.

**Not a fix: re-aligning the designs onto a relaxed input pose** -- *as scoped here, meaning
relaxing the receptor while keeping the ORIGINAL scaffolds.* The 2026-10-30 rebuild does something
different and is not contradicted by this paragraph: it re-excises the motif **from the relaxed
pose itself**, so the scaffold and the alignment target are the same object and the D-vs-L
mismatch disappears by construction. The point below stands only against the narrow version.
The clash is on the design side;
relaxing the input moves K238 and native Y102 by a few tenths of an A, which cannot absorb a 2.6-3.9 A
sidechain displacement, and there is nothing in the input to repair anyway -- it is clean at 2.27 A.

**But one consistency requirement does follow.** The alignment target and the receptor must come from
the same pose. They currently do: `../uOR_scaffolds/{short,medium,long}.pdb` are excised verbatim
from `../uor_complex.pdb` chain A -- 8/8, 10/10 and 14/14 CAs at **0.000 A** deviation. If the uOR
receptor is ever replaced, **the scaffolds must be regenerated from that same replacement pose** or
every design inherits a systematic offset.

So: **every committed `interface.csv` number is an as-is number and therefore suspect**, uOR's
catastrophically and apj's by 15-25 REU. The fix is to port this sbatch's protocol into
`interface.py` (keeping the route-aware `PeptideCyclizeMover` guard) and rescore all four targets.
Run the native control after any change to an input complex -- it costs seconds and it is the only
thing that distinguishes a bad design from a bad receptor.

## Use

```bash
bin/cyclize register apj          # hairpin register report, writes nothing
bin/cyclize setup    apj          # enumerate variants -> inputs/, variants.txt, setup.log
bin/cyclize smoke    apj long__amide__DPRO-GLY --nstruct 3   # one variant in $SCRATCH/cyclize_smoke/
bin/cyclize submit   apj          # one array task per variant + a validation job after
bin/cyclize status   apj          # per-variant closure state and active jobs
bin/cyclize validate apj          # analysis_output/{results,summary}.csv
```

Analysis scripts, run from `cyclize/`:

```bash
python3 scripts/degeneracy.py --csv-dir=scripts   # how many of the passers are actually distinct
python3 scripts/plot_degeneracy.py                # scripts/degeneracy.png
python3 scripts/representatives.py                # <target>/representatives/, one PDB per cluster
ml biology rosetta/3.14
python3 scripts/rescore.py apj uor                # <target>/rescore/rescore.csv
python3 scripts/interface.py apj uor              # <target>/interface/interface.csv
python3 scripts/compare_prebug.py                 # pre/post rm_lower-fix geometry, reports only
```

All of the above over every target at once, as a job rather than in an `sh_dev` session -- the
interface scoring is ~3.6 s/structure and the four targets carry several hundred representatives:

```bash
sbatch -D $PWD -J cyclize_distill -p normal,possu -t 04:00:00 -c 1 --mem 8G \
       -o logs/distill_%j.out scripts/distill.sbatch
```

Order matters inside it: `rescore` and `interface` both read
`<target>/representatives/representatives.csv`, so `representatives` runs first.

The i+1 sweep lives in its own target directories, `apj_i1/` and `uor_i1/`, rather than in `apj/` and
`uor/`. `setup` rewrites `variants.txt` and `inputs/`, so running it in place with a different
`redesign_i1` would destroy the record of the finished 55-variant grid.

`submit` refuses while the target has active jobs, and does not rerun partial/empty variants unless
you pass `--retry` after `status` has shown you why. Each submission is recorded in
`tasks/submissions.tsv`; log `logs/<job>_<i>.out` ran line `i` of that task file.

One array task runs `slurm.per_task` variants (default 20), each running all `nstruct` in a single
Rosetta process. Rosetta's own nstruct loop is far cheaper than a process per structure, but it is not
enough on its own: measured across the first sweep, **one variant takes about a minute**, so a task
per variant was mostly queue latency and sat well under the 10-minute floor that makes a job worth
submitting. One line of the task file is one task; a line holding a single id is still valid, so older
task files replay unchanged. Seeds come from a variant's position in the whole file rather than within
its line, so chunking never collides two RNG streams. A variant that throws is logged and skipped
rather than taking the rest of its task down.

Rosetta exits non-zero whenever *any* structure failed to close, which is normal for constrained
rings, so `runner.state` counts output PDBs rather than trusting the exit code. That count is a bare
`glob('*.pdb')` on `output/<vid>/`, so **no later stage may write PDBs there** -- it would be counted
as closure and corrupt `status`.

## cyclize.json

```json
{"title": "APJ d006 CDR3 macrocycles",
 "complex": "../../apj_complex.pdb",
 "motif_chain": "A",
 "receptor_chain": "B",
 "scaffolds": {"short": "../../apj_scaffolds/short.pdb"},
 "scan": "../../flexddg/apj_nb/analysis_output/alascan_summary.csv"}
```

Optional sections and their defaults (see `cyclize/config.py`; unknown keys are rejected):

- `hotspot`: 1.0 kcal/mol. Residues above it are never released to GenKIC and never mutated to Cys.
- `routes`: `["amide", "ssapp", "ssmut"]`
- `redesign_i1`: false. True sweeps every legal identity at i+1 instead of using `TURN_SEEDS`.
- `protocol`: nstruct 50, closure_attempts 2000, min_solutions 1, scorefxn ref2015, relax_rounds 1,
  rama_cutoff null, seed_rama_cutoff null, coord_dev 0.5, design_i2 false, target_loop 5,
  **omega_trans true**, omega_weight null, extra_flags []
  - **`omega_trans`** (added 2026-10-30, default **on**) emits a GenKIC `set_dihedral` perturber
    fixing every sampled peptide bond to omega 180. This is the fix for the defect that
    invalidated the first design set: with CA pivots the closure solves phi/psi at the pivots and
    leaves **omega free** (see the note above `variants.TARGET_LOOP`), and `PeptideCyclizeMover`
    constrains only the closure bond -- so 47% of the bonds flanking the closure came out >30 deg
    from planar and 50% came out non-proline cis. Because omega is not a solved DOF, setting it
    holds through the closure solve. The disulfide bond is skipped; only peptide bonds are pinned.
    Turn it off only to reproduce the old behaviour.
  - **`target_loop`** raises the number of GenKIC loop residues by releasing additional **cold**
    terminal motif residues (`pick_release` never releases a hotspot). Released residues keep
    their identity -- only their backbone is sampled, nothing is designed -- and the motif still
    gains only the 2 residues the route adds. Raise it only if `omega_trans` alone cannot close.
- `validate`: clash 2.2 A, motif_rmsd 1.0 A
- `slurm`: partition normal,possu, time 02:00:00, mem 2G, max_concurrent 100, per_task 20, analysis_*
- `rosetta_module`: `"biology rosetta/3.14"`

**Use `rosetta/3.14`, not `/home/groups/possu/Rosetta`.** The module has `rosetta_scripts` and
`simple_cycpep_predict` prebuilt and supports `terminal_disulfide` and `thioether_lariat`; the possu
checkout is older (no `thioether`/`1_4_BBMB` crosslinkers) and only has `rosetta_scripts` built. The
possu tree is still the place to *read* source, offline docs and integration tests.

## Layout

- `cyclize/`: `struct` (PDB + geometry), `hairpin` (register), `variants` (enumeration + seed rules),
  `build` (XML generation), `inputs` (setup), `runner` (one variant), `slurm`, `validate`
  (superposition + filters), `analysis`, `__main__` (CLI).
- `tests/test_cyclize.py`: `python3 tests/test_cyclize.py` (~2 s, no Rosetta needed). Pins the
  register, which routes each scaffold is eligible for, the phi-sign turn rules, and the pose index
  arithmetic in the generated XML. A real Rosetta run is not covered -- use `cyclize smoke`.
- Python: `/home/groups/possu/miniconda/bin/python3` (3.8, numpy 1.24), set in `bin/cyclize`;
  override with `CYCLIZE_PYTHON`.

## Known caveats

- **The closure carbonyl oxygen is misplaced in 57% of amide designs, and no metric can see it.**
  At the ring-closure bond of `apj medium__amide__DPRO-GLY__c00`, `O...N` is **0.38 A** against an
  ideal **2.25 A**, with `O=C-N` at **16.3 deg** against 123 -- the oxygen sits on top of the next
  residue's nitrogen. (That structure is the geometry-ranked inspection set's top pick, *not* the
  dG-best `short__amide__DPRO-GLY` named earlier on this page, whose closure is a correct 2.23 A
  and which fails on flanking omega instead.) It is purely a closure defect -- all **2627** internal
  peptide bonds measure 2.23-2.29 A -- and it affects **136 / 239** amide complexes (median
  1.85 A, min 0.14 A), in the raw build output as well as the scored complexes.

  **The closure is paid for by distorting its neighbours.** By ring-distance from the closure
  bond: distance 0 is pristine (0% non-planar, 0% cis), distance 1 is **47% non-planar and 50%
  non-proline cis**, distance 2 is 24% / 43%, and by distance 3 it is clean again.
  `PeptideCyclizeMover` constrains the closure `C-N` and nothing else, so GenKIC buys the ring by
  flipping the flanking bonds cis and twisting them out of plane.

  **Why nothing caught it:** `O` and `N` are a **1-3 pair through `C`**, which `fa_rep` excludes --
  so a 0.14 A overlap costs nothing in score (`corr(closure O...N, total_score) = +0.10`). The
  intra-peptide clash scan skips the closure pair by construction. `validate.py` tests motif RMSD,
  seed phi and receptor distance, none of which involve the oxygen. It was found by eye.

  **Combined with the omega and self-intersection criteria, 0 of 311 scored complexes are clean**
  -- 75 of 239 amides have a correct closure oxygen, 1 of those also avoids non-proline cis next
  to the closure, and none of those has planar flanking bonds. See `../STATUS.md` "Recommendations":
  filtering cannot fix this and the closure protocol needs replacing.

- **The cyclisation bond is never declared in the output PDBs, so every design *looks* open.**
  Rosetta writes no `LINK` or `CONECT` record for the ring closure, and viewers only auto-bond
  sequentially numbered residues -- so the last-residue `C` to residue-1 `N` bond is not drawn and
  the macrocycle renders as a linear peptide whose ends happen to sit 1.33 A apart. **This is a
  display artifact and is the first thing to rule out when a design looks broken on inspection.**
  The underlying geometry is textbook in every amide design checked: closure `C-N` **1.33 A**,
  `CA-C-N` **116.2 deg**, `C-N-CA` **121.7 deg**. `scripts/inspect_set.py` now appends `CONECT`
  records and emits an `inspect/view.pml` that issues explicit PyMOL `bond` commands, since PyMOL
  ignores `CONECT` for `ATOM` records under its default `connect_mode`. Disulfide (`ssapp`,
  `ssmut`) routes close through S-S and are drawn correctly without help.

- **Omega is not planar: the amide bonds are broken, and that is why the rings self-intersect**
  (found 2026-10-30; the most serious problem on this page, and it affects **both** GPCR targets).

  Bond *lengths* are correct -- every peptide bond in `apj short__amide__DPRO-GLY` measures
  1.31-1.34 A, and all bond *angles* are ideal (N-CA-C 108-114, CA-C-N 115-119, C-N-CA 119-126).
  What is wrong is the dihedral. Measuring omega (CA-C-N-CA; validated against experimental 9v32,
  which reads +175 to +179 throughout):

  ```
   1 VAL ->  2 LEU   -75.0   <-- NON-PLANAR
   2 LEU ->  3 GLY  -178.4   trans
   ...
   8 GLU ->  9 DPR    -0.1   cis (before D-Pro -- legitimate)
   9 DPR -> 10 GLY   +45.2   <-- NON-PLANAR
  10 GLY ->  1 VAL  +178.3   trans
  ```

  (One correction in the designs' favour: the companion "self-intersection" claim below used an
  `intra_min >= 3.4 A` cutoff that ignored atom type, and a 2.8-3.0 A N...O contact is a hydrogen
  bond rather than a clash. Of 143 apj structures, **80 genuinely clash** below 2.5 A, 56 are
  polar contacts at 2.5-3.0 A, and 10 are nonpolar at 2.5-3.0 A -- so ~56% self-intersect, not
  100%. The closure findings are unaffected and remain 0 of 311.)

  Across all **314** built complexes in all four targets (314 PDBs on disk; **311** carry a
  scored `interface.csv` row, which is the denominator for the validity counts below):

  | measure | value |
  | ------- | ----- |
  | complexes with >= 1 bond > 30 deg from planar | **281 / 314 (89%)** |
  | bonds non-planar (> 30 deg off) | 476 / 3716 (**13%**) |
  | worst deviation per complex | median 65 deg, **max 90 deg** |
  | cis bonds **not** before proline | **618** (vs 25 legitimate cis-Pro) |

  A 90-degree omega has no partial double-bond character -- it is not a peptide bond. And
  non-proline cis runs ~0.03% in real proteins against **17%** here.

  **With omega free, the chain can fold back in ways a real peptide cannot.** That is what produces
  the self-intersection: `CA8-CA10` **1.87 A**, `CA7-CA10` 2.40 A, `7:C - 10:N` **0.64 A**. GenKIC
  found a valid *closure* that is not a valid *structure* -- and it did so by breaking the amide.

  **Why nothing caught it.** Three independent gaps:
  1. GenKIC runs `loop_bump_check` and `rama_prepro_check` (`build.py:247` and `:252`) but
     **no omega filter** -- omega is not among the quantities it tests.
  2. `FastRelax` cannot repair it: the MoveMap sets `bb="0"` globally and enables `bb="1"` only on
     loop spans, so **motif-residue omega is frozen** at whatever the closure produced, and the
     motif is additionally pinned by CA coordinate constraints.
  3. The scorefunction *did* flag it -- build-time `score.sc` reads `omega` **34.7** and
     `rama_prepro` **40.9** -- but `fa_rep` **3278.8** swamps them in `total_score`, and nothing
     downstream filtered on any of these.

  The cheap diagnostic fix is to enable omega minimization on motif residues (`bb="1"` plus the
  existing CA constraints, rather than freezing the backbone) and add a post-GenKIC omega filter.
  Whether that leaves any surviving closures is the open question in `../STATUS.md` under
  "Recommendations" -- if a 2-residue turn cannot close with planar amides, this is a dead end
  rather than a bug, and the `i + j = 209` register is what has to give.

  Across all 143 scored apj structures the worst non-adjacent intra-peptide contact has median
  **2.43 A**, p25 **1.39 A**, minimum **0.38 A**, and the largest is **2.99 A**. Per the
  correction above, count genuine self-intersection as the **80 of 143** below 2.5 A; the 56
  structures whose closest contact is a 2.5-3.0 A N...O pair are making hydrogen bonds.

  It is in the builds, not the scoring: `apj/output/short__amide__DPRO-GLY/short_0001.pdb` has the
  identical geometry and its build-time `score.sc` reads `total_score` 3358.7 with **`fa_rep`
  3278.8**. `FastRelax` cannot relieve it because the motif is pinned and the ring cannot open.
  This is the origin of the multi-thousand-REU `total_score` values seen throughout the interface
  tables above. Decomposed on one complex: receptor alone -270.5 (`fa_rep` 350.8), **peptide alone
  +3415.8 (`fa_rep` +3337.1)**.

  **The amide bond was the obvious suspect and it is innocent -- do not re-investigate it.** Scoring
  the same complex three ways gives, to the decimal:

  | treatment | total |
  | --------- | ----- |
  | raw read, no cyclization | 4009.6 |
  | `PeptideCyclizeMover` only (what `interface.py` does) | **3236.5** |
  | `ModifyVariantType` x2 + `PeptideCyclizeMover` | **3236.5** |

  `PeptideCyclizeMover` strips the terminus variants itself, so the `ModifyVariantType` pair that
  `build.py` emits (lines 181-183, guarded on `strip_termini`, with a comment saying it is
  required) being **absent** from `interface.py`'s XML is harmless. Cyclization is working and is
  worth 773 REU.

  **`dG_separated` is structurally blind to it.** Intra-peptide energy appears in both the complex
  and the separated-chains term and cancels exactly, so a self-intersecting peptide pays no dG
  penalty while its collapsed mass buries extra surface and gains `fa_atr`. The ranking was
  selecting *for* the defect -- `corr(intra, total_score) = -0.61`, and all seven top designs by dG
  sit at 0.38-0.69 A. Never rank on `dG_separated` alone again.

  The disulfides are less bad but not clean: `short__ssapp__DCYD-DCYD` runs 1.72-1.93 A. When
  checking `ssapp` by hand, note the Cys pair is non-adjacent in sequence so the S-S bond can
  masquerade as a clash -- it is not the problem here (SG-SG is a legitimate 2.37-2.52 A; the
  offenders are backbone overlaps such as `2:O - 10:CB` at 1.93 A).

  **Root cause in this package:** `validate.py` gates on motif RMSD, seed phi, and
  design-versus-*receptor* distance via `min_receptor_distance`. **It never checks the design
  against itself.** Combined with that function's motif-skipping (see the chirality section), there
  are two independent blind spots in the same validator, and each one let a whole target's worth of
  bad structures through. Add an intra-peptide clash gate before rebuilding anything.

  Open design question this raises: a 2-residue turn may be too short to span the hairpin gap at
  all, in which case the closures are geometrically infeasible rather than badly sampled -- which
  would tie straight back to the `i + j = 209` register that forces even-length motifs and 0- or
  2-residue closures. Settle that before spending compute on a rebuild.

- **`ssapp` stripped a terminus it never bonds, costing a ~10 REU phantom on residue 1. FIXED
  2026-10-29.** `build.xml` emitted `rm_upper`/`rm_lower` unconditionally, so an `ssapp` protocol
  removed `LOWER_TERMINUS_VARIANT` from residue 1 and then closed the ring through S-S instead,
  leaving residue 1 with nothing bonded to its N. Rosetta evaluates its rama at phi = 0, a forbidden
  region, and invents about **+10 REU that is not real geometry** (measured 9.993 on
  `long__ssapp__CYD-CYD`). The movers are now emitted only when the ring closes through a peptide
  bond (`close['atom1'] != 'SG'`), and all 25 disulfide variants were rerun. Pre-fix structures are
  kept under `<target>/output_prebug/`; `scripts/compare_prebug.py` measures what actually changed.
  **Settled 2026-10-29 (job 50173650): the artifact moved scores, not geometry** -- psi(1) shifted
  by +3.18 deg on apj and -3.48 deg on uOR, against a within-variant spread of ~100 deg, and the
  S-S distance by <=0.01 A. The reasoning that made it look like it might have moved geometry, kept
  because it is the right thing to check next time: phi(1) is fixed garbage that minimisation
  cannot touch, but psi(1) is real, free and inside the relax MoveMap, so relax *should* have
  optimised it against a bogus phi -- and residue 1 is the prepended Cys, whose backbone places SG.
  When rescoring, read the PDB as-is. **The fix landed the same day**: the two movers are guarded
  on the closure chemistry in `build.py` (`strip_termini`).
- **`ssmut` produces a strained disulfide** (measured 1.70 A S-S against an ideal 2.05) because the
  motif is pinned, so the bond cannot relax into place. Loosening `protocol.coord_dev` or allowing the
  two Cys backbones to move would fix it, at the cost of some motif rigidity.
- `design_i2` is accepted in the config but **not yet implemented**: i+2 of the turn is currently
  always taken from the seed pair rather than designed. Chirality is only load-bearing at i+1, so
  designing i+2 freely among L residues is the obvious next step -- and now that `redesign_i1` exists
  it should be built the same way, by enumeration through `sweep_seeds`, not by adding a packer.
- **~~The uOR input model is globally strained~~, *and* the Y102-K238 clash is real -- in the
  built complexes, not in the input.** *(The "globally strained" half is superseded: the strain
  signal was D residues scored against L tables. The clash half stands and its mechanism is the
  chirality mismatch -- see "Actual root cause" above.)* The two were conflated twice, in both directions. Settled
  2026-10-29:
  - In `../uor_complex.pdb` the closest Y102-K238 heavy-atom pair is **3.39 A** (CD1-CD) and the
    whole A/B interface has nothing under 2.0 A (closest 2.27 A). The input pose is clean.
  - In `uor/interface/complexes/` the clash is back and universal: **68 of 68** complexes have an
    interchain pair under 2.0 A, mean closest approach **1.19 A**, and the transplanted motif Tyr
    sits **1.08 A** from LYS238:CD. apj is 26 of 55 at a mean 1.80 A. So `validate.superimpose_on_motif`
    reintroduces it -- the 2.2 A clash filter is measured on the motif, not on every transplanted
    sidechain.

  **Job 50216944 answered the rest on 2026-10-30: the clash is fully recoverable and the model is
  not the main culprit.** A CA-restrained FastRelax clears all six uOR test structures (0.77-1.45 A
  -> 2.45-2.88 A) and flips `dG_separated` from +367/+913 to **-14 to -33**, while the three apj
  controls barely move on geometry. So the ~400 REU was `interface.py` never relieving backbone.
  Two consequences: the relax arm is *not* neutral (apj improves 15-25 REU too, and its top two
  reorder), and the relief is bought by moving the peptide up to 1.26 A CA RMSD, so **post-relax
  motif RMSD is the missing measurement**.

  The **native**-pair control -- uOR `dG_separated` +34 to +45 and `total_score` +180 to +209 against
  apj's -6 to -19 and ~-240 -- was read as a **floor** under every uOR number. That reading did
  **not** survive: the relaxed designs score below the supposed floor, so the control was itself an
  as-is measurement dominated by clash. Its real size is unknown until it is re-run under relax. The
  ~~model is still strained on independent evidence (`rama_prepro` 33-37 vs apj's 10-14)~~ --
  **that reading did not survive either**: the 33-37 gap is the D residues being scored against L
  Ramachandran tables, and a CA-restrained relax of chain A takes it to 9.50 with no refitting. Any
  residual strain is worth perhaps 10-25 REU, not 400. See "Interface score" above.
- **The RNG seed is the variant's position in the task file, so resubmitting can silently rerun the
  identical trajectory.** `slurm.run_task` computes `seed = first + i + 1` from the variant's offset
  in the whole task file and `runner.command` passes `-constant_seed -jran $((1000 + seed))`, which
  makes Rosetta fully deterministic. The intent (in the code comment) is that each variant keeps a
  distinct stream however the file is chunked; the side effect is that **the same variant at the same
  position in a later submission draws the same seed and reproduces its output bit-for-bit**.

  This bit once, on 2026-10-29: the uOR retry and the uOR `rm_lower` rerun had task files sharing
  their first eight entries, so those eight variants came out 50/50 byte-identical and looked like a
  failed rerun. They were not -- distinct inodes, genuinely recomputed, identical because the seed and
  the sampling were both unchanged. apj escaped only because its retry file omitted one variant and
  shifted all the rest. Useful accident: those eight are the project's one *paired* pre/post
  comparison, and they show the terminus fix changed nothing about the sampling.

  If you want a controlled comparison, pin the seed deliberately; if you want an independent sample,
  change the ordering. Do not assume a resubmission gives you either.
- Hotspot coverage is only 60-78% of the alanine-scan signal for these motifs; the rest sits in CDR1
  and CDR2 (APJ loses R27 + Y33 = 5.6 REU, uOR loses E54 + L32 = 3.0 REU) and no CDR3-only macrocycle
  can recover it. See `../flexddg/README.md`.
