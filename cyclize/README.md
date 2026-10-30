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
  are.** Clustering the passing structures on N/CA/C RMSD (`scripts/degeneracy.py`) gave **112
  distinct conformations at 1.0 A** -- apj 51, uOR 61. The variants with the best pass rates are the
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
`rosetta_loop` from `input.json`. Output is `<target>/rescore/rescore.csv`; the 112 representatives
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
enough to put apj into sensible territory -- the 51 representatives run **-25.2 to positive**, with
`short__amide__DPRO-GLY` best at -25.2 REU over 1429 A^2 of buried surface. It is a **sidechain**
repack only; backbone clash from the superposition is not relieved, which is intended, since
`validate` has already rejected anything closer than 2.2 A. Cost is ~3.6 s/structure, so the 51
representatives take about three minutes.

**uOR interface numbers are not usable.** Its 61 representatives come out at **+370 to +413 REU**
against apj's -25. The cause is *not* one bad contact, and the Y102-K238 clash recorded under Known
caveats does not explain it -- measured in the current `uor_complex.pdb`, Y102:CD1 to K238:CD is
**3.39 A**, an ordinary van der Waals contact, with nothing below 3.3 A between them. Either that
1.09 A figure refers to an earlier file or it was repaired at some point; it is not in the input
being used. Treat the Known-caveats entry as unverified.

The real finding is worse, and comes from scoring the **native** motif:receptor pair -- the actual
nanobody CDR3 in its actual position, a real binder -- through the identical protocol:

| control | `dG_separated` | `ddg` | `total_score` |
| ------- | -------------- | ----- | ------------- |
| apj short / medium / long  | -5.7 / -12.9 / -19.1 | -4.8 / -15.5 / -13.3 | ~-240 |
| uOR short / medium / long  | **+34.4 / +41.7 / +45.2** | +97.3 / +91.6 / +93.8 | **+180 to +209** |

**The native uOR complex scores as a non-binder.** apj's natives sit at -5.7 to -19.1 REU with a
total score near -240; uOR's are positive on every measure. That is a global property of the
Coot-built model, not a single repairable contact, and it sets a floor under every uOR number: no
design can score well against a receptor whose own native complex does not. The designs are a
further ~330 REU worse again, which that floor does not explain and which has not been chased down.

So: **apj interface numbers are usable, uOR's are not**, and the fix is a better uOR input model
rather than a local repair. Run this control after any change to an input complex -- it costs
seconds and it is the only thing that distinguishes a bad design from a bad receptor.

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
  extra_flags []
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

- **`ssapp` stripped a terminus it never bonds, costing a ~10 REU phantom on residue 1. FIXED
  2026-10-29.** `build.xml` emitted `rm_upper`/`rm_lower` unconditionally, so an `ssapp` protocol
  removed `LOWER_TERMINUS_VARIANT` from residue 1 and then closed the ring through S-S instead,
  leaving residue 1 with nothing bonded to its N. Rosetta evaluates its rama at phi = 0, a forbidden
  region, and invents about **+10 REU that is not real geometry** (measured 9.993 on
  `long__ssapp__CYD-CYD`). The movers are now emitted only when the ring closes through a peptide
  bond (`close['atom1'] != 'SG'`), and all 25 disulfide variants were rerun. Pre-fix structures are
  kept under `<target>/output_prebug/`; `scripts/compare_prebug.py` measures what actually changed.
  Whether the artifact moved geometry or only the reported score is an **inference** until that
  comparison runs: phi(1) is fixed garbage that minimisation cannot touch, but psi(1) is real, free,
  and inside the relax MoveMap, so relax should have optimised it against a bogus phi -- and
  residue 1 is the prepended Cys, whose backbone places SG. When rescoring, read
  the PDB as-is. Fixing it means guarding those two movers on `route == 'amide'`.
- **`ssmut` produces a strained disulfide** (measured 1.70 A S-S against an ideal 2.05) because the
  motif is pinned, so the bond cannot relax into place. Loosening `protocol.coord_dev` or allowing the
  two Cys backbones to move would fix it, at the cost of some motif rigidity.
- `design_i2` is accepted in the config but **not yet implemented**: i+2 of the turn is currently
  always taken from the seed pair rather than designed. Chirality is only load-bearing at i+1, so
  designing i+2 freely among L residues is the obvious next step -- and now that `redesign_i1` exists
  it should be built the same way, by enumeration through `sweep_seeds`, not by adding a packer.
- **The uOR input model is globally strained, and the "1.09 A Y102-K238 clash" does not reproduce.**
  The original note here recorded a 1.09 A heavy-atom contact between motif Y102 and receptor K238.
  Re-measured 2026-10-29 in `../uor_complex.pdb`, the closest Y102-K238 heavy-atom pair is **3.39 A**
  (CD1-CD), entirely normal, so that figure either refers to an earlier file or was repaired. What is
  real and worse: the **native** uOR motif:receptor pair scores `dG_separated` +34 to +45 and
  `total_score` +180 to +209, where apj's natives score -6 to -19 and ~-240. The uOR model scores as
  a non-binder in its own native geometry, which puts a floor under every uOR number. See "Interface
  score" above. A better input model is needed, not a local repair.
- Hotspot coverage is only 60-78% of the alanine-scan signal for these motifs; the rest sits in CDR1
  and CDR2 (APJ loses R27 + Y33 = 5.6 REU, uOR loses E54 + L32 = 3.0 REU) and no CDR3-only macrocycle
  can recover it. See `../flexddg/README.md`.
