# Cyclising CDR3 hairpin motifs into macrocyclic peptides

`cyclize/` is the package; each target is a directory with a `cyclize.json`. It takes the CDR3 motifs
selected from the flex ddG alanine scans (`../flexddg/`) and closes them into cyclic peptides by amide
or disulfide, holding the ddG-validated binding conformation rigid.

| dir     | complex                | scaffolds                   | receptor | scan source          |
| ------- | ---------------------- | --------------------------- | -------- | -------------------- |
| `apj/`  | `../apj_complex.pdb`   | `../apj_scaffolds/*.pdb`    | chain B  | `../flexddg/apj_nb/` |
| `uor/`  | `../uor_complex.pdb`   | `../uOR_scaffolds/*.pdb`    | chain B  | `../flexddg/uor_nb/` |

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

The motif is located **by position in the design, not by residue number**. Routes that add residues
make Rosetta renumber the pose 1..N, but `ssmut` adds nothing so the input PDB numbering survives
(101-108). A numeric range matches one and not the other, which made uOR `ssmut` report its own motif
hotspot Y102 as an "added residue" clashing with the receptor and fail all 50 structures. For `ssmut`
the reported receptor distance is `inf`, correctly: it adds no residues, so there is nothing to check.

## Clash threshold

`validate.clash` defaults to **2.2 A, not 3.0**. The native APJ interface packs motif residues to
2.44-2.80 A of the receptor (F104 2.44, E108 2.65, K110 2.70), so a 3.0 A cutoff rejects the native
structure itself. 2.2 A is below anything observed natively, so it flags only real overlap.

## Use

```bash
bin/cyclize register apj          # hairpin register report, writes nothing
bin/cyclize setup    apj          # enumerate variants -> inputs/, variants.txt, setup.log
bin/cyclize smoke    apj long__amide__DPRO-GLY --nstruct 3   # one variant in $SCRATCH/cyclize_smoke/
bin/cyclize submit   apj          # one array task per variant + a validation job after
bin/cyclize status   apj          # per-variant closure state and active jobs
bin/cyclize validate apj          # analysis_output/{results,summary}.csv
```

`submit` refuses while the target has active jobs, and does not rerun partial/empty variants unless
you pass `--retry` after `status` has shown you why. Each submission is recorded in
`tasks/submissions.tsv`; log `logs/<job>_<i>.out` ran line `i` of that task file.

One array task is one variant running all `nstruct` in a single Rosetta process -- Rosetta's own
nstruct loop is far cheaper than a process per structure and keeps each task above the 10-minute floor
that makes a job worth submitting. Rosetta exits non-zero whenever *any* structure failed to close,
which is normal for constrained rings, so `runner.state` counts output PDBs rather than trusting the
exit code.

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
- `protocol`: nstruct 50, closure_attempts 2000, min_solutions 1, scorefxn ref2015, relax_rounds 1,
  rama_cutoff null, seed_rama_cutoff null, coord_dev 0.5, design_i2 false, extra_flags []
- `validate`: clash 2.2 A, motif_rmsd 1.0 A
- `slurm`: partition normal,possu, time 02:00:00, mem 2G, max_concurrent 100, analysis_*
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

- **`ssmut` produces a strained disulfide** (measured 1.70 A S-S against an ideal 2.05) because the
  motif is pinned, so the bond cannot relax into place. Loosening `protocol.coord_dev` or allowing the
  two Cys backbones to move would fix it, at the cost of some motif rigidity.
- `design_i2` is accepted in the config but **not yet implemented**: i+2 of the turn is currently
  always taken from the seed pair rather than designed. Chirality is only load-bearing at i+1, so
  designing i+2 freely among L residues is the obvious next step.
- **The uOR input model has a 1.09 A heavy-atom contact** between motif Y102 and receptor K238. That
  is a genuine clash in the Coot-built structure, not a design artifact -- `ssmut` designs reproduce the
  native backbone at 0.00 A RMSD and inherit it. Worth resolving in the input before trusting fine
  distance filtering on uOR.
- Hotspot coverage is only 60-78% of the alanine-scan signal for these motifs; the rest sits in CDR1
  and CDR2 (APJ loses R27 + Y33 = 5.6 REU, uOR loses E54 + L32 = 3.0 REU) and no CDR3-only macrocycle
  can recover it. See `../flexddg/README.md`.
