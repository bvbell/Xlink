# Archived: the 2026-10-30 invalid design set

Everything the `cyclize` pipeline produced before 2026-10-30 is archived and **superseded**. None
of it is chemically valid. This note records what was archived, where it is, and exactly why it
failed, so the set can be re-examined without re-deriving the diagnosis.

## Where it is

```
$OAK/bvbell_archive/xlink_cyclize_invalid_geometry_20261030.tar.gz     72 MB, 10274 files
```

A single tarball deliberately: OAK was at 93% of its inode quota when this was written, and
`$GROUP_HOME` at 99% of volume, so loose files were not an option. Contents, for all four targets
(`apj/`, `uor/`, `apj_i1/`, `uor_i1/`):

```
<target>/output/            the raw GenKIC + FastRelax trajectories
<target>/interface/         built two-chain complexes and interface.csv
<target>/representatives/   clustering representatives
<target>/analysis_output/   summary.csv and per-variant closure counts
<target>/cyclize.json       the config that produced them
<target>/variants.txt       the variant list
```

`*.silent` files were excluded (regenerable, and large). The originals remain under
`cyclize/<target>/` in `$SCRATCH` until the 90-day purge takes them; the tarball is the durable
copy.

## What it contained

| set | variants | structures | passed the old validation |
| --- | --- | --- | --- |
| original grid (`apj/`, `uor/`) | 55 | 2750 | 1205 |
| i+1 sweep (`apj_i1/`, `uor_i1/`) | 126 | 6300 | 2483 |

Scored interface complexes: 311 (143 apj + 168 uOR).

## Why all of it is invalid

**0 of 311 scored complexes are chemically valid.** Three independent defects, none of which any
filter in the pipeline tested for.

### 1. The closure carbonyl oxygen is misplaced (amide routes)

At the ring-closure bond the `O···N` distance runs as low as **0.14 Å** against an ideal 2.25 Å,
with `O=C–N` angles down to **16°** against 123°. The oxygen sits on top of the next residue's
nitrogen. **136 / 239 amide complexes (57%)**, median 1.85 Å.

It is purely a closure defect: all **2627** internal peptide bonds measure 2.23–2.29 Å.

Invisible to scoring because **`O` and `N` are a 1-3 (geminal) pair through `C`**, and `fa_rep`
excludes 1-3 pairs — `corr(closure O···N, total_score) = +0.10`. Invisible to validation because
`validate.py` tests motif RMSD, seed phi and receptor distance, none of which involve the oxygen.

### 2. Cis and non-planar peptide bonds flanking the closure

| ring-distance from closure | n | non-planar >30° | cis (non-Pro) |
| --- | --- | --- | --- |
| 0 (the closure bond) | 239 | 0% | 0% |
| **1 (flanking)** | 478 | **47%** | **50%** |
| 2 | 478 | 24% | 43% |
| 3 | 478 | 1% | 0% |
| 4+ | 1193 | 0% | 0% |

Cause: GenKIC with CA pivots solves phi/psi at the pivots and leaves **omega free** — the note
above `variants.TARGET_LOOP` states this outright — and `PeptideCyclizeMover` constrains only the
closure bond, which is why row 0 is pristine and row 1 is not. The closure was bought by
distorting its neighbours.

Non-proline cis runs ~0.03% of bonds in real proteins; here it is 17%.

### 3. uOR only — the motif is D-amino acids

Residues **102–107** of the nanobody in `../uor_complex.pdb` have inverted backbone chirality
(improper N–CA–C–CB −50° to −92°, against +57.6 ± 1.3 for every other residue; 0 inverted in the
receptor, the apj nanobody, or experimental 9v32). Bond lengths are normal, so it is a
mirror-imaged backbone — a local Coot artifact; the source filename contains "flip" twice.

They are labelled with standard **L** three-letter codes, so Rosetta assigns L residue types on
read and rebuilds CB from ideal L geometry. Designs are therefore built L and superimposed back
onto D coordinates: CB lands 2.0–2.7 Å off and the sidechain points into the receptor. That is the
68-of-68 interface clash.

This propagated into all three uOR scaffolds and into the `flexddg/uor_nb/` alanine scan that
selected the motif in the first place.

## What is NOT invalidated

- **The flex ddG alanine scans.** They operate on the parent complexes, not the macrocycles.
  TROP2, CLDN18.2, CLDN6 and `apj_nb` all stand. `uor_nb` is compromised by defect 3 and must be
  re-run on the relaxed model (see the standing decision below).
- **The method.** The `i + j = 209` register analysis, route logic, seed-legality rules, motif
  pinning, the validation harness, and the scoring work on `rama_prepro` and interface terms are
  all sound and are reused by the rebuild.

## The gates that would have caught it

Added to `validate.py` on 2026-10-30 (`c369694`) and regression-tested against this archive:
**they reject 314 of 314 with zero false passes.** By gate — cis 311, intra 288, omega 270,
amide closure O 143, S-S distance 48, S-S chi3 14. Full table in `../cyclize/README.md` under
`validate`.

That regression is the point of keeping the archive reachable: any future change to the gates can
be re-checked against a set that is known-bad in five different ways.

## What replaced it

See `../STATUS.md` → "The 2026-10-30 rebuild". In short: relax both input complexes, re-excise the
motifs from the relaxed poses with identical residue numbering, realign designs to those same
relaxed models, and pin omega trans in GenKIC (`protocol.omega_trans`, `build.py`). A small
geometry test gates any production run.

**Standing decision (2026-10-30): the relaxed structures are the design basis.** The chirality
repair is therefore not a prerequisite — relaxing performs the L conversion Rosetta was already
doing silently, and taking the motif and the alignment target from the same relaxed pose removes
the mismatch. The accepted cost is that the uOR motif conformation is Rosetta's rather than the
experimental one (sidechain tips move up to 6.70 A). **The one uOR item this does not settle** is
that `../flexddg/uor_nb/`, the scan that selected the motif, ran on the D-residue model and must
be re-run on the relaxed one before its hotspots are trusted.
