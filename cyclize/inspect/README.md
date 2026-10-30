# Manual inspection set (built 2026-10-30)

Regenerate with `python3 scripts/inspect_set.py --n 6`. Selection criteria and metric definitions
are documented in that script's docstring; the numbers are in `manifest.csv`.

```
inspect/
  apj/{good,bad}/     6 each, from apj + apj_i1
  uor/{good,bad}/     6 each, from uor + uor_i1
  reference/          native complexes + the scaffolds the motifs were excised from
  manifest.csv        every metric for every copied structure
```

Chain `P` is the design, chain `B` the receptor.

## Load with `view.pml` — otherwise every ring looks broken

```bash
cd /scratch/users/bvbell/xlink/cyclize/inspect && pymol view.pml
```

**Rosetta writes no `LINK` or `CONECT` record for the cyclisation**, and viewers only auto-bond
sequentially numbered residues. The ring-closing bond from the last residue's `C` to residue 1's
`N` is therefore never drawn, and a perfectly good macrocycle renders as an **open chain whose two
ends happen to sit 1.33 Å apart**. That is a display artifact, not a defect.

The closure geometry is in fact textbook in every amide design in this set:

| | measured | ideal |
| --- | --- | --- |
| closure `C–N` | **1.33 Å** | 1.33 |
| `CA–C–N` | **116.2°** | 116.2 |
| `C–N–CA` | **121.7°** | 121.7 |

`view.pml` issues the explicit `bond` commands. The copied PDBs also carry `CONECT` records for
the ring bond, which ChimeraX and most other viewers honour, but PyMOL ignores `CONECT` for `ATOM`
records under its default `connect_mode`, which is why the script exists. The two `ssapp`
structures close through S–S instead and PyMOL draws disulfides on its own.

**Judge the closure only after loading this way.** The real defects are the ones below.

## "Good" is relative. Nothing here is clean.

Ranking is on geometry, **not** on `dG_separated` — that metric cancels intra-peptide energy and
therefore rewards self-intersection, which is how the worst structures came to top the score
ranking in the first place. See `../../STATUS.md`.

## What to look at, in order

**1. Omega planarity — the primary defect, both targets.** Measure the CA–C–N–CA dihedral around
each ring. It should be ~180° (trans) or ~0° (cis, and legitimate only before proline).

- apj "good" bottoms out at **25°**, uOR "good" at ~10° — **and those are exactly the input
  maxima** (apj nanobody 25.2°, uOR 10.3°). They are inherited pre-proline bonds from the
  scaffold, they are normal (experimental 9v32 reaches 21.1°), and a design sitting at that value
  has introduced **no** non-planarity of its own. `apj/good/apj__medium__amide__DPRO-GLY__c00.pdb`
  is the clearest case.
- All "bad" structures sit at **~90°**, a fully broken amide with no partial double-bond character.
  No input bond anywhere exceeds 30°, so everything above that threshold was introduced by the
  build — and it lands at the **turn and the ring-closure seam**, not in the motif interior. That
  localisation is the thing to look for.

In PyMOL: `set_dihedral` readouts, or just look down the backbone for a carbonyl twisted out of
the plane of its neighbours.

**2. Self-intersection.** `intra_min` is the closest non-adjacent intra-peptide heavy-atom contact.
Normal packing is ≥3.4 Å. The worst here is **0.38 Å**; the best in the entire apj set is 2.99 Å.
Look for the chain passing through its own ring — in `apj_i1 long__amide__DALA-GLY__c01` two
non-adjacent alpha carbons are 1.87 Å apart.

**3. uOR only — the D-amino-acid motif.** Residues 102–107 of the nanobody in
`reference/uor_complex.pdb` are mirror-imaged (improper N–CA–C–CB negative). Superimpose any uOR
design onto the native motif and watch CB: backbone matches to 0.1–0.5 Å while CB is 2.0–2.7 Å
away. CB is not a rotamer degree of freedom, so that displacement can only be chirality.

> **Caveat on the `n_dres` column.** It counts *all* backbone D residues, and the i+1 sweeps
> (`apj_i1`, `uor_i1`) deliberately enumerate D residues at the turn — `DTHR`, `DPRO`, `DTRP` and
> so on. So a non-zero `n_dres` is only evidence of a defect for uOR, where it reflects the
> inherited motif, and only above the count expected from the variant's own name. apj's values
> of 1 are the designed turn residue and are correct.

**4. Interface clash.** `inter_min` is the closest design-to-receptor contact; <2.0 Å is a clash.
Every uOR structure here is clashing (1.0–1.5 Å) — that is the chirality defect expressing itself
at the interface. apj's "good" rows reach 2.3–2.5 Å, which is normal contact distance.

## The one structure worth most attention

`apj/good/apj__medium__amide__DPRO-GLY__c00.pdb` is the only copied structure with **zero
non-proline cis bonds**, one legitimate cis-Pro, `intra_min` 2.68 Å, `inter_min` 2.45 Å and motif
RMSD 0.17 Å. Its worst bond is 25° off planar, but that is the inherited pre-Pro bond from the scaffold — it
introduces no non-planarity of its own, which makes it the existence proof that a 2-residue turn
*can* close with planar amides. Use it as the baseline for judging whether a repaired pipeline is
actually producing better structures.

Note it is **not** the design that topped the old `dG_separated` ranking — those are all in `bad/`.
