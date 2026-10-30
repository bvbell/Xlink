# CLDN6 / 9v32-Fab chimera prep

**Date:** 2026-10-27
**Output:** `/scratch/users/bvbell/xlink/cldn6_9v32fab_complex.pdb`
**Scripts:** `/scratch/users/bvbell/xlink/scripts/0{1,2,3,4}_*.py`
**Intermediates:** `/scratch/users/bvbell/xlink/work/`

## What was built

CLDN6 (`cldn6.pdb`, chain A, 187 res) was structurally superimposed onto the
CLDN18.2 antigen of 9v32 (`9v32_complex.pdb` chain A, 171 res modeled, numbered
1-195) with **PyMOL 2.5.3 `super`** (whole-domain, all-atom, default settings;
module `viz pymol/2.5.3`). The resulting 4x4 transform was applied to every atom
of `cldn6.pdb`; 9v32 chains H and L were copied **byte-for-byte** (verified:
RMSD to original = 0.000000 A for both). Transform verified against PyMOL's own
output at 0.0005 A (PDB rounding only); det(R) = 1.0000 (proper rotation).

Chain order A (CLDN6, native numbering 1-187), H (2-219), L (2-216); TER between
chains; END at the end.

## Bottom line

The **superposition is good** (1.6 A over 763 atoms, 38% seq id, all three
methods agree to within 0.6 A). The **chimera as built is NOT directly usable
for flex ddG** without repacking/minimization: there are **56 heavy-atom pairs
below 2.0 A**, including several below 1.0 A (worst 0.45 A), concentrated in the
CLDN6 beta1-beta2 loop (res 39-44) against light-chain CDR-L1/CDR-L3. This is
not a method artifact -- all four transforms give 52-57 sub-2.0-A pairs.

---

## (a) Superposition quality

| method | RMSD (refined) | atom pairs | cycles | RMSD (pre-refine) | pre-refine pairs |
|---|---|---|---|---|---|
| `super` (all-atom, default) | **1.598 A** | 763 | 5 | 3.320 A | 954 |
| `super` (CA only) | 1.388 A | 124 CA | 5 | 2.925 A | 157 CA |
| `align` (sequence-driven) | 1.600 A | 815 | 5 | 4.315 A | 1034 |
| `cealign` | 3.365 A | 168 positions (no refinement) | - | - | - |

`cealign` reports no outlier-rejection, so its 3.37 A is directly comparable to
`super`'s **pre-refinement** 3.32 A, not to the 1.60 A.

**Method sensitivity (CA RMSD between the resulting CLDN6 placements):**

| | |
|---|---|
| super vs align | 0.248 A |
| super vs cealign | 0.619 A |
| align vs cealign | 0.524 A |

All three methods land CLDN6 in essentially the same place. The answer is **not**
method-sensitive.

**Sequence identity:** 139 residue pairs in the refined `super` alignment,
53 identical = **38.1%**.

## (b) Steric clashes, transformed CLDN6 vs Fab (H+L), heavy atoms

| cutoff | chimera | original 9v32 (reference) |
|---|---|---|
| < 2.0 A | **56** | 0 |
| < 2.5 A | **109** | 1 |
| < 3.0 A | **173** | 3 |

19 distinct CLDN6-Fab residue pairs come within 3.0 A.

Worst 15 residue pairs (minimum heavy-atom distance):

```
    CLDN6              Fab                 dist(A)   atoms
    A_GLN44            L_ASN31                0.45   NE2 -- CG
    A_ILE40            L_TYR98                0.68   CG2 -- OH
    A_VAL41            L_TYR98                0.72   O   -- CB
    A_VAL42            L_TYR98                0.79   N   -- CD1
    A_VAL41            L_LEU29                0.92   CG2 -- CA
    A_THR59            H_ASN103               1.09   CA  -- OD1
    A_VAL42            L_TYR38                1.40   CG2 -- CD2
    A_GLN57            L_TYR38                1.41   CD  -- OH
    A_ALA43            L_TYR98                1.87   N   -- O
    A_SER58            H_ASN103               2.04   C   -- OD1
    A_GLN61            H_ARG101               2.15   NE2 -- C
    A_GLN61            H_GLY102               2.17   NE2 -- N
    A_GLN63            L_TRP56                2.25   NE2 -- CH2
    A_THR59            H_GLY102               2.28   CG2 -- C
    A_SER39            L_SER28                2.34   O   -- O
```

Backbone-backbone overlaps (V41 O / V42 N / A43 N / T59 CA / S58 C) cannot be
fixed by repacking alone -- backbone minimization or loop remodeling is needed.

Method sensitivity of the clashes (same transform family, different method):

```
method       <2.0     <2.5     <3.0   minDist   CLDN6 iface res @5A
super          56      109      173      0.45        18
super_CA       57      101      164      0.50        18
align          55      105      175      0.36        18
cealign        52       95      147      0.41        18
```

## (c) Interface comparison, 5.0 A heavy-atom cutoff

**Antigen side**

ORIGINAL -- CLDN18.2 (9v32 chain A) residues contacting the Fab [15]:
```
P39 V40 T41 A42 V43 N45 V54 E56 S57 S58 G59 F60 E62 T167 Y169
```

CHIMERA -- CLDN6 residues contacting the Fab [18]:
```
F35 G37 N38 S39 I40 V41 V42 A43 Q44 V46 V55 Q57 S58 T59 G60 Q61 Q63 Q156
```

**Fab side**

ORIGINAL -- Fab residues contacting CLDN18.2 [18]:
```
H_W33 H_S99 H_W100 H_R101 H_G102 H_N103
L_N31 L_S32 L_G33 L_N34 L_K36 L_Y38 L_W56 L_D97 L_Y98 L_S99 L_Y100 L_F102
```

CHIMERA -- Fab residues contacting CLDN6 [26]:
```
H_S99 H_W100 H_R101 H_G102 H_N103 H_S104
L_Q27 L_S28 L_L29 L_L30 L_N31 L_S32 L_G33 L_N34 L_K36 L_N37 L_Y38 L_Y55
L_W56 L_G74 L_N96 L_D97 L_Y98 L_S99 L_Y100 L_F102
```

SHARED [17]:
```
H_S99 H_W100 H_R101 H_G102 H_N103
L_N31 L_S32 L_G33 L_N34 L_K36 L_Y38 L_W56 L_D97 L_Y98 L_S99 L_Y100 L_F102
```

ONLY in ORIGINAL (lost in the chimera) [1]: `H_W33`

ONLY in CHIMERA (new) [9]:
```
H_S104 L_Q27 L_S28 L_L29 L_L30 L_N37 L_Y55 L_G74 L_N96
```

**Epitope mapped through the structural alignment**

```
    CLDN18.2     CLDN6        contacts Fab in chimera?
    P39          (unaligned)  -
    V40          (unaligned)  -
    T41          (unaligned)  -
    A42          (unaligned)  -
    V43          (unaligned)  -
    N45          V46          yes
    V54          V55          yes
    E56          Q57          yes
    S57          S58          yes
    S58          (unaligned)  -
    G59          G60          yes
    F60          Q61          yes
    E62          Q63          yes
    T167         Q156         yes
    Y169         R158         no
```
CLDN6 contact residues with no counterpart in the CLDN18.2 epitope [10]:
`F35 G37 N38 S39 I40 V41 V42 A43 Q44 T59`

(The "unaligned" CLDN18.2 39-43 and 58 were dropped by `super`'s outlier
rejection, not genuinely unpaired; CLDN18.2 39-43 corresponds to CLDN6 ~40-44,
which is exactly the clashing loop.)

## (d) Extracellular loop superposition

Topology verified by DSS on each structure:

- CLDN18.2 (9v32 A): TM1 H3-26; ECL1 beta-strands 30-32, 44-47, 51-55, 62-65 +
  ECH; TM2 ~81-101; TM3 ~117-137; ECL2 ~137-168; TM4 ~169-195.
  **Gaps: 105-111 and 149-165** -- the ECL2 tip is unmodeled in 9v32.
- CLDN6: TM1 H4-26; ECL1 strands 30-35, 44-48, 52-57, 60-66; TM2 H75-101;
  TM3 H111-~136; ECL2 ~137-157 (short helix 154-156); TM4 H162-184. **No gaps.**

Via the alignment, CLDN18.2 ECL1 27-80 maps to CLDN6 27-81, and CLDN18.2
ECL2 137-168 maps to CLDN6 136-157.

Per-region CA RMSD on the `super` alignment pairs:

```
    region                                           nCA     RMSD   maxdev
    TM1        (18.2  3- 26 / CLDN6  3- 26)           18     1.08     1.67
    ECL1       (18.2 27- 80 / CLDN6 27- 81)           41     1.94     3.85
      ECL1 beta-sheet  (18.2 29- 66)                  26     1.48     2.68
      ECL1 helix ECH   (18.2 67- 80)                  13     2.72     3.85
    TM2        (18.2 81-101)                          21     1.76     5.34
    ICL/TM3    (18.2 102-136)                         24     1.81     3.55
    ECL2       (18.2 137-168, 149-165 unmodeled)       8     3.29     4.21
    TM4        (18.2 169-195)                         27     1.36     3.81
    ALL aligned CA                                   139     1.81     5.34
```

**ECL1 beta-sheet superimposes well (1.48 A)** -- this is the core of the
epitope. **ECL1 ECH (2.72 A) and ECL2 (3.29 A, only 8 aligned CA) diverge.**
ECL2 is the weakest region and is also the one where 9v32 is missing 149-165,
so the ECL2 part of the epitope is poorly constrained in both structures.

Per-residue CA deviation at the epitope:

```
    CLDN18.2     CLDN6          dev(A)
    N45          V46              1.59
    V54          V55              0.70
    E56          Q57              2.04
    S57          S58              2.53
    G59          G60              2.06
    F60          Q61              2.29
    E62          Q63              0.37
    T167         Q156             2.09
    Y169         R158             1.85
```

## (e) Parser hazards in the chimera

| check | result |
|---|---|
| altlocs | 0 |
| insertion codes | 0 |
| HETATM records | 0 |
| nonstandard residue names | none |
| atoms without element column | 0 |
| hydrogens | 0 |
| residues missing backbone atoms (N/CA/C/O) | 0 |
| residues missing side-chain heavy atoms | 0 |

Clean. Two cosmetic notes:

- 19 lines are 80 chars because they carry a formal-charge field in columns
  79-80 (`N1+`, `O1-` on Lys NZ / Glu OE2 etc.), inherited from the inputs.
  Standard PDB, harmless to column parsers.
- The segID column (73-76) is `A` for CLDN6 but `B`/`C` for chains H/L (copied
  verbatim from 9v32). Chain IDs in column 22 are correct (A/H/L). Only matters
  if something is run with Rosetta's `-in:use_segid`.

## (f) Backbone chain breaks (C(i)-N(i+1) > 2.0 A)

```
  chain A : residues 1-187 (187 modeled); 0 break(s)
  chain H : residues 2-219 (207 modeled); 2 break(s)
      131 -> 139 : C-N = 13.75 A  (numbering gap of 7)
      188 -> 193 : C-N = 12.50 A  (numbering gap of 4)
  chain L : residues 2-216 (215 modeled); 0 break(s)
```

Both chain-H breaks are pre-existing disordered loops in 9v32 (constant domain),
far from the interface. CLDN6 chain A is now gap-free, which is an improvement
over 9v32 chain A (which had gaps at 105-111 and 149-165).

## Recommendations before running flex ddG

1. Relieve the clashes. The sub-1 A backbone overlaps in CLDN6 39-44 vs CDR-L1 /
   CDR-L3 need backbone relaxation (constrained FastRelax with coordinate
   constraints, or Rosetta's `-relax:constrain_relax_to_start_coords`), not just
   side-chain repacking.
2. Re-run the interface analysis after relaxation; the 5.0 A contact sets above
   will shift, and the alanine-scan position list should be taken from the
   relaxed structure.
3. For a matched comparison against the CLDN18.2 run, restrict the alanine scan
   to the **17 shared Fab positions**; `H_W33` is present only in the original
   and the 9 chimera-only positions only in the chimera.
4. Treat any ECL2-side result with caution: 9v32 is missing CLDN18.2 149-165,
   the ECL2 region superimposes at 3.3 A, and CLDN18.2 Y169's counterpart
   (CLDN6 R158) drops out of the interface entirely.

## Reproduce

```bash
cd /scratch/users/bvbell/xlink
ml viz pymol/2.5.3
pymol -cq scripts/01_superpose.py                       # 3 methods + transforms
/home/groups/possu/miniconda/bin/python3 scripts/02_build_chimera.py
/home/groups/possu/miniconda/bin/python3 scripts/03_diagnostics.py
/home/groups/possu/miniconda/bin/python3 scripts/04_method_sensitivity.py
```
