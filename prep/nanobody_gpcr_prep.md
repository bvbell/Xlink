# Nanobody / GPCR complex prep for Rosetta flex ddG

Prepared 2026-10-27 on Sherlock compute node sh02-02n20 (job 49640222).
No Slurm jobs were submitted; nothing under `flexddg/` was touched.

## Files

| role | path |
|---|---|
| INPUT 1 source (AF3 mmCIF) | `/scratch/users/bvbell/xlink/d006_lmpnn_fkd2_run7_af3 (1).cif` |
| INPUT 1 output | `/scratch/users/bvbell/xlink/prep/d006_complex.pdb` |
| INPUT 2 source (Coot PDB) | `/scratch/users/bvbell/xlink/Rosetta_9_real_trunc_model_1_flip_fit_J126map-coot-0_flip.pdb` |
| INPUT 2 output | `/scratch/users/bvbell/xlink/prep/rosetta9_complex.pdb` |
| reader/writer library | `/scratch/users/bvbell/xlink/scripts/pdbio.py` |
| conversion driver | `/scratch/users/bvbell/xlink/scripts/prep_structures.py` |
| diagnostics | `/scratch/users/bvbell/xlink/scripts/diagnose_complex.py` |

Run with `/home/groups/possu/miniconda/bin/python3` (Biopython 1.83 / numpy 1.24 env).
A hand-written mmCIF `_atom_site` parser and a strict column PDB writer were used
rather than gemmi/Biopython, so that auth chain IDs and auth residue numbers pass
through untouched and the output columns are exactly the ones a strict parser wants
(22 chain, 23-26 resnum, 27 icode, 18-20 resname, 31-54 xyz, 55-60 occ, 61-66 B,
77-78 element, `TER` between chains, `END` at the end).

## Conversion summary

Both conversions were clean and required no repair:

* Single model in both sources; nothing dropped.
* No altlocs anywhere (every altloc field blank), so nothing dropped for altloc.
* No HETATM records, no insertion codes, no hydrogens, no waters/ligands.
* **No atom has occupancy 0.00 in either file.** Every atom is 1.00. Nothing had to be
  forced to 1.00, so there is nothing to warn about on that front.
* No side chains stripped, no atoms added, no minimisation, no repacking.
* Atom counts: INPUT 1 3717 atoms (A 929 + B 2788), INPUT 2 3264 atoms (A 963 + B 2301).
* Every record is padded to the 80-column PDB width; both files parse under
  `Bio.PDB.PDBParser(PERMISSIVE=0)` with `warnings.simplefilter("error")` and yield
  exactly one model with chains A and B at the expected residue/atom counts.
* For INPUT 2 the emitted `ATOM` records are **byte-for-byte identical** to the source
  over columns 1-80 (`diff` is empty), i.e. the Coot file was already conformant and the
  "normalisation" is a no-op apart from the regenerated `TER`/`END`.

## Sequence identification (not requested, but it anchors the topology calls)

* INPUT 1 chain B = human apelin receptor **APJ / APLNR**, modelled residues 1-351 where
  model residue *n* = APLNR residue *n+17* up to model 212, and *n-24* from model 267 on.
  Independent confirmation: the two native APJ disulfides are both present and both land
  exactly where they should (model C85-C164 = APLNR C102(3.25)-C181 ECL2; model
  C2-C305 = APLNR C19 N-term - C281 ECL3).
* INPUT 1 chain B residues **213-266** are a 54-residue **rubredoxin** replacing ICL3
  (`MKKYTCTVCGYIYNPEDGDPDNGVNPGTDFKDIPDDWVCPLCGVGKDQFEEVEE`, two CXXCG iron motifs).
* INPUT 2 chain B = human **mu-opioid receptor (OPRM1)**, model residue *n* = OPRM1
  residue *n+65* throughout. Confirmed by D99-R100-Y101 = DRY (D164/R165/Y166),
  N267-P268-x-x-Y271 = NPxxY (N332/Y336), and the C75-C152 disulfide = OPRM1
  C140(3.25)-C217(ECL2). The construct runs from the top of TM1 to the end of helix 8;
  both termini are truncated.
* Both chains A share the same VHH framework (`EVQLVESGGGLVQPGGSLRLSCAAS...WGQGTLVTVSS`)
  and differ only in the CDRs - a designed-library pair.

## INPUT 1 - d006_complex.pdb (AF3 model; chain A nanobody d006, chain B apelin receptor APJ with rubredoxin in ICL3)

`prep/d006_complex.pdb`

### 1. Per-chain composition

| chain | residues | heavy atoms | all atoms | first | last | numbering gaps |
|---|---|---|---|---|---|---|
| A | 123 | 929 | 929 | 1 | 123 | none |
| B | 351 | 2788 | 2788 | 1 | 351 | none |

### 2. A/B interface, 5.0 A heavy-atom cutoff

Chain A (nanobody) interface residues: **30**

| residue | min dist to B (A) |
|---|---|
| A_E1 | 3.27 |
| A_R27 | 2.89 |
| A_A29 | 3.05 |
| A_G30 | 4.10 |
| A_V32 | 3.66 |
| A_Y33 | 2.82 |
| A_A35 | 4.84 |
| A_Y37 | 3.08 |
| A_L47 | 3.88 |
| A_K53 | 3.85 |
| A_N73 | 2.87 |
| A_A74 | 4.06 |
| A_K75 | 4.81 |
| A_N76 | 2.69 |
| A_Q96 | 4.47 |
| A_V98 | 3.22 |
| A_K100 | 2.80 |
| A_V101 | 3.57 |
| A_L102 | 3.25 |
| A_G103 | 3.50 |
| A_F104 | 2.44 |
| A_I105 | 3.09 |
| A_F106 | 3.21 |
| A_P107 | 3.42 |
| A_E108 | 2.65 |
| A_L109 | 3.34 |
| A_K110 | 2.70 |
| A_P111 | 4.03 |
| A_Y112 | 4.41 |
| A_W113 | 3.27 |

Chain B (GPCR) interface residues: **47**

| residue | min dist to A (A) |
|---|---|
| B_E1 | 2.69 |
| B_C2 | 3.06 |
| B_E3 | 3.85 |
| B_Y4 | 3.05 |
| B_D6 | 2.80 |
| B_Y18 | 3.28 |
| B_W68 | 3.21 |
| B_Y71 | 4.35 |
| B_T72 | 3.91 |
| B_Y76 | 3.91 |
| B_S89 | 3.66 |
| B_I92 | 3.85 |
| B_F93 | 3.29 |
| B_R151 | 2.65 |
| B_L156 | 3.22 |
| B_E157 | 3.34 |
| B_N158 | 2.82 |
| B_Q163 | 3.67 |
| B_C164 | 4.31 |
| B_Y165 | 3.58 |
| B_M166 | 3.72 |
| B_Y168 | 2.98 |
| B_V174 | 3.27 |
| B_S175 | 4.39 |
| B_E177 | 3.04 |
| B_W178 | 3.27 |
| B_E181 | 2.70 |
| B_L184 | 3.56 |
| B_G185 | 3.69 |
| B_S188 | 3.89 |
| B_Y288 | 2.44 |
| B_V291 | 3.99 |
| B_K292 | 3.41 |
| B_Y295 | 2.93 |
| B_M296 | 4.93 |
| B_G298 | 2.89 |
| B_S299 | 3.22 |
| B_W303 | 4.27 |
| B_C305 | 3.26 |
| B_D308 | 2.99 |
| B_L309 | 4.67 |
| B_L311 | 3.61 |
| B_M312 | 3.42 |
| B_F315 | 3.25 |
| B_P316 | 3.69 |
| B_T319 | 3.09 |
| B_Y323 | 3.24 |

Compact lists:

- A: `A_E1, A_R27, A_A29, A_G30, A_V32, A_Y33, A_A35, A_Y37, A_L47, A_K53, A_N73, A_A74, A_K75, A_N76, A_Q96, A_V98, A_K100, A_V101, A_L102, A_G103, A_F104, A_I105, A_F106, A_P107, A_E108, A_L109, A_K110, A_P111, A_Y112, A_W113`
- B: `B_E1, B_C2, B_E3, B_Y4, B_D6, B_Y18, B_W68, B_Y71, B_T72, B_Y76, B_S89, B_I92, B_F93, B_R151, B_L156, B_E157, B_N158, B_Q163, B_C164, B_Y165, B_M166, B_Y168, B_V174, B_S175, B_E177, B_W178, B_E181, B_L184, B_G185, B_S188, B_Y288, B_V291, B_K292, B_Y295, B_M296, B_G298, B_S299, B_W303, B_C305, B_D308, B_L309, B_L311, B_M312, B_F315, B_P316, B_T319, B_Y323`

### 3. Which receptor face does the nanobody bind?

Membrane normal approximated as PC1 of chain-B CA coordinates excluding the fusion domain (PC1 captures 60.8% of the variance). Sign fixed so that chain-B residue 1 (N-terminal end of TM1, extracellular by topology) is at positive projection.

| quantity | projection on normal (A) |
|---|---|
| chain B N-term CA (res 1) | +32.1 |
| chain B C-term CA (res 351) | -18.9 |
| all chain B CA, mean | -0.0 |
| all chain B CA, range | -28.8 .. +32.1 |
| chain B *interface* CA, mean | +17.6 |
| chain B interface CA, range | +2.9 .. +32.1 |
| whole nanobody (chain A) CA, mean | +39.3 |
| fusion domain CA, mean | -37.4 |

Interface centroid sits +17.6 A relative to the receptor centroid along the normal -> **EXTRACELLULAR face**.

Per-interface-residue projections (positive = extracellular side):

`B_E1:+32, B_C2:+30, B_E3:+30, B_Y4:+27, B_D6:+23, B_Y18:+6, B_W68:+11, B_Y71:+15, B_T72:+15, B_Y76:+20, B_S89:+9, B_I92:+5, B_F93:+4, B_R151:+15, B_L156:+28, B_E157:+30, B_N158:+32, B_Q163:+21, B_C164:+18, B_Y165:+18, B_M166:+18, B_Y168:+20, B_V174:+27, B_S175:+26, B_E177:+21, B_W178:+21, B_E181:+16, B_L184:+11, B_G185:+11, B_S188:+6, B_Y288:+7, B_V291:+11, B_K292:+12, B_Y295:+16, B_M296:+15, B_G298:+19, B_S299:+20, B_W303:+21, B_C305:+25, B_D308:+21, B_L309:+22, B_L311:+16, B_M312:+18, B_F315:+13, B_P316:+13, B_T319:+8, B_Y323:+3`

### 4. Fusion domain

Residues B 213-266 (54 residues), sequence `MKKYTCTVCGYIYNPEDGDPDNGVNPGTDFKDIPDDWVCPLCGVGKDQFEEVEE`.

Minimum heavy-atom distance from the fusion domain to chain A: **39.24 A** -> NO contact with the nanobody at 5.0 A.

### 5. Inter-chain steric clashes (A-B heavy atoms < 2.0 A)

Total clashing heavy-atom pairs: **0**

### 6. Backbone chain breaks (C(i)-N(i+1) > 2.0 A)

None. Every consecutive residue pair in both chains is bonded.

### 7. Residue-level sanity

- HETATM records: 0
- altloc characters present in source: [' ']
- atoms dropped for altloc != blank/A: 0
- insertion codes: 0
- models in source: [1] (only model 1 retained)
- atoms with occupancy 0.00: 0
- hydrogens present: 0

Chain A:

- nonstandard residue names: none
- residues missing backbone atoms: none
- residues with incomplete side chains: none

Chain B:

- nonstandard residue names: none
- residues missing backbone atoms: none
- residues with incomplete side chains: none

### 8. Disulfides (SG-SG < 2.5 A)

- B_C85 - B_C164 : 2.08 A  (intra-chain B)
- A_C22 - A_C95 : 2.13 A  (intra-chain A)
- B_C2 - B_C305 : 2.23 A  (intra-chain B)

- free cysteines: B_C102, B_C203, B_C218, B_C221, B_C251, B_C254, B_C284, B_C318, B_C320, B_C327, B_C344


---

## INPUT 2 - rosetta9_complex.pdb (Coot-built experimental model; chain A nanobody, chain B mu-opioid receptor)

`prep/rosetta9_complex.pdb`

### 1. Per-chain composition

| chain | residues | heavy atoms | all atoms | first | last | numbering gaps |
|---|---|---|---|---|---|---|
| A | 123 | 963 | 963 | 1 | 123 | none |
| B | 286 | 2301 | 2301 | 1 | 286 | none |

### 2. A/B interface, 5.0 A heavy-atom cutoff

Chain A (nanobody) interface residues: **24**

| residue | min dist to B (A) |
|---|---|
| A_A29 | 3.11 |
| A_Q30 | 2.63 |
| A_Y31 | 4.26 |
| A_L32 | 2.73 |
| A_W33 | 3.79 |
| A_Y37 | 2.55 |
| A_L47 | 4.10 |
| A_Y52 | 4.13 |
| A_S53 | 2.70 |
| A_E54 | 2.91 |
| A_K56 | 3.20 |
| A_L96 | 4.29 |
| A_V98 | 2.67 |
| A_V100 | 2.88 |
| A_S101 | 4.16 |
| A_Y102 | 2.48 |
| A_N103 | 2.95 |
| A_F104 | 3.44 |
| A_F105 | 3.34 |
| A_F106 | 2.87 |
| A_R107 | 3.27 |
| A_S108 | 3.78 |
| A_Y109 | 3.78 |
| A_E111 | 2.27 |

Chain B (GPCR) interface residues: **33**

| residue | min dist to A (A) |
|---|---|
| B_Q59 | 4.07 |
| B_Y63 | 4.65 |
| B_D82 | 3.65 |
| B_Y83 | 4.01 |
| B_M86 | 4.78 |
| B_R146 | 2.27 |
| B_Q147 | 2.55 |
| B_D151 | 3.78 |
| B_C152 | 4.94 |
| B_T153 | 3.23 |
| B_L154 | 3.03 |
| B_F156 | 2.87 |
| B_E164 | 3.02 |
| B_L167 | 4.56 |
| B_K168 | 2.95 |
| B_V171 | 3.44 |
| B_F172 | 4.73 |
| B_I231 | 4.54 |
| B_H232 | 3.88 |
| B_V235 | 3.45 |
| B_K238 | 3.39 |
| B_T242 | 2.63 |
| B_I243 | 2.89 |
| B_P244 | 2.79 |
| B_E245 | 2.70 |
| B_T246 | 2.91 |
| B_T247 | 3.09 |
| B_Q249 | 2.48 |
| B_T250 | 3.78 |
| B_W253 | 3.03 |
| B_H254 | 3.27 |
| B_I257 | 4.00 |
| B_Y261 | 3.34 |

Compact lists:

- A: `A_A29, A_Q30, A_Y31, A_L32, A_W33, A_Y37, A_L47, A_Y52, A_S53, A_E54, A_K56, A_L96, A_V98, A_V100, A_S101, A_Y102, A_N103, A_F104, A_F105, A_F106, A_R107, A_S108, A_Y109, A_E111`
- B: `B_Q59, B_Y63, B_D82, B_Y83, B_M86, B_R146, B_Q147, B_D151, B_C152, B_T153, B_L154, B_F156, B_E164, B_L167, B_K168, B_V171, B_F172, B_I231, B_H232, B_V235, B_K238, B_T242, B_I243, B_P244, B_E245, B_T246, B_T247, B_Q249, B_T250, B_W253, B_H254, B_I257, B_Y261`

### 3. Which receptor face does the nanobody bind?

Membrane normal approximated as PC1 of chain-B CA coordinates (PC1 captures 64.1% of the variance). Sign fixed so that chain-B residue 1 (N-terminal end of TM1, extracellular by topology) is at positive projection.

| quantity | projection on normal (A) |
|---|---|
| chain B N-term CA (res 1) | +20.3 |
| chain B C-term CA (res 286) | -15.1 |
| all chain B CA, mean | -0.0 |
| all chain B CA, range | -39.0 .. +35.0 |
| chain B *interface* CA, mean | +16.7 |
| chain B interface CA, range | +4.1 .. +35.0 |
| whole nanobody (chain A) CA, mean | +41.9 |

Interface centroid sits +16.7 A relative to the receptor centroid along the normal -> **EXTRACELLULAR face**.

Per-interface-residue projections (positive = extracellular side):

`B_Q59:+16, B_Y63:+21, B_D82:+9, B_Y83:+8, B_M86:+4, B_R146:+32, B_Q147:+35, B_D151:+25, B_C152:+22, B_T153:+21, B_L154:+20, B_F156:+21, B_E164:+18, B_L167:+12, B_K168:+12, B_V171:+7, B_F172:+6, B_I231:+8, B_H232:+8, B_V235:+13, B_K238:+17, B_T242:+22, B_I243:+22, B_P244:+24, B_E245:+24, B_T246:+23, B_T247:+21, B_Q249:+18, B_T250:+19, B_W253:+14, B_H254:+14, B_I257:+10, B_Y261:+5`

### 5. Inter-chain steric clashes (A-B heavy atoms < 2.0 A)

Total clashing heavy-atom pairs: **0**

### 6. Backbone chain breaks (C(i)-N(i+1) > 2.0 A)

None. Every consecutive residue pair in both chains is bonded.

### 7. Residue-level sanity

- HETATM records: 0
- altloc characters present in source: [' ']
- atoms dropped for altloc != blank/A: 0
- insertion codes: 0
- models in source: [1] (only model 1 retained)
- atoms with occupancy 0.00: 0
- hydrogens present: 0

Chain A:

- nonstandard residue names: none
- residues missing backbone atoms: none
- residues with incomplete side chains: none

Chain B:

- nonstandard residue names: none
- residues missing backbone atoms: none
- residues with incomplete side chains: none

### 8. Disulfides (SG-SG < 2.5 A)

- B_C75 - B_C152 : 1.88 A  (intra-chain B)
- A_C22 - A_C95 : 2.11 A  (intra-chain A)

- free cysteines: B_C14, B_C94, B_C105, B_C125, B_C170, B_C186, B_C227, B_C256, B_C265, B_C281, B_C286


### 4. Fusion domain

Not applicable - INPUT 2 chain B is a plain (N- and C-terminally truncated) mu-opioid
receptor with no inserted domain. Its ICL3 (model ~196-210, OPRM1 ~261-275) is native
and short.

---

## Topology mapping of the two epitopes

Converting the interface chain-B residues into receptor numbering/topology
(INPUT 1 offsets +17 / -24, INPUT 2 offset +65):

**INPUT 1, APJ epitope**

| segment | model residues |
|---|---|
| N-terminus (extracellular) | E1, C2, E3, Y4, D6, Y18 |
| TM2 extracellular half | W68, Y71, T72, Y76 |
| TM3 top / orthosteric pocket (3.29-3.33) | S89, I92, F93 |
| ECL2 (incl. the C164 disulfide Cys) + ECL2/TM5 junction | R151, L156, E157, N158, Q163, C164, Y165, M166, Y168, V174, S175, E177, W178, E181, L184, G185, S188 |
| TM6 extracellular half | Y288, V291, K292, Y295, M296, G298, S299 |
| ECL3 (incl. the C305 disulfide Cys) | W303, C305, D308, L309 |
| TM7 extracellular half, down to ~7.43 | L311, M312, F315, P316, T319, Y323 |

**INPUT 2, mu-opioid epitope**

| segment | model residues (OPRM1 numbering) |
|---|---|
| TM2 extracellular half / ECL1 | Q59 (Q124), Y63 (Y128) |
| TM3 orthosteric pocket (3.32-3.36) | D82 (**D147, the orthosteric Asp3.32**), Y83 (Y148), M86 (M151) |
| ECL2 (incl. the C152 disulfide Cys) | R146 (R211), Q147 (Q212), D151 (D216), C152 (C217), T153 (T218), L154 (L219), F156 (F221) |
| TM5 extracellular half | E164 (E229), L167 (L232), K168 (K233), V171 (V236), F172 (F237) |
| TM6 extracellular half | I231 (I296), H232 (**H297, 6.52**), V235 (V300), K238 (K303) |
| ECL3 | T242-T250 (T307-T315) |
| TM7 extracellular half | W253 (W318), H254 (H319), I257 (I322), Y261 (**Y326, 7.43**) |

In **neither** complex does a single ICL1, ICL2, ICL3, helix-8 or C-terminal residue
appear at the interface. Both epitopes are purely extracellular, and in both cases the
long CDR3 reaches down past the extracellular vestibule into the orthosteric pocket
(the deepest interface residues sit only ~3-5 A above the membrane mid-plane:
APJ F93/Y323, MOR M86/Y261 = M151/Y326).

## Additional checks beyond the requested list

* **Intra-chain** heavy-atom contacts < 2.0 A (|i-j| > 1 in sequence): **zero** in INPUT 1
  (both chains) and **one** in INPUT 2 (the C75-C152 SG-SG disulfide at 1.88 A, i.e. not a
  real clash but a short bond, see below). So neither model has internal packing failures.
* **Peptide-bond geometry**: C(i)-N(i+1) is 1.31-1.36 A in INPUT 1 and 1.29-1.35 A in
  INPUT 2. Nothing even approaches the 2.0 A break threshold; both chains are fully
  continuous.
* **Side-chain completeness**: expected-vs-observed heavy-atom counts match exactly for
  all four chains (INPUT 1 A 928+1 OXT / B 2787+1 OXT; INPUT 2 A 963 / B 2301, no OXT).
  Despite the expectation for a Coot model, **no residue in either structure is missing a
  single heavy atom**.

## Problems, bluntly

1. **INPUT 1 chain B carries a buried Lys where class A GPCRs have the invariant
   W6.48 toggle switch.** The TM6 motif reads `...T280 F281 A282 L283 C284 K285 M286
   P287...`, i.e. `C-K-M-P` where APJ (and essentially every class A receptor) has the
   `CWxP` motif with Trp6.48. Lys285 is buried in the TM6/TM7 core (16 residues within
   5 A: N95, S99, F193, F281, L283, C284, M286, P287, Y288, H289, S322, Y323, N325,
   S326) with no acidic partner. This is almost certainly a LigandMPNN redesign artifact
   (`lmpnn` is in the filename) and it is a genuine model-quality defect in the receptor
   core. It is 10.0 A from the nanobody and its neighbour Y323 *is* an interface residue,
   so it will not dominate an interface ddG but it does sit one shell away from the
   epitope. Check whether your actual construct has this mutation. The model's C-terminal
   `...CTSML-L348-L349-M350-G351` also departs from the APJ palmitoylation `CC` motif;
   that one is 32-39 A from the interface and harmless.
2. **INPUT 1's rubredoxin has four cysteines converging on an empty iron site.** SG-SG
   distances C218-C221 3.52, C218-C254 3.59, C221-C251 3.47, C221-C254 3.15,
   C251-C254 3.53, C218-C251 4.09 A. There is **no Fe (or any HETATM) in the file**, so
   Rosetta will score four buried, desolvated, unsatisfied thiolates. All are >39 A from
   the nanobody, so this is pure additive noise in the total score rather than an
   interface problem - but if you ever compare absolute scores across structures, be
   aware of it. If flex ddG gets unhappy, deleting residues 213-266 is safe: the fusion's
   minimum heavy-atom distance to chain A is 39.24 A.
3. **INPUT 2's C75-C152 disulfide is over-compressed at 1.88 A** (a normal SS bond is
   2.03-2.05 A). Coot/real-space refinement squeezed it. Rosetta will detect it as a
   disulfide and relax it, but expect a nonzero starting strain term. By contrast
   INPUT 1's C2-C305 is **long at 2.23 A** and C85-C164 is 2.08 A; A_C22-A_C95 is
   2.13 A (INPUT 1) / 2.11 A (INPUT 2).
4. **No HETATM/ligand in either file.** INPUT 2 is an opioid receptor model with an empty
   orthosteric pocket that the nanobody CDR3 partly occupies - fine for A/B ddG, but if
   the real experiment has an agonist bound, this model does not represent it.
5. **No hydrogens and no OXT in INPUT 2.** Rosetta will add both; harmless.
6. Nothing else. No gaps in numbering, no chain breaks, no altlocs, no insertion codes,
   no nonstandard residues, no zero occupancies, no inter-chain clashes in either
   structure. These two files are unusually clean inputs.
