#!/usr/bin/env python3
"""
Diagnostics for the CLDN6 / 9v32-Fab chimera.

    /home/groups/possu/miniconda/bin/python3 scripts/03_diagnostics.py

Reports: clashes, interface comparison (original vs chimera), per-region
superposition RMSD, parser hazards, backbone chain breaks.
"""
import os
import numpy as np
from scipy.spatial import cKDTree

BASE = "/scratch/users/bvbell/xlink"
WORK = os.path.join(BASE, "work")
ORIG = os.path.join(BASE, "9v32_complex.pdb")
CHIM = os.path.join(BASE, "cldn6_9v32fab_complex.pdb")

STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
          "THR TRP TYR VAL".split())

# expected heavy-atom side-chain counts (no OXT), for completeness checking
SC_ATOMS = {
    'ALA': ['CB'], 'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'ASN': ['CB', 'CG', 'OD1', 'ND2'], 'ASP': ['CB', 'CG', 'OD1', 'OD2'],
    'CYS': ['CB', 'SG'], 'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
    'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'], 'GLY': [],
    'HIS': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE': ['CB', 'CG1', 'CG2', 'CD1'], 'LEU': ['CB', 'CG', 'CD1', 'CD2'],
    'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'], 'MET': ['CB', 'CG', 'SD', 'CE'],
    'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'PRO': ['CB', 'CG', 'CD'],
    'SER': ['CB', 'OG'], 'THR': ['CB', 'OG1', 'CG2'],
    'TRP': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
    'VAL': ['CB', 'CG1', 'CG2'],
}

out = []


def log(s=""):
    print(s)
    out.append(s)


class Atom(object):
    __slots__ = ("rec", "serial", "name", "altloc", "resn", "chain", "resi",
                 "icode", "xyz", "elem")


def read_pdb(path):
    atoms = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        a = Atom()
        a.rec = line[0:6].strip()
        a.serial = line[6:11].strip()
        a.name = line[12:16].strip()
        a.altloc = line[16]
        a.resn = line[17:20].strip()
        a.chain = line[21]
        a.resi = int(line[22:26])
        a.icode = line[26]
        a.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        a.elem = line[76:78].strip() if len(line) > 77 else ""
        atoms.append(a)
    return atoms


def heavy(atoms):
    return [a for a in atoms
            if (a.elem.upper() != "H" if a.elem else not a.name.startswith("H"))
            and a.elem.upper() != "D"]


def tag(a):
    return "%s_%s%d" % (a.chain, a.resn, a.resi)


def pair_contacts(setA, setB, cutoff):
    """Return list of (dist, atomA, atomB) for all heavy-atom pairs < cutoff."""
    ca = np.array([a.xyz for a in setA])
    cb = np.array([b.xyz for b in setB])
    ta, tb = cKDTree(ca), cKDTree(cb)
    hits = ta.query_ball_tree(tb, cutoff)
    res = []
    for i, js in enumerate(hits):
        for j in js:
            d = float(np.linalg.norm(ca[i] - cb[j]))
            res.append((d, setA[i], setB[j]))
    return res


orig = read_pdb(ORIG)
chim = read_pdb(CHIM)

oA = heavy([a for a in orig if a.chain == "A"])
oF = heavy([a for a in orig if a.chain in "HL"])
cA = heavy([a for a in chim if a.chain == "A"])
cF = heavy([a for a in chim if a.chain in "HL"])

log("=" * 74)
log("CHIMERA DIAGNOSTICS  (CLDN6 + 9v32 Fab)")
log("=" * 74)
log("original 9v32: chain A %d heavy atoms, Fab H+L %d" % (len(oA), len(oF)))
log("chimera      : chain A (CLDN6) %d heavy atoms, Fab H+L %d" % (len(cA), len(cF)))

# ------------------------------------------------------------------ (b) clashes
log()
log("-" * 74)
log("(b) STERIC CLASHES: transformed CLDN6 (chain A) vs Fab (H+L)")
log("-" * 74)
cl = pair_contacts(cA, cF, 3.0)
for cut in (2.0, 2.5, 3.0):
    n = sum(1 for d, _, _ in cl if d < cut)
    log("  heavy-atom pairs < %.1f A : %d" % (cut, n))

# worst residue pairs, ranked by their minimum atom-atom distance
best = {}
for d, a, b in cl:
    k = (tag(a), tag(b))
    if k not in best or d < best[k][0]:
        best[k] = (d, a.name, b.name)
log()
log("  worst 15 clashing residue pairs (min heavy-atom distance):")
log("    %-18s %-18s %8s   %s" % ("CLDN6", "Fab", "dist(A)", "atoms"))
for (ka, kb), (d, na, nb) in sorted(best.items(), key=lambda kv: kv[1][0])[:15]:
    log("    %-18s %-18s %8.2f   %s -- %s" % (ka, kb, d, na, nb))
nres = len(best)
log()
log("  %d distinct CLDN6-Fab residue pairs come within 3.0 A." % nres)

# same measurement on the ORIGINAL complex, for reference
log()
log("  reference -- same counts for the ORIGINAL 9v32 CLDN18.2/Fab interface:")
clo = pair_contacts(oA, oF, 3.0)
for cut in (2.0, 2.5, 3.0):
    log("    heavy-atom pairs < %.1f A : %d"
        % (cut, sum(1 for d, _, _ in clo if d < cut)))

# ------------------------------------------------------- (c) interface comparison
log()
log("-" * 74)
log("(c) INTERFACE COMPARISON, 5.0 A heavy-atom cutoff")
log("-" * 74)


def iface(ag, fab):
    hits = pair_contacts(ag, fab, 5.0)
    ag_res, fab_res = {}, {}
    for d, a, b in hits:
        ag_res.setdefault(a.resi, a.resn)
        fab_res.setdefault((b.chain, b.resi), b.resn)
    return ag_res, fab_res


o_ag, o_fab = iface(oA, oF)
c_ag, c_fab = iface(cA, cF)

AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
          'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
          'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}


def fmt_ag(d):
    return " ".join("%s%d" % (AA3to1.get(d[r], 'X'), r) for r in sorted(d))


def fmt_fab(d):
    return " ".join("%s_%s%d" % (ch, AA3to1.get(d[(ch, r)], 'X'), r)
                    for ch, r in sorted(d, key=lambda k: (k[0], k[1])))


log()
log("ANTIGEN SIDE")
log("  ORIGINAL -- CLDN18.2 (9v32 chain A) residues contacting the Fab [%d]:"
    % len(o_ag))
log("    " + fmt_ag(o_ag))
log()
log("  CHIMERA  -- CLDN6 residues contacting the Fab [%d]:" % len(c_ag))
log("    " + fmt_ag(c_ag))

log()
log("FAB SIDE")
log("  ORIGINAL -- Fab residues contacting CLDN18.2 [%d]:" % len(o_fab))
log("    " + fmt_fab(o_fab))
log()
log("  CHIMERA  -- Fab residues contacting CLDN6 [%d]:" % len(c_fab))
log("    " + fmt_fab(c_fab))

shared = set(o_fab) & set(c_fab)
only_o = set(o_fab) - set(c_fab)
only_c = set(c_fab) - set(o_fab)
log()
log("  SHARED Fab positions (in both) [%d]:" % len(shared))
log("    " + fmt_fab({k: o_fab[k] for k in shared}))
log()
log("  ONLY in ORIGINAL (lost in the chimera) [%d]:" % len(only_o))
log("    " + (fmt_fab({k: o_fab[k] for k in only_o}) or "(none)"))
log()
log("  ONLY in CHIMERA (new) [%d]:" % len(only_c))
log("    " + (fmt_fab({k: c_fab[k] for k in only_c}) or "(none)"))

# antigen-side difference, expressed through the structural alignment
respairs = []
for line in open(os.path.join(WORK, "respairs_super.txt")):
    t_i, t_n, m_i, m_n = line.split()
    respairs.append((int(t_i), t_n, int(m_i), m_n))
t2m = dict((t_i, m_i) for t_i, _, m_i, _ in respairs)
m2t = dict((m_i, t_i) for t_i, _, m_i, _ in respairs)
log()
log("  ANTIGEN-side epitope mapped through the structural alignment")
log("  (CLDN18.2 epitope residue -> structurally equivalent CLDN6 residue;")
log("   'contact?' = does that CLDN6 residue contact the Fab in the chimera)")
log("    %-12s %-12s %-10s" % ("CLDN18.2", "CLDN6", "contact?"))
for r in sorted(o_ag):
    m = t2m.get(r)
    if m is None:
        log("    %-12s %-12s %-10s" % ("%s%d" % (AA3to1.get(o_ag[r], 'X'), r),
                                       "(unaligned)", "-"))
    else:
        mn = [p[3] for p in respairs if p[2] == m][0]
        log("    %-12s %-12s %-10s" % ("%s%d" % (AA3to1.get(o_ag[r], 'X'), r),
                                       "%s%d" % (AA3to1.get(mn, 'X'), m),
                                       "yes" if m in c_ag else "no"))
newc = [r for r in sorted(c_ag) if m2t.get(r) not in o_ag]
log("  CLDN6 contact residues with no counterpart in the CLDN18.2 epitope [%d]:"
    % len(newc))
log("    " + (" ".join("%s%d" % (AA3to1.get(c_ag[r], 'X'), r) for r in newc) or "(none)"))

# ------------------------------------------------------- (d) per-region RMSD
log()
log("-" * 74)
log("(d) PER-REGION SUPERPOSITION QUALITY (CA, on the `super` alignment pairs)")
log("-" * 74)
oca = dict(((a.resi), np.array(a.xyz)) for a in orig
           if a.chain == "A" and a.name == "CA")
cca = dict(((a.resi), np.array(a.xyz)) for a in chim
           if a.chain == "A" and a.name == "CA")

REGIONS = [
    ("TM1        (18.2  3- 26 / CLDN6  3- 26)", 3, 26),
    ("ECL1       (18.2 27- 80 / CLDN6 27- 81)", 27, 80),
    ("  ECL1 beta-sheet  (18.2 29- 66)", 29, 66),
    ("  ECL1 helix ECH   (18.2 67- 80)", 67, 80),
    ("TM2        (18.2 81-101)", 81, 101),
    ("ICL/TM3    (18.2 102-136)", 102, 136),
    ("ECL2       (18.2 137-168, 149-165 unmodeled)", 137, 168),
    ("TM4        (18.2 169-195)", 169, 195),
]
log("    %-46s %5s %8s %8s" % ("region", "nCA", "RMSD", "maxdev"))
allv = []
for name, lo, hi in REGIONS:
    dv = []
    for t_i, _, m_i, _ in respairs:
        if lo <= t_i <= hi and t_i in oca and m_i in cca:
            dv.append(np.linalg.norm(oca[t_i] - cca[m_i]))
    if not dv:
        log("    %-46s %5s %8s %8s" % (name, 0, "-", "-"))
        continue
    dv = np.array(dv)
    if not name.startswith("  "):
        allv.extend(dv)
    log("    %-46s %5d %8.2f %8.2f" % (name, len(dv), np.sqrt((dv ** 2).mean()),
                                       dv.max()))
allv = np.array(allv)
log("    %-46s %5d %8.2f %8.2f" % ("ALL aligned CA", len(allv),
                                   np.sqrt((allv ** 2).mean()), allv.max()))

log()
log("  per-residue CA deviation in the Fab-contacting epitope region:")
log("    %-12s %-12s %8s" % ("CLDN18.2", "CLDN6", "dev(A)"))
for r in sorted(o_ag):
    m = t2m.get(r)
    if m is None or r not in oca or m not in cca:
        log("    %-12s %-12s %8s" % ("%s%d" % (AA3to1.get(o_ag[r], 'X'), r),
                                     "(unaligned)", "-"))
        continue
    mn = [p[3] for p in respairs if p[2] == m][0]
    log("    %-12s %-12s %8.2f" % ("%s%d" % (AA3to1.get(o_ag[r], 'X'), r),
                                   "%s%d" % (AA3to1.get(mn, 'X'), m),
                                   np.linalg.norm(oca[r] - cca[m])))

# --------------------------------------------------- (e) parser hazards
log()
log("-" * 74)
log("(e) PARSER HAZARDS IN THE CHIMERA")
log("-" * 74)
alt = [a for a in chim if a.altloc != " "]
ico = [a for a in chim if a.icode != " "]
het = [a for a in chim if a.rec == "HETATM"]
nonstd = sorted(set(a.resn for a in chim if a.resn not in STD))
noelem = [a for a in chim if not a.elem]
log("  altlocs                : %d" % len(alt))
log("  insertion codes        : %d" % len(ico))
log("  HETATM records         : %d" % len(het))
log("  nonstandard residues   : %s" % (", ".join(nonstd) if nonstd else "none"))
log("  atoms w/o element col  : %d" % len(noelem))
log("  hydrogens              : %d" % (len(chim) - len(heavy(chim))))

byres = {}
for a in chim:
    byres.setdefault((a.chain, a.resi), []).append(a)
missing_bb, missing_sc = [], []
for (ch, ri), ats in sorted(byres.items()):
    names = set(a.name for a in ats)
    resn = ats[0].resn
    mbb = [n for n in ("N", "CA", "C", "O") if n not in names]
    if mbb:
        missing_bb.append("%s %s%d missing %s" % (ch, resn, ri, ",".join(mbb)))
    msc = [n for n in SC_ATOMS.get(resn, []) if n not in names]
    if msc:
        missing_sc.append("%s %s%d missing %s" % (ch, resn, ri, ",".join(msc)))
log("  residues missing backbone atoms : %d" % len(missing_bb))
for s in missing_bb:
    log("      " + s)
log("  residues missing side-chain atoms : %d" % len(missing_sc))
for s in missing_sc:
    log("      " + s)

# --------------------------------------------------- (f) chain breaks
log()
log("-" * 74)
log("(f) BACKBONE CHAIN BREAKS  (C(i)-N(i+1) > 2.0 A)")
log("-" * 74)
for ch in ("A", "H", "L"):
    res = sorted(set(a.resi for a in chim if a.chain == ch))
    C = dict((a.resi, np.array(a.xyz)) for a in chim
             if a.chain == ch and a.name == "C")
    N = dict((a.resi, np.array(a.xyz)) for a in chim
             if a.chain == ch and a.name == "N")
    breaks = []
    for i in range(len(res) - 1):
        r1, r2 = res[i], res[i + 1]
        if r1 in C and r2 in N:
            d = float(np.linalg.norm(C[r1] - N[r2]))
            if d > 2.0:
                breaks.append((r1, r2, d, r2 - r1))
    log("  chain %s : residues %d-%d (%d modeled); %d break(s)"
        % (ch, res[0], res[-1], len(res), len(breaks)))
    for r1, r2, d, gap in breaks:
        log("      %d -> %d : C-N = %.2f A  (%s)"
            % (r1, r2, d, "numbering gap of %d" % (gap - 1) if gap > 1
               else "consecutive numbering"))

with open(os.path.join(WORK, "diagnostics.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
print("\nwrote " + os.path.join(WORK, "diagnostics.txt"))
