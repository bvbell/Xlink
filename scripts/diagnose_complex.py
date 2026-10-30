#!/usr/bin/env python3
"""Diagnostics for a two-chain nanobody(A)/GPCR(B) complex PDB.

Reports: per-chain composition and numbering gaps, A/B interface at a
heavy-atom cutoff, membrane-normal analysis of which receptor face is bound,
steric clashes, backbone chain breaks, residue-level completeness, disulfides.
"""
import os
import sys
import argparse
from collections import OrderedDict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdbio import read_pdb, AA3, RES_HEAVY, BACKBONE  # noqa: E402


class Residue(object):
    def __init__(self, chain, num, icode, name):
        self.chain, self.num, self.icode, self.name = chain, num, icode, name
        self.atoms = []

    @property
    def one(self):
        return AA3.get(self.name, "X")

    @property
    def tag(self):
        return "%s_%s%d%s" % (self.chain, self.one, self.num,
                              self.icode.strip())

    def coord(self, aname):
        for a in self.atoms:
            if a.name == aname:
                return np.array([a.x, a.y, a.z])
        return None

    def heavy(self):
        return [a for a in self.atoms if not a.is_hydrogen]


def build(atoms):
    chains = OrderedDict()
    for a in atoms:
        ch = chains.setdefault(a.chain, OrderedDict())
        key = (a.resnum, a.icode)
        if key not in ch:
            ch[key] = Residue(a.chain, a.resnum, a.icode, a.resname)
        ch[key].atoms.append(a)
    return OrderedDict((c, list(v.values())) for c, v in chains.items())


def arr(residues):
    """Return (coords Nx3, residue-index per atom) over heavy atoms."""
    xyz, idx = [], []
    for i, r in enumerate(residues):
        for a in r.heavy():
            xyz.append((a.x, a.y, a.z))
            idx.append(i)
    return np.asarray(xyz, float), np.asarray(idx, int)


def pair_min(resA, resB, cutoff):
    """Min heavy-atom distance per residue, plus all close atom pairs."""
    xa, ia = arr(resA)
    xb, ib = arr(resB)
    d = np.sqrt(((xa[:, None, :] - xb[None, :, :]) ** 2).sum(-1))
    minA = np.full(len(resA), np.inf)
    minB = np.full(len(resB), np.inf)
    np.minimum.at(minA, ia, d.min(axis=1))
    np.minimum.at(minB, ib, d.min(axis=0))
    close = np.argwhere(d < cutoff)
    return minA, minB, d, close, ia, ib


def numbering_gaps(res):
    gaps = []
    for i in range(len(res) - 1):
        a, b = res[i], res[i + 1]
        if b.num != a.num + 1:
            gaps.append((a.num, b.num, b.num - a.num - 1))
    return gaps


def chain_breaks(res, cut=2.0):
    out = []
    for i in range(len(res) - 1):
        a, b = res[i], res[i + 1]
        c, n = a.coord("C"), b.coord("N")
        if c is None or n is None:
            out.append((a, b, None))
            continue
        d = float(np.linalg.norm(c - n))
        if d > cut:
            out.append((a, b, d))
    return out


def completeness(res):
    nonstd, miss_bb, miss_sc = [], [], []
    for r in res:
        if r.name not in RES_HEAVY:
            nonstd.append(r)
            continue
        have = set(a.name for a in r.heavy())
        exp = set(RES_HEAVY[r.name])
        mb = [a for a in BACKBONE if a not in have]
        ms = sorted(exp - have - set(BACKBONE))
        extra = sorted(have - exp - {"OXT"})
        if mb:
            miss_bb.append((r, mb))
        if ms or extra:
            miss_sc.append((r, ms, extra))
    return nonstd, miss_bb, miss_sc


def disulfides(res_all, cut=2.5):
    sgs = [(r, r.coord("SG")) for r in res_all
           if r.name == "CYS" and r.coord("SG") is not None]
    out = []
    for i in range(len(sgs)):
        for j in range(i + 1, len(sgs)):
            d = float(np.linalg.norm(sgs[i][1] - sgs[j][1]))
            if d < cut:
                out.append((sgs[i][0], sgs[j][0], d))
    return sorted(out, key=lambda t: t[2])


def membrane_axis(res, exclude=()):
    """First PC of CA coords = approximate membrane normal for a 7TM bundle."""
    ex = set(exclude)
    ca = np.array([r.coord("CA") for r in res
                   if r.coord("CA") is not None and r.num not in ex])
    c = ca.mean(0)
    u, s, vt = np.linalg.svd(ca - c, full_matrices=False)
    return c, vt[0], (s ** 2) / (s ** 2).sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--label", default=None)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--clash", type=float, default=2.0)
    ap.add_argument("--fusion", default=None,
                    help="start-end residue range of a fusion domain in B")
    a = ap.parse_args()
    label = a.label or os.path.basename(a.pdb)

    atoms, notes = read_pdb(a.pdb)
    ch = build(atoms)
    A, B = ch["A"], ch["B"]
    out = []
    P = out.append

    P("## %s" % label)
    P("")
    P("`%s`" % a.pdb)
    P("")

    # ---- 1. per chain ------------------------------------------------
    P("### 1. Per-chain composition")
    P("")
    P("| chain | residues | heavy atoms | all atoms | first | last | numbering gaps |")
    P("|---|---|---|---|---|---|---|")
    for cid, res in ch.items():
        g = numbering_gaps(res)
        gs = "none" if not g else "; ".join("%d->%d (%d missing)" % x for x in g)
        P("| %s | %d | %d | %d | %d | %d | %s |" % (
            cid, len(res), sum(len(r.heavy()) for r in res),
            sum(len(r.atoms) for r in res), res[0].num, res[-1].num, gs))
    P("")

    # ---- 2. interface ------------------------------------------------
    minA, minB, d, close, ia, ib = pair_min(A, B, a.cutoff)
    ifA = [(A[i], minA[i]) for i in range(len(A)) if minA[i] < a.cutoff]
    ifB = [(B[i], minB[i]) for i in range(len(B)) if minB[i] < a.cutoff]
    P("### 2. A/B interface, %.1f A heavy-atom cutoff" % a.cutoff)
    P("")
    P("Chain A (nanobody) interface residues: **%d**" % len(ifA))
    P("")
    P("| residue | min dist to B (A) |")
    P("|---|---|")
    for r, m in ifA:
        P("| %s | %.2f |" % (r.tag, m))
    P("")
    P("Chain B (GPCR) interface residues: **%d**" % len(ifB))
    P("")
    P("| residue | min dist to A (A) |")
    P("|---|---|")
    for r, m in ifB:
        P("| %s | %.2f |" % (r.tag, m))
    P("")
    P("Compact lists:")
    P("")
    P("- A: `%s`" % ", ".join(r.tag for r, _ in ifA))
    P("- B: `%s`" % ", ".join(r.tag for r, _ in ifB))
    P("")

    # ---- 3. which face -----------------------------------------------
    fus = None
    if a.fusion:
        lo, hi = [int(x) for x in a.fusion.split("-")]
        fus = list(range(lo, hi + 1))
    cB, axis, var = membrane_axis(B, exclude=fus or ())

    def proj(p):
        return float(np.dot(np.asarray(p) - cB, axis))

    # orient: receptor residue 1 is at the extracellular start of TM1
    nterm = proj(B[0].coord("CA"))
    if nterm < 0:
        axis = -axis
    zall = np.array([proj(r.coord("CA")) for r in B
                     if r.coord("CA") is not None and
                     (fus is None or r.num not in fus)])
    zif = np.array([proj(r.coord("CA")) for r, _ in ifB
                    if r.coord("CA") is not None])
    zA = np.array([proj(r.coord("CA")) for r in A if r.coord("CA") is not None])

    P("### 3. Which receptor face does the nanobody bind?")
    P("")
    P("Membrane normal approximated as PC1 of chain-B CA coordinates%s "
      "(PC1 captures %.1f%% of the variance). Sign fixed so that chain-B "
      "residue %d (N-terminal end of TM1, extracellular by topology) is at "
      "positive projection." % (
          " excluding the fusion domain" if fus else "", 100 * var[0], B[0].num))
    P("")
    P("| quantity | projection on normal (A) |")
    P("|---|---|")
    P("| chain B N-term CA (res %d) | %+.1f |" % (B[0].num, proj(B[0].coord("CA"))))
    P("| chain B C-term CA (res %d) | %+.1f |" % (B[-1].num, proj(B[-1].coord("CA"))))
    P("| all chain B CA, mean | %+.1f |" % zall.mean())
    P("| all chain B CA, range | %+.1f .. %+.1f |" % (zall.min(), zall.max()))
    P("| chain B *interface* CA, mean | %+.1f |" % zif.mean())
    P("| chain B interface CA, range | %+.1f .. %+.1f |" % (zif.min(), zif.max()))
    P("| whole nanobody (chain A) CA, mean | %+.1f |" % zA.mean())
    if fus:
        zf = np.array([proj(r.coord("CA")) for r in B
                       if r.num in fus and r.coord("CA") is not None])
        P("| fusion domain CA, mean | %+.1f |" % zf.mean())
    P("")
    side = "EXTRACELLULAR" if zif.mean() > zall.mean() else "INTRACELLULAR"
    P("Interface centroid sits %+.1f A relative to the receptor centroid "
      "along the normal -> **%s face**." % (zif.mean() - zall.mean(), side))
    P("")
    P("Per-interface-residue projections (positive = extracellular side):")
    P("")
    P("`" + ", ".join("%s:%+.0f" % (r.tag, proj(r.coord("CA")))
                      for r, _ in ifB if r.coord("CA") is not None) + "`")
    P("")

    # ---- 4. fusion domain ----------------------------------------------
    if fus:
        fres = [r for r in B if r.num in fus]
        mA, mF, _, fclose, _, _ = pair_min(A, fres, a.cutoff)
        dmin = min(mF)
        P("### 4. Fusion domain")
        P("")
        P("Residues B %d-%d (%d residues), sequence `%s`."
          % (fus[0], fus[-1], len(fres), "".join(r.one for r in fres)))
        P("")
        P("Minimum heavy-atom distance from the fusion domain to chain A: "
          "**%.2f A** -> %s contact with the nanobody at %.1f A."
          % (dmin, "NO" if dmin >= a.cutoff else "YES, there IS", a.cutoff))
        if dmin < a.cutoff:
            for i, r in enumerate(fres):
                if mF[i] < a.cutoff:
                    P("- %s  %.2f A" % (r.tag, mF[i]))
        P("")

    # ---- 5. clashes ------------------------------------------------------
    xa, iaa = arr(A)
    xb, ibb = arr(B)
    ha = [at for r in A for at in r.heavy()]
    hb = [at for r in B for at in r.heavy()]
    cl = np.argwhere(d < a.clash)
    P("### 5. Inter-chain steric clashes (A-B heavy atoms < %.1f A)" % a.clash)
    P("")
    P("Total clashing heavy-atom pairs: **%d**" % len(cl))
    if len(cl):
        rows = sorted(((float(d[i, j]), ha[i], hb[j]) for i, j in cl))
        P("")
        P("| dist (A) | chain A atom | chain B atom |")
        P("|---|---|---|")
        for dd, u, v in rows[:10]:
            P("| %.2f | A %s%d %s | B %s%d %s |" % (
                dd, u.resname, u.resnum, u.name, v.resname, v.resnum, v.name))
    P("")

    # ---- 6. chain breaks --------------------------------------------------
    P("### 6. Backbone chain breaks (C(i)-N(i+1) > 2.0 A)")
    P("")
    any_break = False
    for cid, res in ch.items():
        br = chain_breaks(res)
        for u, v, dd in br:
            any_break = True
            P("- chain %s %s%d -> %s%d : %s" % (
                cid, u.one, u.num, v.one, v.num,
                "missing C or N atom" if dd is None else "%.2f A" % dd))
    if not any_break:
        P("None. Every consecutive residue pair in both chains is bonded.")
    P("")

    # ---- 7. completeness ---------------------------------------------------
    P("### 7. Residue-level sanity")
    P("")
    P("- HETATM records: %d" % len(notes["hetatm"]))
    P("- altloc characters present in source: %r" % sorted(notes["altlocs_seen"]))
    P("- atoms dropped for altloc != blank/A: %d" % notes["dropped_altloc"])
    P("- insertion codes: %d" % len(notes["icodes"]))
    P("- models in source: %s (only model 1 retained)" % sorted(notes["models_seen"]))
    zocc = [at for at in atoms if at.occ <= 0.0]
    P("- atoms with occupancy 0.00: %d" % len(zocc))
    P("- hydrogens present: %d" % sum(1 for at in atoms if at.is_hydrogen))
    P("")
    ifAset = set(r.tag for r, _ in ifA)
    ifBset = set(r.tag for r, _ in ifB)
    ifset = ifAset | ifBset
    for cid, res in ch.items():
        nonstd, mbb, msc = completeness(res)
        P("Chain %s:" % cid)
        P("")
        P("- nonstandard residue names: %s"
          % ("none" if not nonstd else ", ".join(r.tag for r in nonstd)))
        if mbb:
            P("- residues missing backbone atoms (%d):" % len(mbb))
            for r, m in mbb:
                P("  - %s missing %s%s" % (r.tag, ",".join(m),
                                           "  **AT INTERFACE**"
                                           if r.tag in ifset else ""))
        else:
            P("- residues missing backbone atoms: none")
        if msc:
            P("- residues with incomplete/unexpected side chains (%d):" % len(msc))
            for r, m, extra in msc:
                bits = []
                if m:
                    bits.append("missing " + ",".join(m))
                if extra:
                    bits.append("unexpected " + ",".join(extra))
                P("  - %s %s%s" % (r.tag, "; ".join(bits),
                                   "  **AT INTERFACE**" if r.tag in ifset else ""))
        else:
            P("- residues with incomplete side chains: none")
        P("")

    # ---- 8. disulfides -------------------------------------------------------
    P("### 8. Disulfides (SG-SG < 2.5 A)")
    P("")
    ss = disulfides([r for res in ch.values() for r in res])
    if not ss:
        P("None found.")
    else:
        for u, v, dd in ss:
            kind = "intra-chain %s" % u.chain if u.chain == v.chain \
                else "INTER-CHAIN %s/%s" % (u.chain, v.chain)
            P("- %s - %s : %.2f A  (%s)" % (u.tag, v.tag, dd, kind))
    free = [r.tag for res in ch.values() for r in res
            if r.name == "CYS" and not any(
                r in (u, v) for u, v, _ in ss)]
    P("")
    P("- free cysteines: %s" % (", ".join(free) if free else "none"))
    P("")
    print("\n".join(out))


if __name__ == "__main__":
    main()
