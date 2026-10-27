"""Score runs_prod designs: motif RMSD, backbone clashes, broken bonds.

Usage: python score_prod.py [run_dir] [out_csv]
Outputs are backbone-only (N CA C O); chain A = designed peptide, chain B = target.
"""
import numpy as np, glob, os, sys, csv
BASE = "/scratch/users/bvbell/xlink/motif_scaffold"
REF = f"{BASE}/inputs/9pi9_motif_H3L3_pad.pdb"
RUN = sys.argv[1] if len(sys.argv) > 1 else f"{BASE}/runs_prod"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{BASE}/scores_prod.csv"
BB = (" N  ", " CA ", " C  ")
CLASH = 3.0      # heavy-atom backbone pair closer than this = clash ...
CLASH_NO = 2.5   # ... except N-O pairs, which sit at ~2.9 A in backbone H-bonds
CN_BREAK = 2.0   # peptide C(i)-N(i+1) bond; ideal 1.33
CN_STRETCH = 1.5
CACA_MAX = 4.2   # consecutive CA-CA; ideal 3.8

def parse(p):
    """{(chain,resnum): {atom: xyz}} for ATOM records, in file order."""
    res = {}
    for l in open(p):
        if l.startswith("ATOM"):
            k = (l[21], int(l[22:26]))
            res.setdefault(k, {})[l[12:16]] = np.array([float(l[30+8*i:38+8*i]) for i in range(3)])
    return res

def rmsd(A, B): return float(np.sqrt(((A-B)**2).sum(1).mean()))
def superpose(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    V, S, W = np.linalg.svd((mob-mc).T @ (ref-rc))
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ W))]) @ W
    return R, rc - mc @ R, rmsd((mob-mc) @ R + rc, ref)

def breaks(res, chain):
    """Consecutive residues in chain with a stretched C-N or CA-CA distance."""
    ks = sorted(k for k in res if k[0] == chain)
    out, n_str, cn_max, ca_max = [], 0, 0.0, 0.0
    for a, b in zip(ks, ks[1:]):
        cn = np.linalg.norm(res[a][" C  "] - res[b][" N  "])
        ca = np.linalg.norm(res[a][" CA "] - res[b][" CA "])
        cn_max, ca_max = max(cn_max, cn), max(ca_max, ca)
        if cn > CN_BREAK or ca > CACA_MAX:
            out.append(f"{chain}{a[1]}-{b[1]}:{cn:.2f}")
        elif cn > CN_STRETCH:
            n_str += 1
    return out, n_str, cn_max, ca_max

def coords(res, keys):
    return np.array([xyz for k in keys for xyz in res[k].values()])

def cutoffs(ra, rb, keys_a, keys_b):
    """Per-pair clash cutoff: CLASH_NO for N-O pairs, CLASH otherwise."""
    ea = np.array([a.strip()[0] for k in keys_a for a in ra[k]])
    eb = np.array([a.strip()[0] for k in keys_b for a in rb[k]])
    no = ((ea[:, None] == "N") & (eb[None] == "O")) | ((ea[:, None] == "O") & (eb[None] == "N"))
    return np.where(no, CLASH_NO, CLASH)

ref = parse(REF)
rows = []
for trb in sorted(glob.glob(f"{RUN}/*.trb"), key=lambda p: int(p.split("_")[-1][:-4])):
    pdb = trb[:-4] + ".pdb"
    if not os.path.exists(pdb): continue
    d = np.load(trb, allow_pickle=True)
    des = parse(pdb)
    mot = list(zip(map(tuple, d["con_ref_pdb_idx"]), map(tuple, d["con_hal_pdb_idx"])))
    tgt = list(zip(map(tuple, d["receptor_con_ref_pdb_idx"]), map(tuple, d["receptor_con_hal_pdb_idx"])))

    # superpose on the fixed target, then measure the motif in that frame
    R, t, tfit = superpose(np.array([des[h][a] for r, h in tgt for a in BB]),
                           np.array([ref[r][a] for r, h in tgt for a in BB]))
    Mm = np.array([des[h][a] for r, h in mot for a in BB]) @ R + t
    Mr = np.array([ref[r][a] for r, h in mot for a in BB])
    _, _, mot_int = superpose(Mm, Mr)

    # clashes: peptide-target backbone pairs, and peptide-peptide pairs |i-j|>=3
    pch = {h[0] for _, h in mot}
    pep = sorted(k for k in des if k[0] in pch)
    trg = sorted(k for k in des if k[0] not in pch)
    P, T = coords(des, pep), coords(des, trg)
    D = np.linalg.norm(P[:, None] - T[None], axis=2)
    Cpt, Cpp = cutoffs(des, des, pep, trg), cutoffs(des, des, pep, pep)
    pidx = np.repeat(np.arange(len(pep)), [len(des[k]) for k in pep])
    Dpp = np.linalg.norm(P[:, None] - P[None], axis=2)
    far = np.abs(pidx[:, None] - pidx[None]) >= 3
    n_intra = int(((Dpp < Cpp) & far).sum() // 2)
    clash_res = sorted({pep[i][1] for i in pidx[(D < Cpt).any(1)]})

    pbr, n_str, cn_max, ca_max = breaks(des, "A")
    tbr, _, _, _ = breaks(des, "B")
    linker = [x for x in d["sampled_mask"][1].split("/") if x[0].isdigit()]
    rows.append(dict(
        design=os.path.basename(pdb)[:-4], pep_len=len(pep),
        linker=linker[0].split("-")[0] if linker else "",
        plddt_mean=round(float(np.mean(d["plddt"][-1])) if "plddt" in d else float("nan"), 3),
        tgt_fit=round(tfit, 3),
        motif_rmsd_docked=round(rmsd(Mm, Mr), 3), motif_rmsd_internal=round(mot_int, 3),
        min_dist=round(float(D.min()), 2),
        n_clash_inter=int((D < Cpt).sum()), n_clash_intra=n_intra,
        clash_res=";".join(map(str, clash_res)),
        max_CN=round(cn_max, 2), max_CACA=round(ca_max, 2),
        n_stretched=n_str, n_break=len(pbr) + len(tbr), breaks=";".join(pbr + tbr)))

if not rows: sys.exit("no designs found in " + RUN)
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

col = lambda k: np.array([r[k] for r in rows], float)
ok = (col("n_break") == 0) & (col("n_clash_inter") == 0) & (col("n_clash_intra") == 0)
print(f"{len(rows)} designs -> {OUT}")
for k in ("motif_rmsd_docked", "motif_rmsd_internal", "tgt_fit", "min_dist", "max_CN", "max_CACA"):
    v = col(k); print(f"  {k:20s} min {v.min():.2f}  median {np.median(v):.2f}  max {v.max():.2f}")
print(f"  broken bonds:       {int((col('n_break') > 0).sum())} designs (C-N > {CN_BREAK} or CA-CA > {CACA_MAX})")
print(f"  stretched C-N:      {int((col('n_stretched') > 0).sum())} designs ({CN_STRETCH}-{CN_BREAK} A)")
print(f"  inter-chain clash:  {int((col('n_clash_inter') > 0).sum())} designs")
print(f"  intra-pep clash:    {int((col('n_clash_intra') > 0).sum())} designs")
print(f"  clean (no break/clash): {int(ok.sum())}")
for r in [r for r in rows if r["n_break"]][:10]: print("  break:", r["design"], r["breaks"])
