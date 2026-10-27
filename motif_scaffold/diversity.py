"""Structural diversity of the RFdiffusion binder designs.

The target (chain B) is held fixed in every design, so all binders already
share a frame.  Two complementary measures over the binder (chain A):

  1. TM-align all-vs-all (superposition-based, handles the 2-6 linker length
     differences).  TM-score normalized by average length (-a T).  Note the
     binders are only ~14-18 residues, where TM-score's d0 bottoms out at
     0.5 A, so it is strict: ~1.5 A RMSD already gives TM ~0.7.
  2. In-frame CA RMSD (no superposition) between designs with the same linker
     length -- this captures placement on the target, which TM-align ignores.

Designs are clustered (average linkage on 1-TM) at --tm_cut.

Usage:  python diversity.py [--dir runs_prod] [--tm_cut 0.8] [--limit N]
Outputs (in --out): diversity_pairs.npz, diversity_clusters.csv, printed summary.
"""
import argparse, glob, os, re, subprocess, tempfile, csv
from multiprocessing import Pool
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

BASE = "/scratch/users/bvbell/xlink/motif_scaffold"
TMALIGN = "/home/groups/possu/miniconda/envs/allatom_design/bin/TMalign"
TM_RE = re.compile(r"TM-score= ([\d.]+) \(if normalized by average")
RMSD_RE = re.compile(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+)")

def binder_ca(pdb, chain="A"):
    lines, ca = [], []
    for l in open(pdb):
        if l.startswith("ATOM") and l[21] == chain:
            lines.append(l)
            if l[12:16] == " CA ":
                ca.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return lines, np.array(ca)

def tmalign(pair):
    i, j, a, b = pair
    out = subprocess.run([TMALIGN, a, b, "-a", "T"], capture_output=True, text=True).stdout
    tm, r = TM_RE.search(out), RMSD_RE.search(out)
    return i, j, float(tm.group(1)), float(r.group(2))

def ncpu():
    try: return len(os.sched_getaffinity(0))
    except AttributeError: return os.cpu_count()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=f"{BASE}/runs_prod")
    ap.add_argument("--out", default=None, help="default: --dir")
    ap.add_argument("--tm_cut", type=float, default=0.8)
    ap.add_argument("--limit", type=int, default=None, help="use first N designs (testing)")
    args = ap.parse_args()
    out = args.out or args.dir

    pdbs = sorted(glob.glob(f"{args.dir}/*.pdb"),
                  key=lambda p: int(re.search(r"_(\d+)\.pdb$", p).group(1)))[:args.limit]
    names = [os.path.basename(p)[:-4] for p in pdbs]
    n = len(pdbs)

    tmp = tempfile.mkdtemp(prefix="div_", dir=os.environ.get("L_SCRATCH") or None)
    files, cas, linker = [], [], []
    for p, nm in zip(pdbs, names):
        lines, ca = binder_ca(p)
        f = f"{tmp}/{nm}.pdb"; open(f, "w").writelines(lines + ["END\n"])
        files.append(f); cas.append(ca)
        trb = np.load(p[:-4] + ".trb", allow_pickle=True)
        linker.append(int(trb["sampled_mask"][-1].split("/")[1].split("-")[0]))
    linker = np.array(linker)

    # 1. TM-align all-vs-all
    tm = np.eye(n); rms = np.zeros((n, n))
    jobs = [(i, j, files[i], files[j]) for i in range(n) for j in range(i + 1, n)]
    print(f"{n} designs, {len(jobs)} TM-align pairs on {ncpu()} CPU(s)")
    with Pool(ncpu()) as pool:
        for k, (i, j, t, r) in enumerate(pool.imap_unordered(tmalign, jobs, chunksize=256)):
            tm[i, j] = tm[j, i] = t; rms[i, j] = rms[j, i] = r
            if (k + 1) % 20000 == 0: print(f"  {k + 1}/{len(jobs)}", flush=True)

    # 2. in-frame CA RMSD within each linker length
    inframe = np.full((n, n), np.nan)
    for L in np.unique(linker):
        idx = np.where(linker == L)[0]
        X = np.stack([cas[i] for i in idx])
        d = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1).mean(-1))
        inframe[np.ix_(idx, idx)] = d

    # clustering on 1 - TM
    D = 1 - tm; np.fill_diagonal(D, 0)
    clus = fcluster(linkage(squareform(D, checks=False), "average"),
                    t=1 - args.tm_cut, criterion="distance") if n > 1 else np.ones(1, int)

    np.savez(f"{out}/diversity_pairs.npz", names=names, tm=tm, rmsd_aligned=rms,
             rmsd_inframe=inframe, linker=linker, cluster=clus)
    with open(f"{out}/diversity_clusters.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "linker", "cluster", "mean_tm_to_others"])
        for i in range(n):
            w.writerow([names[i], linker[i], clus[i], f"{(tm[i].sum() - 1) / max(n - 1, 1):.3f}"])

    iu = np.triu_indices(n, 1)
    sizes = np.bincount(clus)[1:]; sizes = np.sort(sizes[sizes > 0])[::-1]
    print(f"\nTM-score (avg-len norm): mean {tm[iu].mean():.3f}  median {np.median(tm[iu]):.3f}"
          f"  min {tm[iu].min():.3f}")
    print(f"aligned RMSD (A):        mean {rms[iu].mean():.2f}  median {np.median(rms[iu]):.2f}")
    print("linker length: " + "  ".join(f"{L}:{(linker == L).sum()}" for L in np.unique(linker)))
    for L in np.unique(linker):
        m = linker == L
        v = inframe[np.ix_(m, m)][np.triu_indices(m.sum(), 1)]
        if len(v): print(f"  in-frame CA RMSD, linker {L}: mean {v.mean():.2f}  median {np.median(v):.2f} A")
    print(f"clusters at TM>={args.tm_cut}: {len(sizes)}  (largest: {list(sizes[:10])})")
    print(f"wrote {out}/diversity_pairs.npz, {out}/diversity_clusters.csv")

if __name__ == "__main__":
    main()
