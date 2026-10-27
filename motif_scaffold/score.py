import numpy as np, glob, os, sys, csv, re
REF="/scratch/users/bvbell/xlink/motif_scaffold/inputs/9pi9_motif_H3L3_pad.pdb"
BB=(" N  "," CA "," C  ")

def parse(p):
    bb={}; allat={}
    for l in open(p):
        if not l.startswith("ATOM"): continue
        k=(l[21],int(l[22:26])); a=l[12:16]
        xyz=np.array([float(l[30+8*i:38+8*i]) for i in range(3)])
        if a in BB: bb.setdefault(k,{})[a]=xyz
        if l[76:78].strip()!="H": allat.setdefault(k,[]).append(xyz)
    return bb,allat

def rmsd(A,B): return float(np.sqrt(((A-B)**2).sum(1).mean()))
def superpose(mob,ref):
    """Return (R,t) mapping mob->ref, plus fit rmsd."""
    mc,rc=mob.mean(0),ref.mean(0)
    V,S,W=np.linalg.svd((mob-mc).T@(ref-rc))
    d=np.sign(np.linalg.det(V@W))
    R=V@np.diag([1,1,d])@W
    return R, rc-mc@R, rmsd((mob-mc)@R+rc, ref)

refbb,_=parse(REF); rows=[]
for trb in sorted(glob.glob("runs/*/*.trb")):
    pdb=trb[:-4]+".pdb"
    if not os.path.exists(pdb): continue
    cfg=trb.split("/")[1]
    d=np.load(trb,allow_pickle=True)
    ref_idx=[tuple(x) for x in d["con_ref_pdb_idx"]]
    hal_idx=[tuple(x) for x in d["con_hal_pdb_idx"]]
    bb,allat=parse(pdb)
    mot_p=list(zip(ref_idx,hal_idx))
    tgt_p=list(zip([tuple(x) for x in d["receptor_con_ref_pdb_idx"]],
                   [tuple(x) for x in d["receptor_con_hal_pdb_idx"]]))
    mask="/".join(d["sampled_mask"]) if "sampled_mask" in d else ""
    segs=[x for x in mask.split("/") if x and x!="0"]
    realized="+".join(x.split("-")[0] if re.match(r"^\d+-\d+$",x) and x.split("-")[0]==x.split("-")[1]
                      else ("M" if re.match(r"^[A-Z]",x) else x) for x in segs[1:])
    if not tgt_p or not mot_p: continue
    # superpose on the fixed target, then measure the motif in that frame
    Tm=np.array([bb[h][a]    for r,h in tgt_p for a in BB])
    Tr=np.array([refbb[r][a] for r,h in tgt_p for a in BB])
    R,t,tfit=superpose(Tm,Tr)
    Mm=np.array([bb[h][a]    for r,h in mot_p for a in BB])@R+t
    Mr=np.array([refbb[r][a] for r,h in mot_p for a in BB])
    _,_,mot_aln=superpose(Mm,Mr)       # internal geometry only
    pch={h[0] for _,h in mot_p}
    pep=[k for k in bb if k[0] in pch]; tgt=[k for k in bb if k[0] not in pch]
    P=np.vstack([np.array(allat[k]) for k in pep]); T=np.vstack([np.array(allat[k]) for k in tgt])
    D=np.linalg.norm(P[:,None,:]-T[None,:,:],axis=2)
    nif=sum(1 for k in pep if np.linalg.norm(np.array(allat[k])[:,None,:]-T[None,:,:],axis=2).min()<5.0)
    ca=np.array([bb[k][" CA "] for k in pep])
    rows.append(dict(config=cfg, design=os.path.basename(pdb)[:-4], pep_len=len(pep),
        n_motif=len(mot_p), tgt_fit=round(tfit,3), realized=realized,
        motif_rmsd_docked=round(rmsd(Mm,Mr),3), motif_rmsd_internal=round(mot_aln,3),
        min_dist=round(float(D.min()),2), n_contact=nif,
        rg=round(float(np.sqrt(((ca-ca.mean(0))**2).sum(1).mean())),2)))
if not rows: print("no designs yet"); sys.exit()
with open("scores.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); [w.writerow(r) for r in rows]
print(f"{len(rows)} designs -> scores.csv")
for r in rows[:4]: print("  ",r)
