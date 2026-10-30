#!/usr/bin/env python3
"""Clash / interface counts for the CLDN6 placement from each PyMOL method."""
import os, numpy as np
from scipy.spatial import cKDTree
BASE="/scratch/users/bvbell/xlink"; WORK=os.path.join(BASE,"work")

def read(path, chains=None):
    out=[]
    for l in open(path):
        if l.startswith(("ATOM","HETATM")) and (chains is None or l[21] in chains):
            out.append((l[21],int(l[22:26]),l[17:20].strip(),l[12:16].strip(),
                        (float(l[30:38]),float(l[38:46]),float(l[46:54]))))
    return out

fab=read(os.path.join(BASE,"9v32_complex.pdb"),"HL")
mob=read(os.path.join(BASE,"cldn6.pdb"))
fxyz=np.array([a[4] for a in fab]); ftree=cKDTree(fxyz)

print("%-10s %8s %8s %8s   %-12s %-12s" % ("method","<2.0","<2.5","<3.0","minDist","CLDN6 iface res @5A"))
for name,fn in [("super","transform_super.txt"),("super_CA","transform_super_ca.txt"),
                ("align","transform_align.txt"),("cealign","transform_cealign.txt")]:
    M=[float(x) for x in open(os.path.join(WORK,fn)).read().split()]
    R=np.array(M).reshape(4,4)[:3,:3]; T=np.array([M[3],M[7],M[11]])
    mx=np.array([a[4] for a in mob])@R.T+T
    mt=cKDTree(mx)
    counts=[]
    for cut in (2.0,2.5,3.0):
        counts.append(sum(len(j) for j in mt.query_ball_tree(ftree,cut)))
    d,_=mt.query(fxyz,k=1)
    iface=set()
    for i,js in enumerate(mt.query_ball_tree(ftree,5.0)):
        if js: iface.add(mob[i][1])
    print("%-10s %8d %8d %8d   %-12.2f %-12d" % (name,counts[0],counts[1],counts[2],d.min(),len(iface)))
