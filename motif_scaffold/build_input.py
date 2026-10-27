SRC="/scratch/users/bvbell/xlink/flexddg_ala/9pi9_clean.pdb"
OUT="/scratch/users/bvbell/xlink/motif_scaffold/inputs/9pi9_H3L3_E.pdb"
# margin around motifs so the indexing sweep can extend boundaries
KEEP={"A":(95,112), "B":(86,100), "E":(175,251)}
n=0; kept={}
with open(OUT,"w") as o:
    for l in open(SRC):
        if not l.startswith("ATOM"): continue
        ch=l[21]; r=int(l[22:26]); alt=l[16]
        if ch not in KEEP: continue
        lo,hi=KEEP[ch]
        if not (lo<=r<=hi): continue
        if alt not in (" ","A"): continue
        o.write(l[:16]+" "+l[17:]); n+=1
        kept.setdefault(ch,set()).add(r)
    o.write("TER\nEND\n")
print(f"wrote {OUT}  ({n} atoms)")
for ch in "ABE":
    rs=sorted(kept.get(ch,[]))
    print(f"  chain {ch}: {len(rs)} res, {rs[0]}-{rs[-1]}, contiguous={rs==list(range(rs[0],rs[-1]+1))}")
