SRC="/scratch/users/bvbell/xlink/flexddg_ala/9pi9_clean.pdb"
KEEP={"E":range(175,252),"A":range(101,108),"B":range(91,96)}
out=[]; seen={}; alt=set(); ins=set()
for l in open(SRC):
    if not l.startswith("ATOM"): continue
    ch=l[21]; r=int(l[22:26]); icode=l[26]
    if ch not in KEEP or r not in KEEP[ch]: continue
    if icode.strip(): ins.add((ch,r,icode))
    a=l[16]
    if a.strip(): alt.add(a)
    if a.strip() and a!="A": continue          # keep altloc A only
    out.append(l[:16]+" "+l[17:])              # blank the altloc column
    seen.setdefault(ch,set()).add(r)

# write, TER between chains, preserving E,A,B order
with open("inputs/9pi9_motif_H3L3.pdb","w") as f:
    prev=None; n=0
    for l in out:
        if prev and l[21]!=prev: f.write("TER\n")
        f.write(l); prev=l[21]; n+=1
    f.write("TER\nEND\n")
print("altlocs seen:", sorted(alt) or "none")
print("insertion codes:", sorted(ins) or "none")
for ch in ("E","A","B"):
    s=sorted(seen.get(ch,[]))
    missing=[r for r in KEEP[ch] if r not in seen.get(ch,())]
    print(f"  chain {ch}: {len(s)} res  {min(s)}-{max(s)}  missing={missing or 'none'}")
print("atoms written:", n)
