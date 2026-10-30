from pymol import cmd
cmd.load("/scratch/users/bvbell/xlink/9v32_complex.pdb","cplx")
cmd.create("tgtA","cplx and chain A and polymer")
cmd.load("/scratch/users/bvbell/xlink/cldn6.pdb","mob")
cmd.dss()
aa={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
for obj in ("tgtA","mob"):
    rows=[]
    cmd.iterate(obj+" and name CA","rows.append((int(resi),resn,ss))",space={"rows":rows,"int":int})
    print("=== %s n=%d range %d-%d"%(obj,len(rows),rows[0][0],rows[-1][0]))
    print("SEQ "+"".join(aa.get(r[1],'X') for r in rows))
    print("RESI "+",".join(str(r[0]) for r in rows))
    segs=[]
    for ri,rn,ss in rows:
        s=ss if ss else 'L'
        if segs and segs[-1][0]==s and ri==segs[-1][2]+1: segs[-1][2]=ri
        else: segs.append([s,ri,ri])
    print("SS  "+" ".join("%s%d-%d"%(s,a,b) for s,a,b in segs if s!='L'))
    nums=[r[0] for r in rows]
    print("GAPS "+str([(nums[i],nums[i+1]) for i in range(len(nums)-1) if nums[i+1]!=nums[i]+1]))
