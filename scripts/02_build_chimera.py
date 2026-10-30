#!/usr/bin/env python3
"""
Build the CLDN6 / 9v32-Fab chimera.

Applies the PyMOL `super` transform (written by 01_superpose.py) to every atom
of cldn6.pdb and concatenates it with 9v32 chains H and L, whose coordinates
are copied verbatim from the input file.

    /home/groups/possu/miniconda/bin/python3 scripts/02_build_chimera.py
"""
import os

BASE = "/scratch/users/bvbell/xlink"
WORK = os.path.join(BASE, "work")
OUT = os.path.join(BASE, "cldn6_9v32fab_complex.pdb")

with open(os.path.join(WORK, "transform_super.txt")) as fh:
    M = [float(x) for x in fh.read().split()]
assert len(M) == 16, M

R = [M[0:3], M[4:7], M[8:11]]
T = [M[3], M[7], M[11]]
print("rotation:")
for row in R:
    print("   %10.6f %10.6f %10.6f" % tuple(row))
print("translation: %10.4f %10.4f %10.4f" % tuple(T))
print("last row (should be 0 0 0 1): %s" % (M[12:16],))


def xform(x, y, z):
    return (R[0][0] * x + R[0][1] * y + R[0][2] * z + T[0],
            R[1][0] * x + R[1][1] * y + R[1][2] * z + T[1],
            R[2][0] * x + R[2][1] * y + R[2][2] * z + T[2])


def coord_fmt(v):
    s = "%8.3f" % v
    assert len(s) == 8, (v, s)
    return s


out_lines = []
serial = 0


def emit(line):
    """Re-serialize an ATOM line with the next atom serial number."""
    global serial
    serial += 1
    s = serial if serial <= 99999 else 99999
    out_lines.append(line[:6] + "%5d" % s + line[11:].rstrip("\n").rstrip())


# ---- chain A : transformed CLDN6 (native numbering 1-187) ----
nA = 0
with open(os.path.join(BASE, "cldn6.pdb")) as fh:
    for line in fh:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        nx, ny, nz = xform(x, y, z)
        newline = (line[:21] + "A" + line[22:30]
                   + coord_fmt(nx) + coord_fmt(ny) + coord_fmt(nz)
                   + line[54:].rstrip("\n").rstrip())
        emit(newline)
        nA += 1
last = out_lines[-1]
serial += 1
out_lines.append("TER   %5d      %s %s%4s" % (serial, last[17:20], "A", last[22:26]))

# ---- chains H and L : verbatim 9v32 coordinates ----
counts = {"H": 0, "L": 0}
for ch in ("H", "L"):
    with open(os.path.join(BASE, "9v32_complex.pdb")) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != ch:
                continue
            emit(line)
            counts[ch] += 1
    last = out_lines[-1]
    serial += 1
    out_lines.append("TER   %5d      %s %s%4s" % (serial, last[17:20], ch, last[22:26]))

out_lines.append("END")

with open(OUT, "w") as fh:
    fh.write("\n".join(out_lines) + "\n")

print("wrote %s" % OUT)
print("  chain A (CLDN6) atoms: %d" % nA)
print("  chain H atoms        : %d" % counts["H"])
print("  chain L atoms        : %d" % counts["L"])
