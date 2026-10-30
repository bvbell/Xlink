#!/usr/bin/env python3
"""Convert the two nanobody/GPCR complexes into flex-ddG-ready PDBs.

  INPUT 1  AF3 mmCIF  -> prep/d006_complex.pdb
  INPUT 2  Coot PDB   -> prep/rosetta9_complex.pdb

Faithful conversion only. Auth chain IDs and auth residue numbers preserved.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pdbio import read_cif_atom_site, read_pdb, write_pdb  # noqa: E402

ROOT = "/scratch/users/bvbell/xlink"
PREP = os.path.join(ROOT, "prep")

JOBS = [
    ("d006 (AF3 mmCIF)",
     os.path.join(ROOT, "d006_lmpnn_fkd2_run7_af3 (1).cif"),
     os.path.join(PREP, "d006_complex.pdb"), "cif"),
    ("rosetta9 (Coot PDB)",
     os.path.join(ROOT,
                  "Rosetta_9_real_trunc_model_1_flip_fit_J126map-coot-0_flip.pdb"),
     os.path.join(PREP, "rosetta9_complex.pdb"), "pdb"),
]


def main():
    os.makedirs(PREP, exist_ok=True)
    for label, src, dst, kind in JOBS:
        reader = read_cif_atom_site if kind == "cif" else read_pdb
        atoms, notes = reader(src)
        fixed = write_pdb(atoms, dst, chain_order=["A", "B"])
        print("=== %s ===" % label)
        print("  src          : %s" % src)
        print("  dst          : %s" % dst)
        print("  atoms kept   : %d" % len(atoms))
        print("  models seen  : %s" % sorted(notes["models_seen"]))
        print("  altlocs seen : %r" % sorted(notes["altlocs_seen"]))
        print("  dropped (altloc/model): %d / %d"
              % (notes["dropped_altloc"], notes["dropped_model"]))
        print("  HETATM       : %d" % len(notes["hetatm"]))
        print("  insertion codes: %d" % len(notes["icodes"]))
        if fixed:
            print("  !! ZERO OCCUPANCY REPAIRED on %d atoms:" % len(fixed))
            for a, old in fixed:
                print("     %s %s%d %s (occ %.2f -> 1.00)"
                      % (a.chain, a.resname, a.resnum, a.name, old))
        else:
            print("  zero-occupancy atoms: none")
        print()


if __name__ == "__main__":
    main()
