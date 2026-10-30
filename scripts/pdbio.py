"""Minimal, dependency-light PDB / mmCIF atom_site reader + strict PDB writer.

Used to prepare nanobody/GPCR complexes for Rosetta flex ddG.
Faithful format conversion only: no renumbering, no atom addition/removal
(other than the documented altloc / model filtering), no minimisation.
"""

from collections import OrderedDict

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Expected heavy (non-H) side-chain + backbone atom sets, excluding OXT.
RES_HEAVY = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2",
            "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2",
            "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}

BACKBONE = ("N", "CA", "C", "O")


class Atom(object):
    __slots__ = ("record", "name", "altloc", "resname", "chain", "resnum",
                 "icode", "x", "y", "z", "occ", "bfac", "element")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    @property
    def is_hydrogen(self):
        return self.element.strip().upper() in ("H", "D")

    @property
    def reskey(self):
        return (self.chain, self.resnum, self.icode)


def _guess_element(name, resname):
    n = name.strip()
    if len(n) == 4 and n[0].isdigit():
        n = n[1:]
    if n[:2].upper() in ("CL", "BR", "ZN", "FE", "MG", "MN", "NA", "SE"):
        return n[:2].upper()
    return n[0].upper()


def read_pdb(path, keep_hetatm=True):
    """Return (atoms, notes). Only MODEL 1 is kept; altloc must be '' or 'A'."""
    atoms = []
    notes = {"models_seen": set(), "altlocs_seen": set(), "dropped_altloc": 0,
             "dropped_model": 0, "hetatm": [], "icodes": []}
    model = 1
    in_model = True
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec == "MODEL ":
                model = int(line[10:14])
                notes["models_seen"].add(model)
                in_model = (model == 1)
                continue
            if rec == "ENDMDL":
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue
            if not in_model:
                notes["dropped_model"] += 1
                continue
            notes["models_seen"].add(model)
            alt = line[16]
            notes["altlocs_seen"].add(alt)
            if alt not in (" ", "A"):
                notes["dropped_altloc"] += 1
                continue
            name = line[12:16]
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])
            icode = line[26]
            element = line[76:78].strip() if len(line) >= 78 else ""
            if not element:
                element = _guess_element(name, resname)
            occ = line[54:60].strip()
            bfac = line[60:66].strip()
            a = Atom(record=rec.strip(), name=name.strip(), altloc=alt,
                     resname=resname, chain=chain, resnum=resnum, icode=icode,
                     x=float(line[30:38]), y=float(line[38:46]),
                     z=float(line[46:54]),
                     occ=float(occ) if occ else 1.0,
                     bfac=float(bfac) if bfac else 0.0,
                     element=element)
            if rec == "HETATM":
                notes["hetatm"].append(a)
                if not keep_hetatm:
                    continue
            if icode != " ":
                notes["icodes"].append(a)
            atoms.append(a)
    return atoms, notes


def read_cif_atom_site(path):
    """Parse the _atom_site loop of a simple (AF3-style) mmCIF."""
    tags = []
    rows = []
    state = 0  # 0 outside, 1 reading tags, 2 reading rows
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if state == 0:
                if s.startswith("_atom_site."):
                    tags = [s]
                    state = 1
                continue
            if state == 1:
                if s.startswith("_atom_site."):
                    tags.append(s)
                    continue
                state = 2
                # fall through to row handling
            if state == 2:
                if s == "#" or s.startswith("loop_") or s.startswith("data_"):
                    break
                if not s:
                    continue
                rows.append(s.split())
    tags = [t.split(".", 1)[1] for t in tags]
    idx = {t: i for i, t in enumerate(tags)}

    def g(row, key, default=None):
        if key not in idx:
            return default
        v = row[idx[key]]
        return default if v in (".", "?") else v

    atoms = []
    notes = {"models_seen": set(), "altlocs_seen": set(), "dropped_altloc": 0,
             "dropped_model": 0, "hetatm": [], "icodes": []}
    for row in rows:
        if len(row) != len(tags):
            raise ValueError("atom_site row has %d fields, expected %d: %r"
                             % (len(row), len(tags), row))
        model = int(g(row, "pdbx_PDB_model_num", "1"))
        notes["models_seen"].add(model)
        if model != 1:
            notes["dropped_model"] += 1
            continue
        alt = g(row, "label_alt_id", None)
        notes["altlocs_seen"].add(" " if alt is None else alt)
        if alt is not None and alt != "A":
            notes["dropped_altloc"] += 1
            continue
        icode = g(row, "pdbx_PDB_ins_code", None) or " "
        chain = g(row, "auth_asym_id") or g(row, "label_asym_id")
        resnum = int(g(row, "auth_seq_id") or g(row, "label_seq_id"))
        resname = g(row, "label_comp_id")
        occ_s = g(row, "occupancy", "1.00")
        b_s = g(row, "B_iso_or_equiv", "0.00")
        a = Atom(record=row[idx["group_PDB"]] if "group_PDB" in idx else "ATOM",
                 name=g(row, "label_atom_id"),
                 altloc=" " if alt is None else alt,
                 resname=resname, chain=chain, resnum=resnum, icode=icode,
                 x=float(g(row, "Cartn_x")), y=float(g(row, "Cartn_y")),
                 z=float(g(row, "Cartn_z")),
                 occ=float(occ_s), bfac=float(b_s),
                 element=(g(row, "type_symbol") or
                          _guess_element(g(row, "label_atom_id"), resname)))
        if a.record == "HETATM":
            notes["hetatm"].append(a)
        if icode != " ":
            notes["icodes"].append(a)
        atoms.append(a)
    return atoms, notes


def _fmt_atom_name(name, element):
    """PDB column 13-16 convention."""
    n = name.strip()
    if len(n) >= 4:
        return n[:4]
    el = element.strip().upper()
    if len(el) == 2 and n.upper().startswith(el):
        return "%-4s" % n
    return " %-3s" % n


def write_pdb(atoms, path, chain_order=None, fix_zero_occ=True):
    """Write strict column-formatted ATOM records, TER between chains, END.

    Returns a list of (atom, old_occ) for atoms whose occupancy was repaired.
    """
    fixed = []
    bychain = OrderedDict()
    for a in atoms:
        bychain.setdefault(a.chain, []).append(a)
    if chain_order:
        order = [c for c in chain_order if c in bychain]
        order += [c for c in bychain if c not in order]
    else:
        order = list(bychain)

    serial = 0
    lines = []
    for ch in order:
        last = None
        for a in bychain[ch]:
            serial += 1
            occ = a.occ
            if fix_zero_occ and occ <= 0.0:
                fixed.append((a, occ))
                occ = 1.00
            lines.append(
                "ATOM  %5d %s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
                "          %2s  " % (
                    serial % 100000, _fmt_atom_name(a.name, a.element),
                    " ", a.resname[:3].rjust(3), a.chain, a.resnum, a.icode,
                    a.x, a.y, a.z, occ, a.bfac, a.element.strip().upper()[:2]
                ))
            last = a
        serial += 1
        lines.append("TER   %5d      %3s %1s%4d%1s" % (
            serial % 100000, last.resname[:3].rjust(3), last.chain,
            last.resnum, last.icode))
    lines.append("END")
    with open(path, "w") as fh:
        # pad every record to the 80-column PDB width
        fh.write("\n".join("%-80s" % ln for ln in lines) + "\n")
    return fixed
