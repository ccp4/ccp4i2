"""What a link dictionary says changes in each monomer when the link forms.

AceDRG does not rewrite a monomer to describe a link. The link dictionary
keeps the monomer's own description (its ``data_comp_<id>`` block, when it
carries one at all) exactly as it was, and adds a modification block,
``data_mod_<id>m1``, listing what differs for a residue joined by this
link: atoms added or deleted, bond orders and lengths changed, the angles,
torsions and planes that follow. Refinement applies the modification to the
residues a LINK record names; a viewer drawing from the monomer block alone
shows the unlinked chemistry. So a vinyl made single by the link still reads
``DOUBLE`` in the monomer block, and single in the modification -- which is
why the report lays the modification out rather than leaving anyone to find
it in the file.

Read with gemmi; no CCP4 needed.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class AtomChange:
    """An atom added to or deleted from the monomer."""
    function: str            # "add" | "delete"
    name: str
    element: str
    bonded_to: str = ""      # the atom it hangs off, when a bond row says


@dataclass
class BondChange:
    """A bond between heavy atoms whose order or length the link changes."""
    function: str            # "add" | "change" | "delete"
    atom_1: str
    atom_2: str
    type: str = ""           # as written: single, double, ...
    dist: str = ""           # Angstrom, as written
    old_type: str = ""       # from the monomer block, when the file has one
    old_dist: str = ""


@dataclass
class MonomerModification:
    comp_id: str
    mod_id: str
    atoms: List[AtomChange] = field(default_factory=list)
    bonds: List[BondChange] = field(default_factory=list)
    retyped: List[str] = field(default_factory=list)   # atoms whose energy type changes
    revised: Dict[str, int] = field(default_factory=dict)  # "angles" -> count, ...


def _rows(block, prefix: str, tags: List[str]):
    table = block.find(prefix, tags)
    if not table:
        return []
    return [[_value(row[i]) for i in range(len(tags))] for row in table]


def _value(raw: str) -> str:
    text = raw.strip()
    if text in (".", "?"):
        return ""
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1]
    return text


def _is_hydrogen(name: str, elements: Dict[str, str]) -> bool:
    element = elements.get(name)
    if element:
        return element.upper() in ("H", "D")
    return bool(name) and name[0] in "HD"


def link_modifications(path) -> List[MonomerModification]:
    """The modifications in a link dictionary, one per monomer, in file order.

    An empty list for a file with no ``_chem_mod`` blocks (or no file).
    """
    import gemmi

    path = Path(path)
    if not path.is_file():
        return []
    doc = gemmi.cif.read(str(path))

    comp_of: Dict[str, str] = {}
    mod_list = doc.find_block("mod_list")
    if mod_list is not None:
        for mod_id, comp_id in _rows(mod_list, "_chem_mod.", ["id", "comp_id"]):
            comp_of[mod_id] = comp_id

    mods = []
    for block in doc:
        if not block.name.startswith("mod_") or block.name == "mod_list":
            continue
        mod_id = block.name[len("mod_"):]
        comp_id = comp_of.get(mod_id, mod_id)
        mod = MonomerModification(comp_id=comp_id, mod_id=mod_id)

        # The monomer's own bonds, for a before-and-after. Only a monomer
        # given as a file travels in the link dictionary; a library monomer
        # is named, not copied.
        elements: Dict[str, str] = {}
        old_bonds: Dict[tuple, tuple] = {}
        comp = doc.find_block("comp_" + comp_id)
        if comp is not None:
            for name, element in _rows(comp, "_chem_comp_atom.", ["atom_id", "type_symbol"]):
                elements[name] = element
            # The library writes the order as "type"; AceDRG's own monomer
            # block as PDBx "value_order".
            for order_tag in ("type", "value_order"):
                for a1, a2, kind, dist in _rows(
                        comp, "_chem_comp_bond.", ["atom_id_1", "atom_id_2", order_tag, "value_dist"]):
                    old_bonds[frozenset((a1, a2))] = (kind, dist)

        atom_rows = _rows(block, "_chem_mod_atom.",
                          ["function", "atom_id", "new_atom_id", "new_type_symbol"])
        bond_rows = _rows(block, "_chem_mod_bond.",
                          ["function", "atom_id_1", "atom_id_2", "new_type", "new_value_dist"])
        for function, name, new_name, element in atom_rows:
            name = new_name if function == "add" else name
            if element:
                elements[name] = element

        for function, name, new_name, element in atom_rows:
            if function == "change":
                mod.retyped.append(name)
                continue
            name = new_name if function == "add" else name
            bonded_to = ""
            for bond_function, a1, a2, _kind, _dist in bond_rows:
                if bond_function == function and name in (a1, a2):
                    bonded_to = a2 if a1 == name else a1
                    break
            mod.atoms.append(AtomChange(function, name, element, bonded_to))

        for function, a1, a2, kind, dist in bond_rows:
            if _is_hydrogen(a1, elements) or _is_hydrogen(a2, elements):
                continue   # an added or deleted hydrogen is told by its atom row
            old_type, old_dist = old_bonds.get(frozenset((a1, a2)), ("", ""))
            mod.bonds.append(BondChange(function, a1, a2, kind, dist, old_type, old_dist))

        for label, prefix in (("angles", "_chem_mod_angle."), ("torsions", "_chem_mod_tor."),
                              ("chiral centres", "_chem_mod_chir.")):
            count = len(_rows(block, prefix, ["function"]))
            if count:
                mod.revised[label] = count
        planes = {row[0] for row in _rows(block, "_chem_mod_plane_atom.", ["plane_id"])}
        if planes:
            mod.revised["planes"] = len(planes)
        mods.append(mod)
    return mods


def describe_bond(bond: BondChange) -> str:
    """One short phrase for a bond change: 'double 1.306 A -> single 1.520 A'."""
    now = " ".join(part for part in (bond.type, bond.dist and bond.dist + " A") if part)
    if bond.function == "delete":
        return "deleted"
    if bond.function == "add":
        return "added, " + now if now else "added"
    was = " ".join(part for part in (bond.old_type, bond.old_dist and bond.old_dist + " A") if part)
    return f"{was} -> {now}" if was else now
