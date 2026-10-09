"""Hydrogens a search model should not carry into refinement.

A search model's hydrogens are its template's. After Chainsaw mutates side
chains to the target's residues their names no longer belong to the residue
they sit on (glycine's HA2/HA3 on what is now a serine), and REFMAC refuses
the model (error 350, an atom its dictionary does not know for that
residue). An agent run (Opus, 2026-10-09) hit this on a Chainsaw model taken
from MrBUMP's working directory and stripped the hydrogens by hand.

A model that has just been placed by molecular replacement never carries
hydrogens worth keeping, so the MR pipelines strip them before Sheetbend and
REFMAC, and the chainsaw task strips them from what it writes. Nothing here
touches a refined model: riding hydrogens from a previous refinement are
legitimate input, and only the dictionary can tell a right name from a
wrong one.

gemmi only; the PDB-record test is stdlib.
"""
import gemmi


def is_hydrogen_record(line: str) -> bool:
    """Whether a PDB ATOM/HETATM record is a hydrogen (or deuterium): by the
    element columns when they are filled, else by the atom name, where a
    hydrogen's name starts in column 14 (" HA ", "1HB ") or fills four
    columns from 13 ("HD21"), and a two-letter element at column 13 ("HG  ")
    is mercury, not a hydrogen."""
    if not line.startswith(("ATOM", "HETATM")):
        return False
    element = line[76:78].strip().upper()
    if element:
        return element in ("H", "D")
    name = line[12:16]
    if len(name) < 4:
        return False
    if name[0] in " 0123456789":
        return name[1].upper() in ("H", "D")
    return name[0].upper() in ("H", "D") and len(name.strip()) == 4


def strip_hydrogens(source, destination) -> int:
    """Write ``source`` (PDB or mmCIF) to ``destination`` without hydrogens;
    returns how many were removed. ``destination`` is written as PDB."""
    structure = gemmi.read_structure(str(source), format=gemmi.CoorFormat.Detect)
    before = structure[0].count_atom_sites() if len(structure) else 0
    structure.remove_hydrogens()
    after = structure[0].count_atom_sites() if len(structure) else 0
    structure.setup_entities()
    structure.write_pdb(str(destination))
    return before - after
