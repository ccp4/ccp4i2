"""The monomers in a model that hold a metal, for MetalCoord.

MetalCoord makes restraints for one monomer code at a time, so servalcat_pipe
needs the codes of every residue holding a metal atom: single ions (ZN, CA,
MG) and metal-containing monomers (HEM, SF4) alike. The Qt interface found
them by walking the model for atoms whose element is a metal; this is that
walk, shared by the file digest the task interface reads and the pipeline's
own fallback when no codes were chosen. The composition digest's ligand list
cannot stand in for it: it leaves out single ions and anything of fewer than
five atoms, which is most metal sites.

gemmi only.
"""
from typing import List

import gemmi


def metal_monomer_codes(structure: gemmi.Structure) -> List[str]:
    """Sorted residue names, in the first model, of residues with at least
    one metal atom."""
    if len(structure) == 0:
        return []
    return sorted({
        residue.name
        for chain in structure[0]
        for residue in chain
        if any(atom.element.is_metal for atom in residue)
    })


def metal_monomer_codes_in_file(path: str) -> List[str]:
    """As metal_monomer_codes, from a coordinate file (PDB or mmCIF). An
    unreadable file has none; gemmi returns an empty structure for some
    files rather than raising, which comes to the same thing."""
    try:
        return metal_monomer_codes(gemmi.read_structure(str(path)))
    except (RuntimeError, ValueError, OSError):
        return []
