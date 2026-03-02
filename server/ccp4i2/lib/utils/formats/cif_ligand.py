import logging
import gemmi
from typing import List, Dict, Optional, Any

logger = logging.getLogger(__name__)


def parse_cif_ligand_summary(cif_file_path: str) -> List[Dict[str, Optional[Any]]]:
    """
    Parse a CIF file and return a summary of ligand information.

    Args:
        cif_file_path (str): Path to the CIF file

    Returns:
        List[Dict[str, Optional[Any]]]: List of summary dictionaries with ligand information
    """
    summary_list = []

    try:
        # Read the CIF file
        doc = gemmi.cif.read_file(cif_file_path)

        if len(doc) == 0:
            return summary_list

        # Look for comp_list block first
        comp_list_block = _find_comp_list_block(doc)

        if comp_list_block:
            summary_list.extend(_extract_from_comp_list(comp_list_block))
        else:
            # Fallback: look for individual comp blocks (like comp_LIG)
            summary_list.extend(_extract_from_individual_comp_blocks(doc))

    except FileNotFoundError:
        raise FileNotFoundError(f"CIF file not found: {cif_file_path}")
    except gemmi.cif.ParsingError as e:
        raise ValueError(f"Error parsing CIF file {cif_file_path}: {e}")
    except Exception as e:
        raise RuntimeError(f"Unexpected error processing CIF file {cif_file_path}: {e}")

    return summary_list


def _find_comp_list_block(doc: gemmi.cif.Document) -> Optional[gemmi.cif.Block]:
    """Find the comp_list block in the CIF document."""
    for block in doc:
        if block.name == "comp_list":
            return block
    return None


def _extract_from_comp_list(
    comp_list_block: gemmi.cif.Block,
) -> List[Dict[str, Optional[Any]]]:
    """Extract ligand summaries from comp_list block."""
    summaries = []

    # Check if the required tags exist
    if "_chem_comp.id" not in comp_list_block.as_string():
        return summaries

    try:
        comp_table = comp_list_block.find(
            "_chem_comp.",
            [
                "id",
                "three_letter_code",
                "name",
                "group",
                "number_atoms_all",
                "number_atoms_nh",
                "desc_level",
            ],
        )

        if comp_table:
            for row in comp_table:
                summary = _create_summary_from_row(row)
                summaries.append(summary)

    except Exception as e:
        print(f"Warning: Error extracting from comp_list block: {e}")

    return summaries


def _extract_from_individual_comp_blocks(
    doc: gemmi.cif.Document,
) -> List[Dict[str, Optional[Any]]]:
    """Extract ligand summaries from individual comp_* blocks."""
    summaries = []

    for block in doc:
        if block.name.startswith("comp_"):
            # Extract ligand ID from block name
            ligand_id = block.name.replace("comp_", "")
            summary = _create_empty_summary()
            summary["id"] = ligand_id
            summary["three_letter_code"] = ligand_id
            summaries.append(summary)

    return summaries


def _create_summary_from_row(row: Dict[str, str]) -> Dict[str, Optional[Any]]:
    """Create a summary dictionary from a comp_table row."""
    summary = _create_empty_summary()

    # Extract basic string fields
    summary["id"] = _safe_get_value(row, "id")
    summary["three_letter_code"] = _safe_get_value(row, "three_letter_code")
    summary["name"] = _clean_quoted_string(_safe_get_value(row, "name"))
    summary["group"] = _safe_get_value(row, "group")
    summary["desc_level"] = _safe_get_value(row, "desc_level")

    # Extract and convert numeric fields
    summary["number_atoms_all"] = _safe_convert_to_int(row, "number_atoms_all")
    summary["number_atoms_nh"] = _safe_convert_to_int(row, "number_atoms_nh")

    return summary


def extract_monomer_atoms_bonds(cif_file_path: str) -> Dict[str, Any]:
    """
    Extract non-hydrogen atom names and bond connectivity from a CIF dictionary.

    Works with both CCP4 monomer library files and custom AceDRG dictionaries.
    Returns data for the first (or only) comp block found.

    Returns:
        Dict with "atoms" (list of atom ID strings) and
        "bonds" (list of {"atom1", "atom2", "type"} dicts).
    """
    try:
        doc = gemmi.cif.read_file(cif_file_path)
    except Exception:
        return {"atoms": [], "bonds": []}

    # Find the comp block — try comp_* blocks, fall back to last block
    block = None
    for b in doc:
        if b.name.startswith("comp_") and b.name != "comp_list":
            block = b
            break
    if block is None and len(doc) > 0:
        block = doc[-1]
    if block is None:
        return {"atoms": [], "bonds": []}

    return _extract_atoms_bonds_from_block(block)


def extract_all_monomers_atoms_bonds(cif_file_path: str) -> Dict[str, Dict[str, Any]]:
    """
    Extract non-hydrogen atom names and bond connectivity for ALL monomers
    in a CIF dictionary file.

    A dictionary file may contain multiple monomer definitions (comp_* blocks).
    This function returns a dict keyed by monomer code, where each value
    contains the atoms and bonds for that monomer.

    Returns:
        Dict mapping monomer code -> {"atoms": [...], "bonds": [...]}.
        If a comp_list block is present, the code is taken from there;
        otherwise it is derived from the block name (comp_CODE -> CODE).
    """
    try:
        doc = gemmi.cif.read_file(cif_file_path)
    except Exception:
        return {}

    monomers: Dict[str, Dict[str, Any]] = {}

    for block in doc:
        if block.name == "comp_list":
            continue
        if block.name.startswith("comp_"):
            code = block.name[len("comp_"):]
            result = _extract_atoms_bonds_from_block(block)
            if result["atoms"]:  # Only include monomers with actual atom data
                monomers[code] = result

    # If no comp_* blocks found, try the last block as a fallback
    if not monomers and len(doc) > 0:
        block = doc[-1]
        if block.name != "comp_list":
            result = _extract_atoms_bonds_from_block(block)
            if result["atoms"]:
                # Use block name as code, or "UNK" as last resort
                code = block.name.replace("comp_", "") if block.name.startswith("comp_") else block.name
                monomers[code] = result

    return monomers


def _extract_atoms_bonds_from_block(block: gemmi.cif.Block) -> Dict[str, Any]:
    """Extract non-hydrogen atoms and bonds from a single CIF comp block."""
    try:
        comp = gemmi.make_chemcomp_from_block(block)
    except Exception:
        return {"atoms": [], "bonds": []}

    # Non-hydrogen atoms
    atoms = [a.id for a in comp.atoms if a.el.name != "H"]
    atom_set = set(atoms)

    # Bonds between non-hydrogen atoms
    bond_type_map = {
        "Single": "single",
        "Double": "double",
        "Triple": "triple",
        "Aromatic": "aromatic",
        "Deloc": "deloc",
        "Metal": "metal",
    }
    bonds = []
    for b in comp.rt.bonds:
        a1, a2 = b.id1.atom, b.id2.atom
        if a1 in atom_set and a2 in atom_set:
            bonds.append({
                "atom1": a1,
                "atom2": a2,
                "type": bond_type_map.get(b.type.name, b.type.name),
            })

    return {"atoms": atoms, "bonds": bonds}


def generate_all_molblocks(cif_file_path: str) -> Dict[str, str]:
    """
    Generate 2D MolBlock strings for all monomers in a CIF dictionary file.

    Uses the cifFileToMolBlocks function which tries SMILES-first conversion,
    falling back to atom coordinate extraction for each monomer block.

    Returns:
        Dict mapping monomer code -> MolBlock V2000 string.
        Monomers that fail conversion are omitted.
    """
    try:
        from ....wrappers.acedrgNew.script.cifToMolBlock import cifFileToMolBlocks
        return cifFileToMolBlocks(cif_file_path)
    except ImportError as e:
        logger.warning("Cannot generate molblocks (missing dependency): %s", e)
        return {}
    except Exception as e:
        logger.warning("Failed to generate molblocks for %s: %s", cif_file_path, e)
        return {}


def _create_empty_summary() -> Dict[str, Optional[Any]]:
    """Create an empty summary dictionary with the expected structure."""
    return {
        "id": None,
        "three_letter_code": None,
        "name": None,
        "group": None,
        "number_atoms_all": None,
        "number_atoms_nh": None,
        "desc_level": None,
    }


def _safe_get_value(row: Dict[str, str], key: str) -> Optional[str]:
    """Safely get a value from the row dictionary."""
    value = row[key]
    return None if value in (None, "", ".") else value


def _clean_quoted_string(value: Optional[str]) -> Optional[str]:
    """Remove surrounding quotes from string values."""
    if value is None:
        return None
    return value.strip("'\"")


def _safe_convert_to_int(row: Dict[str, str], key: str) -> Optional[int]:
    """Safely convert a row value to integer."""
    value = _safe_get_value(row, key)
    if value is None:
        return None

    try:
        return int(value)
    except (ValueError, TypeError):
        return None


# Example usage and testing
if __name__ == "__main__":
    # Test with a sample file path
    try:
        result = parse_cif_ligand_summary("sample_ligand.cif")
        print("Parsed CIF summary:")
        for i, summary in enumerate(result):
            print(f"Ligand {i + 1}: {summary}")
    except Exception as e:
        print(f"Error: {e}")
