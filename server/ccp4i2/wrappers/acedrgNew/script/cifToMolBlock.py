import sys
import tempfile
import logging
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from .setDativeBonds import set_dative_bonds

logger = logging.getLogger(__name__)

# Mapping from CIF bond types to RDKit bond types
BOND_TYPE_MAP = {
    "sing": Chem.BondType.SINGLE,
    "doub": Chem.BondType.DOUBLE,
    "trip": Chem.BondType.TRIPLE,
    "arom": Chem.BondType.AROMATIC,
    "delo": Chem.BondType.ONEANDAHALF,  # Delocalized bonds
    "cova": Chem.BondType.SINGLE,  # Covalent treated as single
    "meta": Chem.BondType.SINGLE,  # Metal bonds treated as single (dative handled separately)
}


def cifBlockToMolBlock(block):
    """
    Convert a single CIF block to MolBlock format.

    Tries SMILES from _pdbx_chem_comp_descriptor first,
    falls back to atom coordinates from _chem_comp_atom.

    Returns empty string on failure.
    """
    # Try SMILES first
    molblock = _smiles_to_molblock(block)
    if molblock:
        return molblock

    # Fall back to atom/bond extraction
    molblock = _atoms_bonds_to_molblock(block)
    if molblock:
        return molblock

    return ""


def cifFileToMolBlocks(input_file):
    """
    Convert all monomers in a CIF dictionary file to MolBlock format.

    Iterates comp_* blocks (skipping comp_list), converts each to a molblock.
    Falls back to the last block if no comp_* blocks are found.

    Returns:
        Dict[str, str]: mapping monomer code -> molblock string.
        Entries with empty molblocks are omitted.
    """
    logger.debug("Converting all monomers in CIF to MolBlocks: %s", input_file)

    try:
        doc = gemmi.cif.read_file(input_file)
    except Exception as e:
        logger.error("Failed to read CIF file %s: %s", input_file, e)
        return {}

    logger.debug("Read CIF file with %d blocks", len(doc))
    result = {}

    for block in doc:
        if block.name == "comp_list":
            continue
        if block.name.startswith("comp_"):
            code = block.name[len("comp_"):]
            molblock = cifBlockToMolBlock(block)
            if molblock:
                result[code] = molblock

    # If no comp_* blocks found, try the last block as a fallback
    if not result and len(doc) > 0:
        block = doc[-1]
        if block.name != "comp_list":
            molblock = cifBlockToMolBlock(block)
            if molblock:
                code = block.name.replace("comp_", "") if block.name.startswith("comp_") else block.name
                result[code] = molblock

    logger.debug("Generated molblocks for %d monomers", len(result))
    return result


def cifFileToMolBlock(input_file):
    """
    Convert a ligand dictionary CIF file to MolBlock format.

    Returns the first monomer's molblock, or empty string on failure.
    """
    result = cifFileToMolBlocks(input_file)
    if result:
        return next(iter(result.values()))
    return ""


def _smiles_to_molblock(block):
    """
    Try to convert a CIF block to MolBlock via SMILES descriptors.

    Returns molblock string or empty string on failure.
    """
    smiles_dict = {}
    searches_SMILES = [
        "_pdbx_chem_comp_descriptor.type",
        "_pdbx_chem_comp_descriptor.descriptor",
    ]

    for search in searches_SMILES:
        i = 0
        for thing in block.find_loop(search):
            if i not in smiles_dict:
                smiles_dict[i] = {}
            smiles_dict[i][search] = thing
            i += 1

    smiles_list = list(smiles_dict.values())
    if not smiles_list:
        return ""

    logger.debug("Found %d SMILES descriptors in block %s", len(smiles_list), block.name)

    for s in smiles_list:
        if not s.get('_pdbx_chem_comp_descriptor.type', '').startswith("SMILES"):
            continue

        ss = s['_pdbx_chem_comp_descriptor.descriptor'].strip().lstrip('"').rstrip('"')
        logger.debug("Trying SMILES: %s", ss)

        mol = Chem.MolFromSmiles(ss)
        if not mol:
            logger.debug("Standard SMILES parsing failed, trying without sanitization")
            try:
                mol = Chem.MolFromSmiles(ss, sanitize=False)
                mol = set_dative_bonds(mol)
            except Exception as e:
                logger.debug("SMILES with dative bonds failed: %s", e)
                continue
            if not mol:
                logger.debug("SMILES parsing completely failed")
                continue

        try:
            AllChem.Compute2DCoords(mol)
            mol.SetProp("_MolFileChiralFlag", "1")
            molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
            logger.debug("Successfully generated MolBlock from SMILES")
            return molBlock
        except Exception as e:
            logger.debug("Failed to compute 2D coords: %s", e)
            continue

    return ""


def _atoms_bonds_to_molblock(block):
    """
    Try to convert a CIF block to MolBlock via atom coordinates and bonds.

    Returns molblock string or empty string on failure.
    """
    searches = [
        "_chem_comp_atom.comp_id",
        "_chem_comp_atom.atom_id",
        "_chem_comp_atom.type_symbol",
        "_chem_comp_atom.type_energy",
        "_chem_comp_atom.charge",
        "_chem_comp_atom.x",
        "_chem_comp_atom.y",
        "_chem_comp_atom.z",
        "_chem_comp_atom.model_Cartn_x",
        "_chem_comp_atom.model_Cartn_y",
        "_chem_comp_atom.model_Cartn_z",
    ]
    searches_bonds = [
        "_chem_comp_bond.atom_id_1",
        "_chem_comp_bond.atom_id_2",
        "_chem_comp_bond.type",
        "_chem_comp_bond.value_order",
        "_chem_comp_bond.aromatic",
        "_chem_comp_bond.pdbx_aromatic_flag",
    ]

    atom_dict = {}
    bond_dict = {}

    # Extract atoms
    for search in searches:
        i = 0
        for thing in block.find_loop(search):
            if i not in atom_dict:
                atom_dict[i] = {}
            atom_dict[i][search] = thing
            i += 1

    # Extract bonds
    for search in searches_bonds:
        i = 0
        for thing in block.find_loop(search):
            if i not in bond_dict:
                bond_dict[i] = {}
            bond_dict[i][search] = thing
            i += 1

    atom_list = list(atom_dict.values())
    bond_list = list(bond_dict.values())
    logger.debug("Found %d atoms and %d bonds in block %s", len(atom_list), len(bond_list), block.name)

    if not atom_list:
        return ""

    try:
        return _build_molblock_from_atoms_bonds(atom_list, bond_list)
    except Exception as e:
        logger.error("Failed to convert atoms/bonds to MolBlock: %s", e)
        logger.info("Falling back to PDB-based conversion (bond orders may be incorrect)")
        return _fallback_pdb_conversion(atom_list)


def _build_molblock_from_atoms_bonds(atom_list, bond_list):
    """Build an RDKit molecule from atom/bond lists and return its MolBlock."""
    rwmol = Chem.RWMol()
    conformer = Chem.Conformer(len(atom_list))

    # Map atom_id to RDKit atom index
    atom_id_to_idx = {}

    for i, atom in enumerate(atom_list):
        atom_id = atom['_chem_comp_atom.atom_id']
        element = atom['_chem_comp_atom.type_symbol']

        # Normalize element symbol to title case (e.g., "BR" -> "Br", "CL" -> "Cl")
        element = element.strip().capitalize()

        # Create RDKit atom
        rdatom = Chem.Atom(element)

        # Set formal charge if available
        try:
            charge = int(float(atom.get('_chem_comp_atom.charge', 0)))
            rdatom.SetFormalCharge(charge)
        except (ValueError, TypeError):
            pass

        # Add atom to molecule
        idx = rwmol.AddAtom(rdatom)
        atom_id_to_idx[atom_id] = idx

        # Set 3D coordinates if available
        x, y, z = 0.0, 0.0, 0.0
        if '_chem_comp_atom.x' in atom:
            x = float(atom['_chem_comp_atom.x'])
            y = float(atom['_chem_comp_atom.y'])
            z = float(atom['_chem_comp_atom.z'])
        elif '_chem_comp_atom.model_Cartn_x' in atom:
            x = float(atom['_chem_comp_atom.model_Cartn_x'])
            y = float(atom['_chem_comp_atom.model_Cartn_y'])
            z = float(atom['_chem_comp_atom.model_Cartn_z'])

        conformer.SetAtomPosition(idx, (x, y, z))

    # Add bonds with proper bond types
    aromatic_atoms = set()
    for bond in bond_list:
        atom1_id = bond.get('_chem_comp_bond.atom_id_1')
        atom2_id = bond.get('_chem_comp_bond.atom_id_2')

        if atom1_id not in atom_id_to_idx or atom2_id not in atom_id_to_idx:
            logger.warning("Bond references unknown atom: %s - %s", atom1_id, atom2_id)
            continue

        idx1 = atom_id_to_idx[atom1_id]
        idx2 = atom_id_to_idx[atom2_id]

        # Determine bond type
        bond_type_str = bond.get('_chem_comp_bond.type', bond.get('_chem_comp_bond.value_order', 'single'))
        bond_type_key = bond_type_str.lower()[:4] if bond_type_str else 'sing'

        # Check for aromatic flag
        is_aromatic = False
        aromatic_flag = bond.get('_chem_comp_bond.aromatic', bond.get('_chem_comp_bond.pdbx_aromatic_flag', ''))
        if aromatic_flag and aromatic_flag.lower() in ('y', 'yes', 'true', '1'):
            is_aromatic = True

        # Get RDKit bond type
        if is_aromatic or bond_type_key == 'arom':
            rdkit_bond_type = Chem.BondType.AROMATIC
            aromatic_atoms.add(idx1)
            aromatic_atoms.add(idx2)
        else:
            rdkit_bond_type = BOND_TYPE_MAP.get(bond_type_key, Chem.BondType.SINGLE)

        rwmol.AddBond(idx1, idx2, rdkit_bond_type)

    # Mark aromatic atoms
    for idx in aromatic_atoms:
        rwmol.GetAtomWithIdx(idx).SetIsAromatic(True)

    # Add conformer with 3D coordinates
    rwmol.AddConformer(conformer, assignId=True)

    # Try to sanitize the molecule
    mol = rwmol.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        logger.warning("Sanitization failed, trying partial sanitization: %s", e)
        try:
            # Try sanitizing without kekulization for aromatic systems
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        except Exception as e2:
            logger.warning("Partial sanitization also failed: %s", e2)
            # Continue anyway - may still produce usable output

    # Generate 2D coordinates for display
    AllChem.Compute2DCoords(mol)
    mol.SetProp("_MolFileChiralFlag", "1")
    molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
    logger.debug("Successfully generated MolBlock from atom coordinates and bonds")
    return molBlock


def _fallback_pdb_conversion(atom_list):
    """
    Fallback conversion using PDB format when direct RWMol construction fails.
    Note: This loses bond order information.
    """
    try:
        f = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        output_file = f.name
        f.close()

        s = gemmi.Structure()
        m = gemmi.Model("1")
        c = gemmi.Chain("A")
        r = gemmi.Residue()
        r.name = atom_list[0]['_chem_comp_atom.comp_id']
        sid = gemmi.SeqId("1")
        r.seqid = sid
        for atom in atom_list:
            a = gemmi.Atom()
            a.name = atom['_chem_comp_atom.atom_id']
            a.element = gemmi.Element(atom['_chem_comp_atom.type_symbol'])
            a.charge = int(float(atom['_chem_comp_atom.charge']))
            if '_chem_comp_atom.x' in atom:
                pos = gemmi.Position(float(atom['_chem_comp_atom.x']),float(atom['_chem_comp_atom.y']),float(atom['_chem_comp_atom.z']))
            elif '_chem_comp_atom.model_Cartn_x' in atom:
                pos = gemmi.Position(float(atom['_chem_comp_atom.model_Cartn_x']),float(atom['_chem_comp_atom.model_Cartn_y']),float(atom['_chem_comp_atom.model_Cartn_z']))
            else:
                logger.warning("No coordinates found for atom %s", atom.get('_chem_comp_atom.atom_id'))
                continue
            a.pos = pos
            r.add_atom(a)
        c.add_residue(r)
        m.add_chain(c)
        s.add_model(m)
        s.write_pdb(output_file)
        logger.debug("Wrote temporary PDB file: %s", output_file)
        mol = Chem.MolFromPDBFile(output_file)
        if mol is None:
            logger.error("RDKit failed to read PDB file")
            return ""
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_MolFileChiralFlag","1")
        molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
        logger.debug("Successfully generated MolBlock from PDB (bond orders may be incorrect)")
        return molBlock
    except Exception as e:
        logger.error("Fallback PDB conversion failed: %s", e)
        return ""


if __name__ == "__main__":
    input_file = sys.argv[1]
    molBlock = cifFileToMolBlock(input_file)
    print(molBlock)
