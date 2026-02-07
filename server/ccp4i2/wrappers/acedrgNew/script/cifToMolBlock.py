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


def cifFileToMolBlock(input_file):
    """
    Convert a ligand dictionary CIF file to MolBlock format.

    Tries two approaches:
    1. Extract SMILES from _pdbx_chem_comp_descriptor and convert to 2D
    2. Extract atom coordinates from _chem_comp_atom and convert via PDB

    Returns empty string on failure.
    """
    logger.debug("Converting CIF to MolBlock: %s", input_file)

    f = tempfile.NamedTemporaryFile(suffix=".pdb")
    output_file = f.name
    f.close()

    searches_SMILES = ["_pdbx_chem_comp_descriptor.type","_pdbx_chem_comp_descriptor.descriptor"]

    searches = ["_chem_comp_atom.comp_id",
    "_chem_comp_atom.atom_id",
    "_chem_comp_atom.type_symbol",
    "_chem_comp_atom.type_energy",
    "_chem_comp_atom.charge",
    "_chem_comp_atom.x",
    "_chem_comp_atom.y",
    "_chem_comp_atom.z",
    "_chem_comp_atom.model_Cartn_x",
    "_chem_comp_atom.model_Cartn_y",
    "_chem_comp_atom.model_Cartn_z" ]

    searches_bonds = [
    "_chem_comp_bond.atom_id_1",
    "_chem_comp_bond.atom_id_2",
    "_chem_comp_bond.type",
    "_chem_comp_bond.value_order",  # Alternative field name
    "_chem_comp_bond.aromatic",
    "_chem_comp_bond.pdbx_aromatic_flag",  # Alternative field name
    ]

    atom_dict = {}
    found_molecule = False
    found_SMILES = False

    smiles_dict = {}

#First of all, we see if there are any SMILES strings in this file and take the first one if there are.
    try:
        doc = gemmi.cif.read_file(input_file)  # copy all the data from mmCIF file
        logger.debug("Read CIF file with %d blocks", len(doc))
        for block in doc:
            logger.debug("Processing block: %s", block.name)
            for search in searches_SMILES:
                i = 0
                for thing in block.find_loop(search):
                    if not i in smiles_dict:
                        smiles_dict[i] = {}
                    smiles_dict[i][search] = thing
                    i += 1
    except Exception as e:
        logger.error("Failed to read CIF file %s: %s", input_file, e)
        return ""

    smiles_list = list(smiles_dict.values())
    logger.debug("Found %d SMILES descriptors", len(smiles_list))
    mol = None
    for s in smiles_list:
        if s['_pdbx_chem_comp_descriptor.type'].startswith("SMILES"):
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
            if mol:
                try:
                    AllChem.Compute2DCoords(mol)
                    mol.SetProp("_MolFileChiralFlag","1")
                    molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                    logger.debug("Successfully generated MolBlock from SMILES")
                    return molBlock
                except Exception as e:
                    logger.debug("Failed to compute 2D coords: %s", e)
                    continue
            else:  # invalid SMILES
                logger.debug("Invalid SMILES: %s", ss)
                continue

#If no SMILES string, we hope there is at least one molecule, either in PDB or acedrg chem_comp format.
    logger.debug("No SMILES found, trying to extract atom coordinates and bonds")
    bond_dict = {}
    try:
        for block in doc:
            if found_molecule:
                break # We only take the first molecule.
            # Extract atoms
            for search in searches:
                i = 0
                for thing in block.find_loop(search):
                    if not found_molecule:
                        found_molecule = True
                    if not i in atom_dict:
                        atom_dict[i] = {}
                    atom_dict[i][search] = thing
                    i += 1
            # Extract bonds
            for search in searches_bonds:
                i = 0
                for thing in block.find_loop(search):
                    if not i in bond_dict:
                        bond_dict[i] = {}
                    bond_dict[i][search] = thing
                    i += 1
    except Exception as e:
        logger.error("Failed to extract atoms/bonds from CIF: %s", e)
        return ""

    atom_list = list(atom_dict.values())
    bond_list = list(bond_dict.values())
    logger.debug("Found %d atoms and %d bonds in CIF", len(atom_list), len(bond_list))

    if len(atom_list) > 0:
        try:
            # Build molecule using RWMol with explicit atoms and bonds
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

        except Exception as e:
            logger.error("Failed to convert atoms/bonds to MolBlock: %s", e)
            # Fall back to PDB-based approach (loses bond orders but may still work)
            logger.info("Falling back to PDB-based conversion (bond orders may be incorrect)")
            return _fallback_pdb_conversion(atom_list, output_file)

    logger.warning("No atoms found in CIF file %s", input_file)
    return ""


def _fallback_pdb_conversion(atom_list, output_file):
    """
    Fallback conversion using PDB format when direct RWMol construction fails.
    Note: This loses bond order information.
    """
    try:
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
