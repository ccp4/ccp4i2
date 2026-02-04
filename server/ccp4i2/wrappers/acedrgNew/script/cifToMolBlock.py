import sys
import tempfile
import logging
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from .setDativeBonds import set_dative_bonds

logger = logging.getLogger(__name__)


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
    logger.debug("No SMILES found, trying to extract atom coordinates")
    try:
        for block in doc:
            if found_molecule:
                break # We only take the first molecule.
            for search in searches:
                i = 0
                for thing in block.find_loop(search):
                    if not found_molecule:
                        found_molecule = True
                    if not i in atom_dict:
                        atom_dict[i] = {}
                    atom_dict[i][search] = thing
                    i += 1
    except Exception as e:
        logger.error("Failed to extract atoms from CIF: %s", e)
        return ""

    atom_list = list(atom_dict.values())
    logger.debug("Found %d atoms in CIF", len(atom_list))

    if len(atom_list) > 0:
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
            logger.debug("Successfully generated MolBlock from atom coordinates")
            return molBlock
        except Exception as e:
            logger.error("Failed to convert atoms to MolBlock: %s", e)
            return ""

    logger.warning("No atoms found in CIF file %s", input_file)
    return ""

if __name__ == "__main__":
    input_file = sys.argv[1]
    molBlock = cifFileToMolBlock(input_file)
    print(molBlock)
