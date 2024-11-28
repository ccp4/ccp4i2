import sys
import tempfile
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from .setDativeBonds import set_dative_bonds


def cifFileToMolBlock(input_file):

    f = tempfile.NamedTemporaryFile(suffix=".pdb")
    output_file = f.name
    f.close()

    searches_SMILES = ["_pdbx_chem_comp_descriptor.type","_pdbx_chem_comp_descriptor.descriptor"]

    bond_searches = ["_chem_comp_atom.comp_id",
    "_chem_comp_bond.atom_id_1",
    "_chem_comp_bond.atom_id_2",
    "_chem_comp_bond.value_order"
    ]

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
    bond_dict = {}
    found_molecule = False
    found_SMILES = False

    smiles_dict = {}

#First of all, we see if there are any SMILES strings in this file and take the first one if there are.
    try:
        doc = gemmi.cif.read_file(input_file)  # copy all the data from mmCIF file
        for block in doc:
            for search in searches_SMILES:
                i = 0
                for thing in block.find_loop(search):
                    if not i in smiles_dict:
                        smiles_dict[i] = {}
                    smiles_dict[i][search] = thing
                    i += 1
    except Exception as e:
        print("Oops. %s" % e)
        return ""

    smiles_list = list(smiles_dict.values())
    mol = None
    for s in smiles_list:
        if s['_pdbx_chem_comp_descriptor.type'].startswith("SMILES"):
            ss = s['_pdbx_chem_comp_descriptor.descriptor'].strip().lstrip('"').rstrip('"')
            mol = Chem.MolFromSmiles(ss)
            if not mol:
                print("Not successful Chem.MolFromSmiles(ss) for ss =", ss)
                try:
                    mol = Chem.MolFromSmiles(ss, sanitize=False)
                    mol = set_dative_bonds(mol)
                except:
                    print("Not successful try Chem.MolFromSmiles(ss, sanitize=False) or set_dative_bonds(mol) for ss =", ss)
                    continue
                if not mol:
                    print("Not successful Chem.MolFromSmiles(ss, sanitize=False) or set_dative_bonds(mol) for ss =", ss)
                    continue
            if mol:
                try:
                    AllChem.Compute2DCoords(mol)
                    mol.SetProp("_MolFileChiralFlag","1")
                    molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
                    return molBlock
                except:
                    print("Not successful try set_dative_bonds(mol) for ", ss)
                    continue
            else:  # invalid SMILES
                print("Not successful set_dative_bonds(mol) for", ss)
                continue

#If no SMILES string, we hope there is at least one molecule, either in PDB or acedrg chem_comp format.
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
            for search in bond_searches:
                i = 0
                for thing in block.find_loop(search):
                    if not i in bond_dict:
                        bond_dict[i] = {}
                    bond_dict[i][search] = thing
                    i += 1
    except Exception as e:
        print("Oops. %s" % e)
        return ""

    atom_list = list(atom_dict.values())
    bond_list = list(bond_dict.values())

    if len(atom_list) > 0:
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
            a.pos = pos
            r.add_atom(a)
        c.add_residue(r)
        m.add_chain(c)
        s.add_model(m)
        s.write_pdb(output_file)
        mol = Chem.MolFromPDBFile(output_file)
        for b in mol.GetBonds():
            at1 = atom_list[b.GetBeginAtomIdx()]['_chem_comp_atom.atom_id']
            at2 = atom_list[b.GetEndAtomIdx()]['_chem_comp_atom.atom_id']
            for b2 in bond_list:
                if '_chem_comp_bond.atom_id_1' in b2 and '_chem_comp_bond.atom_id_2' in b2:
                    if (b2['_chem_comp_bond.atom_id_1'] == at1 and b2['_chem_comp_bond.atom_id_2'] == at2) or (b2['_chem_comp_bond.atom_id_1'] == at2 and b2['_chem_comp_bond.atom_id_2'] == at1):
                        if b2["_chem_comp_bond.value_order"].upper()[0:4] == "AROM":
                            b.SetBondType(Chem.BondType.AROMATIC)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:4] == "SING":
                            b.SetBondType(Chem.BondType.SINGLE)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:4] == "DOUB":
                            b.SetBondType(Chem.BondType.DOUBLE)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:4] == "TRIP":
                            b.SetBondType(Chem.BondType.TRIPLE)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:4] == "QUAD":
                            b.SetBondType(Chem.BondType.QUADRUPLE)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:4] == "QUIN":
                            b.SetBondType(Chem.BondType.QUINTUPLE)
                        elif b2["_chem_comp_bond.value_order"].upper()[0:3] == "HEX":
                            b.SetBondType(Chem.BondType.HEXTUPLE)
                        elif b2["_chem_comp_bond.value_order"].upper() == "DATIVE":
                            b.SetBondType(Chem.BondType.DATIVE)
                        elif b2["_chem_comp_bond.value_order"].upper() == "DATIVEL":
                            b.SetBondType(Chem.BondType.DATIVEL)
                        elif b2["_chem_comp_bond.value_order"].upper() == "DATIVER":
                            b.SetBondType(Chem.BondType.DATIVER)
                        elif b2["_chem_comp_bond.value_order"].upper() == "DATIVEONE":
                            b.SetBondType(Chem.BondType.DATIVEONE)
                        print("Match!",b.GetBondType(),b2["_chem_comp_bond.value_order"])
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_MolFileChiralFlag","1")
        molBlock = Chem.MolToMolBlock(mol, includeStereo=True, forceV3000=False)
        return molBlock

    return ""

if __name__ == "__main__":
    input_file = sys.argv[1]
    molBlock = cifFileToMolBlock(input_file)
    print(molBlock)
