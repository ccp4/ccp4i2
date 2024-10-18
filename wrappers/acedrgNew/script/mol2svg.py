import os,glob,re,time,sys
import math
import getopt

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def is_transition_metal(at):
    n = at.GetAtomicNum()
    return (n>=22 and n<=29) or (n>=40 and n<=47) or (n>=72 and n<=79)

def set_dative_bonds(mol, fromAtoms=(7,8)):
    """ convert some bonds to dative

    Replaces some single bonds between metals and atoms with atomic numbers in fomAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Returns the modified molecule.

    """
    pt = Chem.GetPeriodicTable()
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        for nbr in metal.GetNeighbors():
            if nbr.GetAtomicNum() in fromAtoms and \
               nbr.GetExplicitValence()>pt.GetDefaultValence(nbr.GetAtomicNum()) and \
               rwmol.GetBondBetweenAtoms(nbr.GetIdx(),metal.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                print("Changing a bond to dative",nbr.GetSymbol(),"->",metal.GetSymbol(),nbr.GetSymbol())
                rwmol.RemoveBond(nbr.GetIdx(),metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(),metal.GetIdx(),Chem.BondType.DATIVE)
    return rwmol

def svgFromMol(mol):
    m2 = set_dative_bonds(mol)
    try:
        Chem.SanitizeMol(m2draw2)
        Chem.Kekulize(m2draw2)
    except:
        print("Warning, could not sanitize/kekulize molecule.")

    d = rdMolDraw2D.MolDraw2DSVG(300, 300)
    d.DrawMolecule(m2)
    d.FinishDrawing()
    p = d.GetDrawingText()
    return p

if __name__ == "__main__":

    """
    #Some examples:
    python3 mol2svg.py --smiles="COc1cc(cc(c1O)OC)[C@@H]2c3cc4c(cc3[C@H]([C@@H]5[C@H]2C(=O)OC5)NC(=O)CC[C@@H]6C[NH2][Pt]([NH2]6)Cl)OCO4"
    python3 mol2svg.py --smiles=smiles="[Na+].[Na+].[O+]#[C][Fe-6]([C]#[O+])([C]#[O+])[C]#[O+]"
    python3 mol2svg.py --smiles=smiles="[BH2]1[H][BH2][H]1"
    """

    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:m:o:", ["smiles=", "molfile=", "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        sys.exit(2)
    verbose = False
    mol = None
    for o, a in opts:
        if o in ("-s", "--smiles"):
            mol = Chem.MolFromSmiles(a, sanitize=False)
        elif o in ("-m", "--molfile"):
            with open(a) as f:
               b = f.read()
               mol = Chem.MolFromMolBlock(b)
        elif o in ("-o", "--output"):
            output = a
        else:
            assert False, "unhandled option: "+o

    p = svgFromMol(mol)

    if output:
        with open(output,"w+") as f:
            f.write(p)
    else:
        print(p)
        
