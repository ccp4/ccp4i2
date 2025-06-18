from pathlib import Path
import gemmi
from . import utils


def check_output(job: Path, code: str, hasLongLigandName=False):
    cifPath = job / f"{code}.cif"
    if not hasLongLigandName:
        pdbPath = job / f"{code}.pdb"
        gemmi.read_pdb(str(pdbPath))
    doc = gemmi.cif.read(str(cifPath))
    gemmi.make_chemcomp_from_block(doc[-1])


def test_from_mol():
    molUrl = "https://www.ebi.ac.uk/chebi/saveStructure.do"
    molUrl += "?defaultImage=true&chebiId=46195&imageId=0"
    with utils.download(molUrl) as molPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "MOL"]
        args += ["--TLC", "PCA"]
        args += ["--MOLIN", molPath]
        with utils.i2run(args) as job:
            check_output(job, "PCA")


def test_from_smiles():
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "SMILES"]
    args += ["--TLC", "HCA"]
    args += ["--SMILESIN", '"CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4"']
    with utils.i2run(args) as job:
        check_output(job, "HCA")


def test_from_sdf_rcsb(sdfA1LU6):
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "MOL"]
    args += ["--TLC", "A1LU6"]
    args += ["--MOLIN", sdfA1LU6]
    with utils.i2run(args) as job:
        check_output(job, "A1LU6", hasLongLigandName=True)


def test_from_cif_rcsb(cifA1LU6):
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "DICT"]
    args += ["--TLC", "A1LU6"]
    args += ["--DICTIN2", cifA1LU6]
    with utils.i2run(args) as job:
        check_output(job, "A1LU6", hasLongLigandName=True)


def test_from_cif_monomer_library():
    molUrl = "https://raw.githubusercontent.com/MonomerLibrary/"
    molUrl += "monomers/refs/heads/master/a/A1LU6.cif"
    with utils.download(molUrl) as molPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "DICT"]
        args += ["--TLC", "A1LU6"]
        args += ["--DICTIN2", molPath]
        args += ["--USE_COORD", "True"]
        with utils.i2run(args) as job:
            check_output(job, "A1LU6", hasLongLigandName=True)


def test_from_mol2():
    molUrl = "https://upjv.q4md-forcefieldtools.org/Tutorial/leap-off/Frag-DA.mol2"
    with utils.download(molUrl) as molPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "MOL2"]
        args += ["--TLC", "LIG"]
        args += ["--MOL2IN", molPath]
        with utils.i2run(args) as job:
            check_output(job, "LIG")


def test_from_cif_rcsb_metal_AF3(cifAF3, cif4dl8):
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "DICT"]
    args += ["--TLC", "AF3"]
    args += ["--DICTIN2", cifAF3]
    args += ["--TOGGLE_METAL", "True"]
    args += ["--METAL_STRUCTURE", cif4dl8]
    with utils.i2run(args) as job:
        check_output(job, "AF3")


def test_from_smiles_atom_name_matching():
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "SMILES"]
    args += ["--TLC", "LIG"]
    args += ["--SMILESIN", '"CO[C@@H]1[C@@H]([C@H](O[C@H]1N2C=NC3=C2N=C(NC3=O)N)COP(=O)(O)O)O"']
    args += ["--ATOMMATCHOPTION", "MONLIBCODE"]
    args += ["--MATCHTLC", "5GP"]
    args += ["--NRANDOM", "5"]
    with utils.i2run(args) as job:
        check_output(job, "LIG")


'''These tests crashes as acedrg or metalCoord crash themselves, not an i2 problem.

def test_from_cif_rcsb_metal_OEX(cifOEX, cif4ub6):
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "DICT"]
    args += ["--TLC", "OEX"]
    args += ["--DICTIN2", cifOEX]
    args += ["--TOGGLE_METAL", "True"]
    args += ["--METAL_STRUCTURE", cif4ub6]
    with utils.i2run(args) as job:
        check_output(job, "OEX")

def test_from_sdf_rcsb_metal_OEX(sdfOEX, cif4ub6):
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "MOL"]
    args += ["--TLC", "OEX"]
    args += ["--MOLIN", sdfOEX]
    args += ["--TOGGLE_METAL", "True"]
    args += ["--METAL_STRUCTURE", cif4ub6]
    with utils.i2run(args) as job:
        check_output(job, "OEX")
'''