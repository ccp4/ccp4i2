from os import environ
from pathlib import Path
import gemmi
from .urls import rcsb_ligand_cif, rcsb_ligand_sdf, rcsb_mmcif
from .utils import download, i2run


def test_from_cif_monomer_library():
    cifPath = str(Path(environ["CLIBD_MON"], "a", "A1LU6.cif"))
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "DICT"]
    args += ["--TLC", "A1LU6"]
    args += ["--DICTIN2", cifPath]
    args += ["--USE_COORD", "True"]
    with i2run(args) as job:
        check_output(job, "A1LU6")


def test_from_cif_rcsb():
    with download(rcsb_ligand_cif("A1LU6")) as cifPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "DICT"]
        args += ["--TLC", "A1LU6"]
        args += ["--DICTIN2", cifPath]
        with i2run(args) as job:
            check_output(job, "A1LU6")


def test_from_cif_rcsb_metal_AF3():
    with download(rcsb_ligand_cif("AF3")) as ligand, download(rcsb_mmcif("4dl8")) as mmcif:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "DICT"]
        args += ["--TLC", "AF3"]
        args += ["--DICTIN2", ligand]
        args += ["--TOGGLE_METAL", "True"]
        args += ["--METAL_STRUCTURE", mmcif]
        with i2run(args) as job:
            check_output(job, "AF3")


def test_from_mol():
    molUrl = "https://www.ebi.ac.uk/chebi/backend/api/public/molfile/46195"
    with download(molUrl) as molPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "MOL"]
        args += ["--TLC", "PCA"]
        args += ["--MOLIN", molPath]
        with i2run(args) as job:
            check_output(job, "PCA")


def test_from_sdf():
    with download(rcsb_ligand_sdf("A1LU6")) as sdfPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "MOL"]
        args += ["--TLC", "A1LU6"]
        args += ["--MOLIN", sdfPath]
        with i2run(args) as job:
            check_output(job, "A1LU6")


def test_from_smiles():
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "SMILES"]
    args += ["--TLC", "HCA"]
    args += ["--SMILESIN", '"CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4"']
    with i2run(args) as job:
        check_output(job, "HCA")


def test_from_mol2():
    url = "https://upjv.q4md-forcefieldtools.org/Tutorial/leap-off/Frag-DA.mol2"
    with download(url) as molPath:
        args = ["LidiaAcedrgNew"]
        args += ["--MOLSMILESORSKETCH", "MOL2"]
        args += ["--TLC", "LIG"]
        args += ["--MOL2IN", molPath]
        with i2run(args) as job:
            check_output(job, "LIG")


def test_from_smiles_atom_name_matching():
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "SMILES"]
    args += ["--TLC", "LIG"]
    args += ["--SMILESIN", '"CO[C@@H]1[C@@H]([C@H](O[C@H]1N2C=NC3=C2N=C(NC3=O)N)COP(=O)(O)O)O"']
    args += ["--ATOMMATCHOPTION", "MONLIBCODE"]
    args += ["--MATCHTLC", "5GP"]
    args += ["--NRANDOM", "5"]
    with i2run(args) as job:
        check_output(job, "LIG")


def check_output(job: Path, code: str):
    if len(code) <= 3:
        gemmi.read_pdb(str(job / f"{code}.pdb"))
    doc = gemmi.cif.read(str(job / f"{code}.cif"))
    gemmi.make_chemcomp_from_block(doc[-1])


# These tests crashes as acedrg or metalCoord crash themselves, not an i2 problem.

# @fixture(name="cif4ub6", scope="module")
# def cif4ub6_fixture():
#     with download(rcsb_mmcif("4ub6")) as path:
#         yield path

# def test_from_cif_rcsb_metal_OEX(cif4ub6):
#     with download(rcsb_ligand_cif("OEX")) as cifOEX:
#         args = ["LidiaAcedrgNew"]
#         args += ["--MOLSMILESORSKETCH", "DICT"]
#         args += ["--TLC", "OEX"]
#         args += ["--DICTIN2", cifOEX]
#         args += ["--TOGGLE_METAL", "True"]
#         args += ["--METAL_STRUCTURE", cif4ub6]
#         with i2run(args) as job:
#             check_output(job, "OEX")

# def test_from_sdf_rcsb_metal_OEX(cif4ub6):
#     with download(rcsb_ligand_sdf("OEX")) as sdfOEX:
#         args = ["LidiaAcedrgNew"]
#         args += ["--MOLSMILESORSKETCH", "MOL"]
#         args += ["--TLC", "OEX"]
#         args += ["--MOLIN", sdfOEX]
#         args += ["--TOGGLE_METAL", "True"]
#         args += ["--METAL_STRUCTURE", cif4ub6]
#         with i2run(args) as job:
#             check_output(job, "OEX")
