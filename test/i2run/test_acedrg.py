from pathlib import Path
import gemmi
from . import utils


def check_output(directory: Path, code: str):
    cifPath = directory / f"{code}.cif"
    pdbPath = directory / f"{code}.pdb"
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
        args += ["--MOLIN", str(molPath)]
        with utils.i2run(args) as directory:
            check_output(directory, "PCA")


def test_from_smiles():
    args = ["LidiaAcedrgNew"]
    args += ["--MOLSMILESORSKETCH", "SMILES"]
    args += ["--TLC", "HCA"]
    args += ["--SMILESIN", "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4"]
    with utils.i2run(args) as directory:
        check_output(directory, "HCA")
