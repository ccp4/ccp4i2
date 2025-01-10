"Testing that tasks are capable of using long ligand names in mmCIF using 8xfm"

from os import environ
from pathlib import Path
from random import choice
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import call
from tempfile import NamedTemporaryFile
from urllib.request import urlretrieve
import xml.etree.ElementTree as ET
import gemmi
import pytest


# Utilty functions and fixtures


def pdbefile(endpoint: str, suffix: str):
    "Download a file from PDBe to a temporary file and return the path"
    server = "https://www.ebi.ac.uk/pdbe"
    with NamedTemporaryFile(suffix=suffix) as tmp_file:
        urlretrieve(f"{server}/{endpoint}", tmp_file.name)
        yield tmp_file.name


@pytest.fixture(name="mmcif", scope="module")
def fixture_mmcif():
    "Download 8xfm.cif to a temporary file and return the path"
    yield from pdbefile("entry-files/download/8xfm.cif", "_8xfm.cif")


@pytest.fixture(name="sfcif", scope="module")
def fixture_sfcif():
    "Download r8xfmsf.ent to a temporary file and return the path"
    yield from pdbefile("entry-files/download/r8xfmsf.ent", "_r8xfmsf.ent")


@pytest.fixture(name="mtz", scope="module")
def fixture_mtz(sfcif):
    "Convert the structure factor CIF to MTZ format and return the path"
    with NamedTemporaryFile(suffix="_8xfm.mtz") as tmp_file:
        call(["gemmi", "cif2mtz", sfcif, tmp_file.name])
        yield tmp_file.name


@pytest.fixture(name="fasta", scope="module")
def fixture_fasta():
    "Download 8xfm.fasta to a temporary file and return the path"
    yield from pdbefile("entry/pdb/8xfm/fasta", "_8xfm.fasta")


def i2_path() -> Path:
    "Find the CCP4I2 installation"
    if "CCP4I2" in environ:
        return Path(environ["CCP4I2"])
    if "CCP4" in environ:
        ccp4 = Path(environ["CCP4"])
        for subdir in (
            "Frameworks/Python.framework/Versions/Current/lib/python3.7",
            "Frameworks/Python.framework/Versions/Current/lib/python3.9",
            "lib/python3.7",
            "lib/python3.9",
            "lib",
        ):
            path = ccp4 / subdir / "site-packages/ccp4i2"
            if path.exists():
                return path
    raise RuntimeError("CCP4I2 installation not found")


def has_residue_name(structure: gemmi.Structure, residue_name: str) -> bool:
    "Determine if any residue in the structure has a specific residue"
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == residue_name:
                    return True
    return False


def i2run(args, outputFilename="XYZOUT.cif"):
    "Run a task with i2run and check the output"
    chars = ascii_letters + digits
    tmp_name = "tmp_" + "".join(choice(chars) for _ in range(10))
    i2run_path = str(i2_path() / "bin" / "i2run")
    command = [i2run_path] + args
    command += ["--projectName", tmp_name]
    command += ["--projectPath", tmp_name]
    command += ["--dbFile", f"{tmp_name}.sqlite"]
    call(command)
    directory: Path = Path(tmp_name, "CCP4_JOBS", "job_1")
    xml_path = str(directory / "diagnostic.xml")
    out_path = str(directory / outputFilename)
    assert len(list(ET.parse(xml_path).iter("errorReport"))) == 0
    structure = gemmi.read_structure(out_path, format=gemmi.CoorFormat.Mmcif)
    assert has_residue_name(structure, "A1LU6")
    rmtree(tmp_name, ignore_errors=True)
    Path(f"{tmp_name}.sqlite").unlink(missing_ok=True)
    Path(f"{tmp_name}.sqlite-shm").unlink(missing_ok=True)
    Path(f"{tmp_name}.sqlite-wal").unlink(missing_ok=True)


# Tests


def test_8xfm_files(mmcif, sfcif, mtz, fasta):
    "Test that the downloaded file exists and can be read by gemmi"
    assert Path(mmcif).exists()
    assert Path(sfcif).exists()
    assert Path(mtz).exists()
    assert Path(fasta).exists()
    assert Path(mmcif).stat().st_size > 0
    assert Path(sfcif).stat().st_size > 0
    assert Path(mtz).stat().st_size > 0
    assert Path(fasta).stat().st_size > 0
    gemmi.cif.read(mmcif)
    gemmi.cif.read(sfcif)
    gemmi.read_mtz_file(mtz)
    gemmi.read_structure(mmcif)


# TODO
# def test_arp_warp_classic(mmcif, mtz, fasta):
#     "Test that arp_warp_classic can handle long ligand names in mmCIF"
#     args = ["arp_warp_classic"]
#     args += ["--AWA_MODELIN", mmcif]
#     args += ["--AWA_FOBS", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--AWA_FREE", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--AWA_SEQIN", f"seqFile={fasta}"]
#     args += ["--AWA_ARP_MODE", "WARPNTRACEMODEL"]


# TODO
# def test_buster(mmcif, mtz):
#     "Test that buster can handle long ligand names in mmCIF"
#     args = ["buster"]
#     args += ["--XYZIN", mmcif]
#     args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--WAT", "OFF"]


def test_coordinate_selector(mmcif):
    "Test that coordinate_selector can handle long ligand names in mmCIF"
    args = ["coordinate_selector"]
    args += ["--XYZIN", f"fullPath={mmcif}", "selection/text=(A1LU6)"]
    i2run(args)


def test_coot_find_waters(mmcif, mtz):
    "Test that coot_find_waters can handle long ligand names in mmCIF"
    args = ["coot_find_waters"]
    args += ["--XYZIN", mmcif]
    args += ["--FPHIIN", f"fullPath={mtz}", "columnLabels=/*/*/[FWT,PHWT]"]
    i2run(args)


def test_coot_rsr_morph(mmcif, mtz):
    "Test that coot_rsr_morph can handle long ligand names in mmCIF"
    args = ["coot_rsr_morph"]
    args += ["--XYZIN", mmcif]
    args += ["--FPHIIN", f"fullPath={mtz}", "columnLabels=/*/*/[FWT,PHWT]"]
    i2run(args)


def test_csymmatch(mmcif):
    "Test that csymmatch can handle long ligand names in mmCIF"
    args = ["csymmatch"]
    args += ["--XYZIN_QUERY", mmcif]
    args += ["--XYZIN_TARGET", mmcif]
    i2run(args)


# TODO
# def test_dr_mr_modelbuild_pipeline(mmcif, mtz, fasta):
#     "Test that dr_mr_modelbuild_pipeline can handle long ligand names in mmCIF"
#     args = ["dr_mr_modelbuild_pipeline"]
#     args += ["--XYZIN", mmcif]
#     args += ["--F_SIGF_IN", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREER_IN", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--ASUIN", f"seqFile={fasta}"]
#     args += ["--BUCC_NCYC", "5"]
#     args += ["--REFMAC_NCYC", "0"]
#     args += ["--BUCCANEER_OR_MODELCRAFT", "MODELCRAFT"]
#     args += ["--MERGED_OR_UNMERGED", "MERGED"]
#     args += ["--NMON", "1"]


# TODO
# def test_lorestr_i2(mmcif, mtz):
#     "Test that lorestr_i2 can handle long ligand names in mmCIF"
#     args = ["lorestr_i2"]
#     args += ["--XYZIN", f"fullPath={mmcif}"]
#     args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]


def test_modelcraft(mmcif, mtz, fasta):
    "Test that modelcraft can handle long ligand names in mmCIF"
    args = ["modelcraft"]
    args += ["--XYZIN", mmcif]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
    args += ["--ASUIN", f"seqFile={fasta}"]
    args += ["--CYCLES", "1"]
    i2run(args)


# TODO
# def test_molrep_pipe(mmcif, mtz, fasta):
#     "Test that molrep_pipe can handle long ligand names in mmCIF"
#     args = ["molrep_pipe"]
#     args += ["--XYZIN", mmcif]
#     args += ["--inputData.F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--inputData.FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--ASUIN", f"seqFile={fasta}"]
#     args += ["--RUNSHEETBEND", "False"]
#     args += ["--REFMAC_NCYC", "0"]
#     args += ["--NMON", "1"]
#     args += ["--SEQ", "n"]  # Don't modify search model


# TODO
# def test_phaser_rnp_pipeline(mmcif, mtz):
#     "Test that phaser_rnp_pipeline can handle long ligand names in mmCIF"
#     args = ["phaser_rnp_pipeline"]
#     args += ["--XYZIN_PARENT", mmcif]
#     args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]


# TODO
# def test_phaser_simple(mmcif, mtz, fasta):
#     "Test that phaser_simple can handle long ligand names in mmCIF"
#     args = ["phaser_simple"]
#     args += ["--XYZIN", mmcif]
#     args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--ASUFILE", f"seqFile={fasta}"]
#     args += ["--F_OR_I", "F"]
#     args += ["--COMP_BY", "ASU"]
#     args += ["--SEARCHSEQUENCEIDENTITY", "1.0"]


def test_prosmart_refmac(mmcif, mtz):
    "Test that prosmart_refmac can handle long ligand names in mmCIF"
    args = ["prosmart_refmac"]
    args += ["--XYZIN", mmcif]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
    args += ["--NCYCLES", "2"]
    args += ["--ADD_WATERS", "True"]
    i2run(args, "CIFFILE.pdb")


def test_servalcat_pipe(mmcif, mtz):
    "Test that servalcat can handle long ligand names in mmCIF"
    args = ["servalcat_pipe"]
    args += ["--XYZIN", mmcif]
    args += ["--DATA_METHOD", "xtal"]
    args += ["--HKLIN", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
    args += ["--NCYCLES", "2"]
    args += ["--F_SIGF_OR_I_SIGI", "F_SIGF"]
    i2run(args, "CIFFILE.pdb")


def test_sheetbend(mmcif, mtz):
    "Test that sheetbend can handle long ligand names in mmCIF"
    args = ["sheetbend"]
    args += ["--XYZIN", mmcif]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
    i2run(args)


# TODO
# def test_shelxemr(mmcif, mtz):
#     "Test that shelxeMR can handle long ligand names in mmCIF"
#     args = ["shelxeMR"]
#     args += ["--XYZIN", mmcif]
#     args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
#     args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
#     args += ["--FSOLVENT", "0.6"]


# TODO
# def test_substituteligand():
#     "Test that SubstituteLigand can handle long ligand names in mmCIF"
