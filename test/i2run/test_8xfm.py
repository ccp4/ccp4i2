"Testing that tasks are capable of using long ligand names in mmCIF using 8xfm"

from contextlib import contextmanager
from email.message import EmailMessage
from os.path import basename
from pathlib import Path
from random import choice
from shutil import rmtree
from string import ascii_letters, digits
from subprocess import call
from multiprocessing import Process
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
import xml.etree.ElementTree as ET

from requests import get, Response
import gemmi
import pytest

from ...core import CCP4I2Runner


# Utility functions and fixtures


def valid_filename_from_response(response: Response):
    """
    Extracts a valid filename from the response headers or URL.
    Ensures that the filename is safe to use by stripping whitespace,
    replacing spaces with underscores, and removing any characters
    that are not alphanumeric, dash, underscore or dot.
    """
    message = EmailMessage()
    for header, value in response.headers.items():
        message[header] = value
    url_name = unquote(basename(urlparse(response.url).path))
    name = message.get_filename() or url_name
    name = name.strip().replace(" ", "_")
    name = "".join(c for c in name if c.isalnum() or c in "-_.")
    return name


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a pathlib.Path object to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    response = get(url, allow_redirects=True, stream=True, timeout=30)
    response.raise_for_status()
    suffix = f"_{valid_filename_from_response(response)}"
    with NamedTemporaryFile(suffix=suffix, delete=False) as temp:
        for chunk in response.iter_content(chunk_size=1_000_000):
            temp.write(chunk)
    path = Path(temp.name).resolve()
    try:
        yield path
    finally:
        path.unlink(missing_ok=True)


def pdbefile(endpoint: str):
    "Download a file from PDBe to a temporary file and return the path"
    server = "https://www.ebi.ac.uk/pdbe"
    with download(f"{server}/{endpoint}") as path:
        yield str(path)


@pytest.fixture(name="mmcif", scope="module")
def fixture_mmcif():
    "Download 8xfm.cif to a temporary file and return the path"
    yield from pdbefile("entry-files/download/8xfm.cif")


@pytest.fixture(name="sfcif", scope="module")
def fixture_sfcif():
    "Download r8xfmsf.ent to a temporary file and return the path"
    yield from pdbefile("entry-files/download/r8xfmsf.ent")


@pytest.fixture(name="mtz", scope="module")
def fixture_mtz(sfcif):
    "Convert the structure factor CIF to MTZ format and return the path"
    temp = NamedTemporaryFile(suffix="_8xfm.mtz", delete=False)
    temp.close()
    call(["gemmi", "cif2mtz", sfcif, temp.name])
    try:
        yield temp.name
    finally:
        Path(temp.name).unlink(missing_ok=True)


@pytest.fixture(name="fasta", scope="module")
def fixture_fasta():
    "Download 8xfm.fasta to a temporary file and return the path"
    yield from pdbefile("entry/pdb/8xfm/fasta")


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
    args = ["i2run"] + args
    args += ["--projectName", tmp_name]
    args += ["--projectPath", tmp_name]
    args += ["--dbFile", f"{tmp_name}.sqlite"]
    process = Process(target=CCP4I2Runner.main, args=(args,))
    process.start()
    process.join()
    directory: Path = Path(tmp_name, "CCP4_JOBS", "job_1")
    xml_path = str(directory / "diagnostic.xml")
    out_path = str(directory / outputFilename)
    assert len(list(ET.parse(xml_path).iter("errorReport"))) == 0
    structure = gemmi.read_structure(out_path, format=gemmi.CoorFormat.Mmcif)
    assert has_residue_name(structure, "A1LU6")
    rmtree(tmp_name, ignore_errors=True)
    for extension in ("sqlite", "sqlite-shm", "sqlite-wal"):
        Path(f"{tmp_name}.{extension}").unlink(missing_ok=True)


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
    args += ["--VALIDATE_IRIS", "False"]
    args += ["--VALIDATE_BAVERAGE", "False"]
    args += ["--VALIDATE_RAMACHANDRAN", "False"]
    args += ["--VALIDATE_MOLPROBITY", "False"]
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
