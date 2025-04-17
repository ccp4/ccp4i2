"Testing that tasks are capable of using long ligand names in mmCIF using 8xfm"

from pathlib import Path
from subprocess import call
from tempfile import NamedTemporaryFile

import gemmi
import pytest

from . import utils


# Utility functions and fixtures


def pdbefile(endpoint: str):
    "Download a file from PDBe to a temporary file and return the path"
    server = "https://www.ebi.ac.uk/pdbe"
    with utils.download(f"{server}/{endpoint}") as path:
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
    with utils.i2run(args) as directory:
        out_path = str(directory / outputFilename)
        structure = gemmi.read_structure(out_path, format=gemmi.CoorFormat.Mmcif)
        assert has_residue_name(structure, "A1LU6")


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


def test_modelcraft(mmcif, mtz, fasta):
    "Test that modelcraft can handle long ligand names in mmCIF"
    args = ["modelcraft"]
    args += ["--XYZIN", mmcif]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FreeR_flag]"]
    args += ["--ASUIN", f"seqFile={fasta}"]
    args += ["--CYCLES", "1"]
    i2run(args)


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
