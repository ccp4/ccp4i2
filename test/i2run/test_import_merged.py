from gemmi import read_mtz_file
from .urls import pdbe_sfcif
from .utils import demoData, download, i2run


def test_2ceu_cif():
    with download(pdbe_sfcif("2ceu")) as cif:
        args = ["import_merged"]
        args += ["--HKLIN", str(cif)]
        args += ["--SPACEGROUP", "I 2 2 2"]
        with i2run(args) as job:
            for name in ["OBS", "FREE"]:
                read_mtz_file(str(job / f"{name}OUT.mtz"))
            assert False, "Check output files"


def test_gamma_mtz():
    args = ["import_merged"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    with i2run(args) as job:
        for name in ["OBS", "FREE"]:
            read_mtz_file(str(job / f"{name}OUT.mtz"))
