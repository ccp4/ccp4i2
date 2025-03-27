from pathlib import Path
from gemmi import CoorFormat, read_mtz_file, read_structure
from .utils import demoData, i2run


def _check_output(job: Path):
    read_structure(str(job / "XYZOUT.cif"), format=CoorFormat.Mmcif)
    for name in ["ABCD", "DIFFPHI", "FPHI"]:
        read_mtz_file(str(job / f"{name}OUT.mtz"))


def test_gamma_ep():
    args = ["modelcraft"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--PHASES", demoData("gamma", "initial_phases.mtz")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--USE_MODEL_PHASES", "False"]
    args += ["--CYCLES", "2"]
    with i2run(args) as job:
        _check_output(job)


def test_gamma_mr():
    args = ["modelcraft"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--CYCLES", "2"]
    with i2run(args) as job:
        _check_output(job)
