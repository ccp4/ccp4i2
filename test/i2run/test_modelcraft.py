from pathlib import Path
import json
from gemmi import CoorFormat, read_mtz_file, read_structure
from .utils import demoData, hasLongLigandName, i2run


def _check_output(job: Path):
    read_structure(str(job / "XYZOUT.cif"), format=CoorFormat.Mmcif)
    for name in ["ABCD", "DIFFPHI", "FPHI"]:
        read_mtz_file(str(job / f"{name}OUT.mtz"))
    with (job / "modelcraft" / "modelcraft.json").open() as json_file:
        results = json.load(json_file)
        assert results["final"]["r_free"] < 0.3


def test_8xfm(cif8xfm, mtz8xfm, seq8xfm):
    args = ["modelcraft"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--ASUIN", f"seqFile={seq8xfm}"]
    args += ["--XYZIN", cif8xfm]
    args += ["--CYCLES", "2"]
    with i2run(args) as job:
        _check_output(job)
        assert hasLongLigandName(job / "XYZOUT.cif")


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
