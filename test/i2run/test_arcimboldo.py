from sys import platform
from pytest import mark
from .utils import demoData, i2run


pytestmark = mark.skipif(platform == "win32", reason="Not supported on Windows")


def test_arcimboldo():
    args = ["arcimboldo"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ARCIMBOLDO_OPTIONS", "LITE"]
    args += ["--MOLECULAR_WEIGHT", "15000"]
    args += ["--N_FRAGMENTS", "1"]
    with i2run(args) as job:
        assert False, "Check output files"
