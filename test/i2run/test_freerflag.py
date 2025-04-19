from pytest import approx
import gemmi
from .utils import demoData, i2run


def test_freerflag():
    args = ["freerflag"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--COMPLETE", "False"]
    args += ["--FRAC", "0.2"]
    with i2run(args) as job:
        mtz = gemmi.read_mtz_file(str(job / "FREEROUT.mtz"))
        column = mtz.rfree_column()
        fraction = list(column).count(0) / len(column)
        assert fraction == approx(0.2, abs=0.01)
