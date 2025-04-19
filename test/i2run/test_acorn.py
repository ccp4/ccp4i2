import gemmi
from .utils import demoData, i2run


def test_acorn():
    args = ["acorn"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        gemmi.read_mtz_file(str(job / "PHSOUT.mtz"))
