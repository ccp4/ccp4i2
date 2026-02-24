import gemmi

from .utils import i2run


def test_comit(mtz8xfm):
    args = ["comit"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--F_PHI_IN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        gemmi.read_mtz_file(str(job / "F_PHI_OUT.mtz"))
