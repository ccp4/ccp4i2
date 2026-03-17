import pytest
from .utils import i2run
import gemmi


@pytest.mark.skip(reason="Zanuda's internal refmac5 fails to parse 8xfm mmCIF model")
def test_8xfm(cif8xfm, mtz8xfm):
    """Test zanuda space group validation with 8xfm."""
    args = ["zanuda"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--XYZIN", cif8xfm]

    with i2run(args) as job:
        gemmi.read_structure(str(job / "XYZOUT.cif"))
        gemmi.read_mtz_file(str(job / "FPHIOUT.mtz"))
        gemmi.read_mtz_file(str(job / "DIFFPHIOUT.mtz"))
