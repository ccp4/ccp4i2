import pytest
import gemmi
from .utils import demoData, i2run


def test_gamma_merge_two_mtz():
    """Test merging two mini-MTZ files from gamma data."""
    mtz1 = demoData("gamma", "merged_intensities_native.mtz")
    mtz2 = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["mergeMtz"]
    args += ["--MINIMTZINLIST", f"fileName={mtz1}"]
    args += ["--MINIMTZINLIST", f"fileName={mtz2}"]
    with i2run(args) as job:
        hklout = job / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job.iterdir())}"
        gemmi.read_mtz_file(str(hklout))
