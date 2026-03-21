import gemmi
from .utils import demoData, i2run


def test_gamma_merge_two_mtz():
    """Test merging two mini-MTZ files from gamma data."""
    mtz = demoData("gamma", "initial_phases.mtz")
    args = ["mergeMtz"]
    args += ["--MINIMTZINLIST", f"filename={mtz}"]
    with i2run(args) as job:
        hklout = job / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job.iterdir())}"
        gemmi.read_mtz_file(str(hklout))
