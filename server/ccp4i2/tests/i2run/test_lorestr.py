import gemmi
from .utils import i2run


def test_lorestr_1h1s(pdb1h1s, mtz1h1s):
    """Test lorestr low-resolution refinement with 1h1s data (no auto homologue fetch)."""
    args = ["lorestr_i2"]
    args += ["--XYZIN", f"fullPath={pdb1h1s}"]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FREE]"]
    args += ["--AUTO", "none"]
    args += ["--CPU", "1"]

    with i2run(args) as job:
        # Check refined model output
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        gemmi.read_pdb(str(xyzout))

        # Check map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"
