import gemmi
from .utils import i2run


def test_zanuda_8xfm(cif8xfm, mtz8xfm):
    """Test zanuda space group validation with 8xfm mmCIF input."""
    args = ["zanuda"]
    args += ["--XYZIN", f"fullPath={cif8xfm}"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]

    with i2run(args) as job:
        # Zanuda writes output as zanuda.pdb / zanuda.mtz (not XYZOUT)
        xyzout = job / "zanuda.pdb"
        assert xyzout.exists(), f"No zanuda.pdb: {list(job.iterdir())}"
        gemmi.read_pdb(str(xyzout))

        # Check split map coefficients
        for name in ["FPHIOUT", "DIFFPHIOUT"]:
            mtz_file = job / f"{name}.mtz"
            assert mtz_file.exists(), f"No {name}: {list(job.iterdir())}"
