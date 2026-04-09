import gemmi
from .utils import demoData, i2run


def test_gamma_model():
    """Test adding fractional coordinates to gamma model."""
    args = ["add_fractional_coords"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        xyzout = job / "XYZOUT.cif"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0, "Output structure has no chains"
