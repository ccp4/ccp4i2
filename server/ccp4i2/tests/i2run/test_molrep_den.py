import gemmi
from .utils import demoData, i2run


def test_gamma_patterson():
    """Test molrep density-mode (Patterson) search with gamma data."""
    args = ["molrep_den"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--PERFORM", "pat"]
    with i2run(args) as job:
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0
