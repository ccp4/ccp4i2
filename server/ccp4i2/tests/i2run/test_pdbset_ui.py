import gemmi
from .utils import demoData, i2run


def test_gamma_model():
    """Test pdbset with a simple CRYST1 keyword."""
    args = ["pdbset_ui"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--EXTRA_PDBSET_KEYWORDS", "CHAIN A"]
    with i2run(args) as job:
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0
