import gemmi
from .utils import demoData, i2run


def test_gamma_chainsaw():
    """Test chainsaw model editing with gamma PIR alignment."""
    args = ["chainsaw"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ALIGNIN", demoData("gamma", "gamma.pir")]
    with i2run(args) as job:
        xyzout = job / "XYZOUT.pdb"
        assert xyzout.exists(), f"No XYZOUT: {list(job.iterdir())}"
        st = gemmi.read_structure(str(xyzout))
        assert len(st[0]) > 0, "Output structure has no chains"
