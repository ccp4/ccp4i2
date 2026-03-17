import gemmi
from .utils import i2run


def test_self_superposition(cif8xfm):
    """Test gesamt pairwise superposition (identity)."""
    args = ["gesamt"]
    args += ["--XYZIN_QUERY", cif8xfm]
    args += ["--XYZIN_TARGET", cif8xfm]

    with i2run(args) as job:
        gemmi.read_structure(str(job / "XYZOUT.cif"))
