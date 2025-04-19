from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm):
    args = ["csymmatch"]
    args += ["--XYZIN_QUERY", cif8xfm]
    args += ["--XYZIN_TARGET", cif8xfm]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
