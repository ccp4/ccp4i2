from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm):
    args = ["coordinate_selector"]
    args += ["--XYZIN", f"fullPath={cif8xfm}", "selection/text=(A1LU6)"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
