from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["sheetbend"]
    args += ["--XYZIN", cif8xfm]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
