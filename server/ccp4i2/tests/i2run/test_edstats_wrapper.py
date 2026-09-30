from .utils import i2run


def test_8xfm(cif8xfm, mtz8xfm):
    """Test edstats with 8xfm structure and map coefficients."""
    args = ["edstats"]
    args += ["--XYZIN", cif8xfm]
    args += ["--FPHIIN1", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    args += ["--FPHIIN2", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[DELFWT,PHDELWT]"]
    args += ["--RES_LOW", "50.0"]
    args += ["--RES_HIGH", "1.3"]

    with i2run(args) as job:
        assert (job / "program.xml").exists(), "No program.xml output"


def test_8xfm_resolution_from_the_map_coefficients(cif8xfm, mtz8xfm):
    """With no resolution given, edstats takes the range of the map
    coefficients (it once refused to run until someone typed it in)."""
    args = ["edstats"]
    args += ["--XYZIN", cif8xfm]
    args += ["--FPHIIN1", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    args += ["--FPHIIN2", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[DELFWT,PHDELWT]"]

    with i2run(args) as job:
        assert (job / "program.xml").exists(), "No program.xml output"
        params = (job / "params.xml").read_text()
        assert "<RES_HIGH>" in params and "<RES_LOW>" in params, params
