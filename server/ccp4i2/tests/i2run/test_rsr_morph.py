from .utils import hasLongLigandName, i2run


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["coot_rsr_morph"]
    args += ["--XYZIN", cif8xfm]
    args += ["--FPHIIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        assert hasLongLigandName(job / "XYZOUT.cif")
