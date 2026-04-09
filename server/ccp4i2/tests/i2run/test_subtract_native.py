from .utils import i2run
import gemmi


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["SubtractNative"]
    args += ["--XYZIN", cif8xfm]
    args += ["--MAPIN", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FWT,PHWT]"]
    with i2run(args) as job:
        gemmi.read_ccp4_map(str(job / "MAPOUT.map"))
