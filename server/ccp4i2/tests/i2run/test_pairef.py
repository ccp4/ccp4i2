from .utils import i2run


def test_pairef(cif8xfm, mtz8xfm):
    args = ["pairef"]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--XYZIN", cif8xfm]
    args += ["--INIRES", "3.0"]
    with i2run(args) as job:
        path = job / "pairef_project" / "PAIREF_cutoff.txt"
        with open(path, encoding="utf-8") as f:
            assert float(f.read()) == 2.95
