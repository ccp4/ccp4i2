from .utils import demoData, i2run


def test_aimless():
    args = ["aimless_pipe"]
    args += ["--UNMERGEDFILES", demoData("gamma", "gamma_native.mtz")]
    with i2run(args) as directory:
        assert directory.exists()
