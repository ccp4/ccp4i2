from .utils import demoData, i2run


def test_glyco_4iid():
    """Test privateer glycan validation with 4iid."""
    args = ["privateer"]
    args += ["--XYZIN", demoData("glyco", "4iid.pdb")]
    args += ["--F_SIGF", f"fullPath={demoData('glyco', '4iid.mtz')}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--SINGLETHREADED", "True"]

    with i2run(args) as job:
        assert (job / "program.xml").exists(), "No program.xml output"
