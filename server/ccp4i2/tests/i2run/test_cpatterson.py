from .utils import demoData, i2run


def test_gamma_native():
    """Test Patterson map calculation from gamma native data."""
    args = ["cpatterson"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    with i2run(args) as job:
        assert (job / "MAPOUT.map").exists(), f"No MAPOUT: {list(job.iterdir())}"
