import xml.etree.ElementTree as ET
from .utils import demoData, i2run


def test_gamma_phase_comparison():
    """Test phase comparison between two sets of phases."""
    args = ["cphasematch"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    args += ["--ABCD1", demoData("gamma", "initial_phases.mtz")]
    args += ["--ABCD2", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        abcdout = job / "ABCDOUT.mtz"
        assert abcdout.exists(), f"No ABCDOUT: {list(job.iterdir())}"
