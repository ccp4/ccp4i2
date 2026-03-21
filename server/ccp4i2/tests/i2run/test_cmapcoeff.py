import gemmi
from .utils import demoData, i2run


def test_gamma_anomalous():
    """Test anomalous map coefficient calculation with gamma Xe data."""
    mtz = demoData("gamma", "merged_intensities_Xe.mtz")
    args = ["cmapcoeff"]
    args += ["--F_SIGF1", mtz]
    args += ["--ABCD1", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        fphiout = job / "FPHIOUT.mtz"
        assert fphiout.exists(), f"No FPHIOUT: {list(job.iterdir())}"
        gemmi.read_mtz_file(str(fphiout))
