import gemmi
from .utils import demoData, i2run


def test_hl_to_phifom():
    """Test chltofom HL→PHI/FOM conversion with gamma initial phases."""
    args = ["chltofom"]
    args += ["--HKLIN", demoData("gamma", "initial_phases.mtz")]

    with i2run(args) as job:
        # Check output MTZ exists
        hklout = job / "HKLOUT.mtz"
        assert hklout.exists(), f"No HKLOUT: {list(job.iterdir())}"

        # Verify output has PHI and FOM columns
        mtz = gemmi.read_mtz_file(str(hklout))
        labels = [c.label for c in mtz.columns]
        assert "PHI" in labels, f"No PHI column in output: {labels}"
        assert "FOM" in labels, f"No FOM column in output: {labels}"
