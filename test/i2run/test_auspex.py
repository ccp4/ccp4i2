from .utils import demoData, i2run


def test_auspex():
    args = ["AUSPEX"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    with i2run(args) as directory:
        for plot in ["F", "FSigF", "I", "ISigI", "SigF", "SigI"]:
            assert (directory / f"{plot}_plot.png").exists()
