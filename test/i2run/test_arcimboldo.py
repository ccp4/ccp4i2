from .utils import demoData, i2run


def test_arcimboldo():
    args = ["arcimboldo"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ARCIMBOLDO_OPTIONS", "LITE"]
    args += ["--MOLECULAR_WEIGHT", "15000"]
    args += ["--N_FRAGMENTS", "1"]
    with i2run(args) as job:
        assert job.exists()
