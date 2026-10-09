import gemmi

from .utils import demoData, i2run


def test_gamma_intensities():
    """Anomalous intensities in, anomalous and mean amplitudes out as data objects."""
    args = ["ctruncate"]
    args += ["--OBSIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]

    with i2run(args) as job:
        fpair = gemmi.read_mtz_file(str(job / "OBSOUT.mtz"))
        assert {"Fplus", "SIGFplus", "Fminus", "SIGFminus"} <= set(fpair.column_labels())
        fmean = gemmi.read_mtz_file(str(job / "OBSOUT_asFMEAN.mtz"))
        assert {"F", "SIGF"} <= set(fmean.column_labels())
