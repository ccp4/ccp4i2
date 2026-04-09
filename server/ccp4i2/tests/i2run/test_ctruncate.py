import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="KeywordExtractor bug: get_merged_metadata is None for ctruncate container")
def test_gamma_intensities():
    """Test ctruncate with gamma anomalous intensities."""
    args = ["ctruncate"]
    args += ["--OBSIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]

    with i2run(args) as job:
        # Check for output MTZ
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from ctruncate: {list(job.iterdir())}"
