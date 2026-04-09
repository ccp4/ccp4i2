from .utils import demoData, i2run


def test_reindex_to_coords():
    """Test pointless reindexing gamma data to match coordinate reference."""
    args = ["pointless_reindexToMatch"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_REF", demoData("gamma", "gamma_model.pdb")]
    args += ["--REFERENCE", "XYZIN_REF"]

    with i2run(args) as job:
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from pointless: {list(job.iterdir())}"
