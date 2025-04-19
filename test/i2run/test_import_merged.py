from .utils import demoData, i2run


def test_import_merged():
    args = ["import_merged"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    with i2run(args) as job:
        assert job.exists()
