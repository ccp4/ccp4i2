from .utils import demoData, i2run


def test_scaleit():
    args = ["scaleit"]
    args += ["--MERGEDFILES", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--MERGEDFILES", demoData("gamma", "merged_intensities_Xe.mtz")]
    with i2run(args) as job:
        pass
