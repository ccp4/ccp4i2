import gemmi
from .utils import demoData, i2run


def test_gamma():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["PHASER.1", "XYZOUT_REFMAC", "XYZOUT_SHEETBEND"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))


def test_no_solution():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("rnase", "rnase_model.pdb")]
    args += ["--RESOLUTION_HIGH", "5.0"]
    with i2run(args) as job:
        assert not (job / "PHASER.1.pdb").exists()
