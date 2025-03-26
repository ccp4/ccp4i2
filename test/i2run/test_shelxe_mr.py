from gemmi import read_mtz_file, read_pdb
from .utils import demoData, i2run


def test_gamma():
    args = ["shelxeMR"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--NTCYCLES", "2"]
    args += ["--NMCYCLES", "2"]
    with i2run(args) as job:
        read_pdb(str(job / "shelxrun.pdb"))
        read_mtz_file(str(job / "FPHIOUT.mtz"))
