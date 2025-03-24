import gemmi
from .utils import demoData, i2run


def test_molrep():
    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as directory:
        for name in ["XYZOUT_MOLREP", "XYZOUT_SHEETBEND", "XYZOUT"]:
            gemmi.read_pdb(str(directory / f"{name}.pdb"))
