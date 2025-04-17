import gemmi
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_molrep():
    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["XYZOUT_MOLREP", "XYZOUT_SHEETBEND", "XYZOUT"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
