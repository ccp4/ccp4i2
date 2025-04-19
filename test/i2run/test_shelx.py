from gemmi import read_mtz_file, read_pdb
from .utils import demoData, i2run


def test_gamma_sad():
    args = ["shelx"]
    args += [
        "--F_SIGFanom",
        f"fullPath={demoData('gamma', 'merged_intensities_Xe.mtz')}",
        "columnLabels=/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
    ]
    args += ["--SEQIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--ATOM_TYPE", "Xe"]
    args += ["--NUMBER_SUBSTRUCTURE", "2"]
    args += ["--WAVELENGTH", "1.54179"]
    args += ["--FPRIME", "-0.79"]
    args += ["--FDPRIME", "7.36"]
    with i2run(args) as job:
        for name in ["n_part", "n_PDBCUR"]:
            read_pdb(str(job / f"{name}.pdb"))
        for name in ["FPHOUT_2FOFC", "FPHOUT_DIFF", "FPHOUT_HL", "FREEROUT"]:
            read_mtz_file(str(job / f"{name}.mtz"))


def test_gamma_siras():
    args = ["shelx"]
    args += [
        "--F_SIGFanom",
        f"fullPath={demoData('gamma', 'merged_intensities_Xe.mtz')}",
        "columnLabels=/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
    ]
    args += [
        "--F_SIGFnative",
        f"fullPath={demoData('gamma', 'merged_intensities_native.mtz')}",
        "columnLabels=/*/*/[Iplus,SIGIplus,Iminus,SIGIminus]"
    ]
    args += ["--SEQIN", demoData("gamma", "gamma.asu.xml")]
    args += ["--NATIVE", "True"]
    args += ["--EXPTYPE", "SIRAS"]
    args += ["--ATOM_TYPE", "Xe"]
    args += ["--NUMBER_SUBSTRUCTURE", "2"]
    args += ["--WAVELENGTH", "1.54179"]
    args += ["--FPRIME", "-0.79"]
    args += ["--FDPRIME", "7.36"]
    with i2run(args) as job:
        read_pdb(str(job / "n_REFMAC5.pdb"))
        for name in ["FPHOUT_2FOFC", "FPHOUT_DIFF", "FPHOUT_HL", "FREEROUT"]:
            read_mtz_file(str(job / f"{name}.mtz"))
