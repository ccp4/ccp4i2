import re
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
        read_pdb(str(job / "n_part.pdb"))
        read_pdb(str(job / "n_PDBCUR.pdb"))
        _check_output(job, min_fom=0.7, max_rwork=0.21, max_rfree=0.25)


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
        _check_output(job, min_fom=0.7, max_rwork=0.22, max_rfree=0.28)


def _check_output(job, min_fom, max_rwork, max_rfree):
    for name in ["FPHOUT_2FOFC", "FPHOUT_DIFF", "FPHOUT_HL", "FREEROUT"]:
        read_mtz_file(str(job / f"{name}.mtz"))
    log = (job / "log.txt").read_text()
    foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
    rworks = [float(x) for x in re.findall(r"R factor .* (0\.\d+)", log)]
    rfrees = [float(x) for x in re.findall(r"R-free .* (0\.\d+)", log)]
    assert max(foms) > min_fom
    assert min(rworks) < max_rwork
    assert min(rfrees) < max_rfree
