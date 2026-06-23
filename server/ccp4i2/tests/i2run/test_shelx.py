import re
from gemmi import read_mtz_file, read_pdb
from .utils import demoData, i2run


def test_substrdet():
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
    args += ["--END_PIPELINE", "substrdet"]
    with i2run(args) as job:
        read_pdb(str(job / "PDBCUR.pdb"))

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
        # Thresholds give headroom for CCP4-version numerical drift: the
        # ccp4-20251105 suite landed R-work ~0.21, the 20260520 suite ~0.242
        # for an identical (correctly phased, FOM ~0.78) result. A genuine
        # phasing/build failure produces R-factors of 0.4+, well above these.
        _check_output(job, min_fom=0.7, max_rwork=0.27, max_rfree=0.30)


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
        # See note in test_gamma_sad: thresholds tolerate CCP4-version drift
        # (20260520 suite lands R-work ~0.244, FOM ~0.765 for a good result).
        _check_output(job, min_fom=0.7, max_rwork=0.27, max_rfree=0.30)


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
