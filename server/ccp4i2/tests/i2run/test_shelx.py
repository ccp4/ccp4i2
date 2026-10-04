import re
import gemmi
from gemmi import read_mtz_file, read_pdb
from .programs import requires_program
from .utils import demoData, i2run


@requires_program("shelxc", "shelxd")
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

@requires_program("shelxc", "shelxd", "shelxe")
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
        _check_output(job, max_rwork=0.27, max_rfree=0.30, min_trace_cc=40)
        # The SHELX route: SHELXE density modification and tracing ran. Only
        # the React interface used to select it, so every i2run/API "shelx"
        # job ran crank2's own route instead, and this test passed on that.
        log = (job / "log.txt").read_text(errors="replace")
        assert "Running shelxe with hand 1" in log, "SHELXE did not run"
        # The model to refine carries the xenon (XYZOUT leaves it out), and
        # carries each atom once
        import xml.etree.ElementTree as ET
        done = ET.parse(job / "program.xml").find("CompleteModel")
        assert done is not None and done.get("error") is None, done.attrib if done is not None else None
        complete = gemmi.read_structure(str(job / "XYZOUT_COMPLETE.pdb"))[0]
        built = gemmi.read_structure(str(job / "n_part.pdb"))[0]
        xenons = [a for ch in complete for r in ch for a in r if a.element.name == "Xe"]
        assert len(xenons) == int(done.get("added")) >= 1
        assert complete.count_atom_sites() == built.count_atom_sites() + len(xenons)


@requires_program("shelxc", "shelxd", "shelxe")
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
        _check_output(job, max_rwork=0.27, max_rfree=0.30, min_trace_cc=40)


def _check_output(job, max_rwork, max_rfree, min_fom=None, min_trace_cc=None):
    """min_fom for crank2's route (its "FOM is" lines); min_trace_cc for the
    SHELX route, whose phasing measure is SHELXE's trace CC (above ~25% is
    usually solved)."""
    for name in ["FPHOUT_2FOFC", "FPHOUT_DIFF", "FPHOUT_HL", "FREEROUT"]:
        read_mtz_file(str(job / f"{name}.mtz"))
    log = (job / "log.txt").read_text()
    foms = [float(x) for x in re.findall(r"FOM is (0\.\d+)", log)]
    rworks = [float(x) for x in re.findall(r"R factor .* (0\.\d+)", log)]
    rfrees = [float(x) for x in re.findall(r"R-free .* (0\.\d+)", log)]
    if min_fom is not None:
        assert max(foms) > min_fom
    if min_trace_cc is not None:
        ccs = [float(x) for x in re.findall(r"best correlation coef\. (\d+\.?\d*)", log)]
        assert ccs, "SHELXE reported no trace CC"
        assert max(ccs) > min_trace_cc
    assert min(rworks) < max_rwork
    assert min(rfrees) < max_rfree
