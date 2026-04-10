import xml.etree.ElementTree as ET
from gemmi import CoorFormat, read_mtz_file, read_structure
from .utils import demoData, hasLongLigandName, i2run


def _check_output(job, anom, expected_cycles, expected_rwork):
    read_structure(str(job / "XYZOUT.pdb"), format=CoorFormat.Pdb)
    read_structure(str(job / "CIFFILE.pdb"), format=CoorFormat.Mmcif)
    for name in ["ABCD", "ANOMFPHI", "DIFANOMFPHI", "DIFFPHI", "FPHI"]:
        if anom or "ANOM" not in name:
            read_mtz_file(str(job / f"{name}OUT.mtz"))
    xml = ET.parse(job / "program.xml")
    cycles = xml.findall(".//RefmacInProgress/Cycle")
    rworks = [float(c.find("r_factor").text) for c in cycles]
    assert len(rworks) == expected_cycles
    assert rworks[-1] < rworks[0]
    assert rworks[-1] < expected_rwork
    assert xml.find(".//Molprobity_score") is not None
    assert xml.find(".//B_factors/all[@chain='All']") is not None
    assert xml.find(".//Ramachandran/Totals") is not None


def test_8xfm(cif8xfm, mtz8xfm):
    args = ["prosmart_refmac"]
    args += ["--XYZIN", cif8xfm]
    args += ["--F_SIGF", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz8xfm}", "columnLabels=/*/*/[FREE]"]
    args += ["--NCYCLES", "2"]
    args += ["--ADD_WATERS", "True"]
    with i2run(args) as job:
        _check_output(job, anom=False, expected_cycles=6, expected_rwork=0.19)
        assert hasLongLigandName(job / "CIFFILE.pdb")


def test_gamma():
    args = ["prosmart_refmac"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--NCYCLES", "4"]
    args += ["--USEANOMALOUS", "True"]
    with i2run(args) as job:
        _check_output(job, anom=True, expected_cycles=5, expected_rwork=0.27)
