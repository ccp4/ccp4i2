import xml.etree.ElementTree as ET
from gemmi import CoorFormat, read_mtz_file, read_structure
from .utils import demoData, i2run


def test_gamma():
    args = ["prosmart_refmac"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--NCYCLES", "4"]
    args += ["--USEANOMALOUS", "True"]
    with i2run(args) as job:
        read_structure(str(job / "XYZOUT.pdb"), format=CoorFormat.Pdb)
        read_structure(str(job / "CIFFILE.pdb"), format=CoorFormat.Mmcif)
        for name in ["ABCD", "ANOMFPHI", "DIFANOMFPHI", "DIFFPHI", "FPHI"]:
            read_mtz_file(str(job / f"{name}OUT.mtz"))
        xml = ET.parse(job / "program.xml")
        cycles = xml.findall(".//RefmacInProgress/Cycle")
        rworks = [float(c.find("r_factor").text) for c in cycles]
        assert len(cycles) == 5
        assert rworks[0] > rworks[-1]
        assert xml.find(".//Molprobity_score") is not None
        assert xml.find(".//B_factors/all[@chain='All']") is not None
        assert xml.find(".//Ramachandran/Totals") is not None
