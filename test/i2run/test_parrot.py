import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def test_parrot():
    args = ["parrot"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ABCD", demoData("gamma", "initial_phases.mtz")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        foms = [float(e.text) for e in xml.findall(".//MeanFOM")]
        assert max(foms) > 0.8
