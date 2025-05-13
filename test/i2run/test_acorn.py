import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def test_acorn():
    args = ["acorn"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        gemmi.read_mtz_file(str(job / "PHSOUT.mtz"))
        tree = ET.parse(job / "program.xml")
        final_cc = float(tree.findall(".//CorrelationCoef")[-1].text)
        assert final_cc > 0.1
