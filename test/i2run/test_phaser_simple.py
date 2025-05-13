import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_gamma():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["PHASER.1", "XYZOUT_REFMAC", "XYZOUT_SHEETBEND"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        for name in ["DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        assert min(rworks) < 0.27
        llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
        assert max(llgs) > 1000


def test_no_solution():
    args = ["phaser_simple"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("rnase", "rnase_model.pdb")]
    args += ["--RESOLUTION_HIGH", "5.0"]
    with i2run(args) as job:
        assert not (job / "PHASER.1.pdb").exists()
