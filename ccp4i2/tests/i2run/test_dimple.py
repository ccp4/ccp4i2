import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def test_dimple():
    args = ["i2Dimple"]  # Plugin is named i2Dimple in registry
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_native.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["FPHIOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
