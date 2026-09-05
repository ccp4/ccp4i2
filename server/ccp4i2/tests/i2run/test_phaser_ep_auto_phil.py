from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


@pytest.mark.order("first")
def test_gamma_xe():
    args = ["phaser_ep_auto_phil"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_HA", demoData("gamma", "heavy_atoms.pdb")]
    args += ["--ASUFILE", demoData("gamma", "gamma.asu.xml")]
    args += ["--COMP_BY", "ASU"]
    args += ["--WAVELENGTH", "1.542"]
    args += ["--LLGC_CYCLES", "20"]
    args += ["--ELEMENTS", "Xe"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "PHASER.1.pdb"))
        for mtz in ("ABCDOUT_1", "MAPOUT_1"):
            gemmi.read_mtz_file(str(job / f"{mtz}.mtz"))
        xml = ET.parse(job / "program.xml")
        assert len(xml.findall("Hands/Hand")) >= 1
        assert float(xml.findtext("Overall/fom")) > 0.3
        assert len(xml.findall("Sites/Site")) >= 4
        assert xml.findtext("Sites/Site/Element") == "XE"
        completion = xml.find("Completion")
        assert len(completion.findall("Cycle")) >= 1
        assert float(completion.get("llg")) > 900
        assert len(xml.findall("Summaries/Summary")) >= 2
        assert "composition" in (job / "working.phil").read_text()
