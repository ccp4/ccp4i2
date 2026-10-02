import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


def _args(sheetbend="False"):
    args = ["phaser_simple_phil"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--RUNREFMAC", "False", "--RUNSHEETBEND", sheetbend]
    return args


@pytest.mark.order("first")
def test_gamma_basic():
    with i2run(_args()) as job:
        gemmi.read_pdb(str(job / "PHASER.1.pdb"))
        record = ET.parse(job / "program.xml").find("PhaserMrResults")
        assert record.findtext("Verdict") == "Single Solution"
        assert max(float(e.text) for e in record.findall("Solutions/Solution/LLG")) > 1000
        assert record.findtext("Solutions/Solution/Components/Component/Name") == "SearchModel"
        # What Phaser expected of the model, as numbers a judge can read
        assert float(record.findtext("ExpectedLLG/Ensemble/eLLG")) > 0
        assert 0 < float(record.findtext("ExpectedLLG/FractionScatteringOfEnsembles")) <= 1


@pytest.mark.order("first")
def test_gamma_sheetbend():
    with i2run(_args(sheetbend="True")) as job:
        gemmi.read_pdb(str(job / "XYZOUT_SHEETBEND.pdb"))
        assert ET.parse(job / "program.xml").find("SheetbendResult") is not None
