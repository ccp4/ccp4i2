from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


def _args(copies=1, comp_by="DEFAULT"):
    args = ["phaser_mr_auto_phil"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ENSEMBLES", "label=model", "use=True", f"number={copies}",
             "pdbItemList/identity_to_target=0.9",
             f"pdbItemList/structure={demoData('gamma', 'gamma_model.pdb')}"]
    args += ["--COMP_BY", comp_by]
    return args


@pytest.mark.order("first")
def test_gamma_single_solution():
    with i2run(_args()) as job:
        gemmi.read_pdb(str(job / "PHASER.1.pdb"))
        for mtz in ("MAPOUT_1", "DIFMAPOUT_1", "PHASEOUT_1"):
            gemmi.read_mtz_file(str(job / f"{mtz}.mtz"))
        xml = ET.parse(job / "program.xml")
        assert xml.findtext("Verdict") == "Single Solution"
        llgs = [float(e.text) for e in xml.findall("Solutions/Solution/LLG")]
        assert max(llgs) > 4000
        assert xml.findtext("Solutions/Solution/Components/Component/Name") == "model"
        assert len(xml.findall("Modules/Module")) >= 10
        assert len(xml.findall("Summaries/Summary")) >= 10
        strategy = xml.find("Strategy")
        assert strategy.get("unparsed") == "0"
        assert strategy.find("Attempt").get("outcome").startswith("placed")
        assert xml.find("Solutions/Solution/UnknownTokens") is None


@pytest.mark.order("first")
def test_gamma_asu_composition():
    args = _args(comp_by="ASU") + ["--ASUFILE", demoData("gamma", "gamma.asu.xml")]
    with i2run(args) as job:
        assert (job / "composition_1.fasta").is_file()
        assert ET.parse(job / "program.xml").findtext("Verdict") == "Single Solution"
        phil = (job / "working.phil").read_text()
        assert "composition" in phil and "sequence_file" in phil
