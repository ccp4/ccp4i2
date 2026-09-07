from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


def _mr_args(task, copies=1):
    args = [task]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ENSEMBLES", "label=model", "use=True", f"number={copies}",
             "pdbItemList/identity_to_target=0.9",
             f"pdbItemList/structure={demoData('gamma', 'gamma_model.pdb')}"]
    args += ["--COMP_BY", "DEFAULT"]
    return args


@pytest.mark.order("first")
def test_rnp_refines_the_solution_of_an_mr_auto_run():
    # The RNP run nests inside the MR run: the harness's test database
    # belongs to the outer context
    with i2run(_mr_args("phaser_mr_auto_phil")) as mr_job:
        solutions = mr_job / "SOLOUT.phaser_sol.pkl"
        assert solutions.is_file()
        mr_xml = ET.parse(mr_job / "program.xml")
        mr_llg = max(float(e.text) for e in mr_xml.findall("Solutions/Solution/LLG"))
        assert mr_xml.findtext("Solutions/Solution/spaceGroup")
        with i2run(_mr_args("phaser_mr_rnp_phil", copies=0) + ["--SOLIN", str(solutions)]) as job:
            gemmi.read_pdb(str(job / "PHASER.1.pdb"))
            xml = ET.parse(job / "program.xml")
            llgs = [float(e.text) for e in xml.findall("Solutions/Solution/LLG")]
            # refined in place, at RNP's own protocol: the same solution, a comparable LLG
            assert llgs and max(llgs) > 1000 and max(llgs) > 0.5 * mr_llg
            assert "MR_RNP" in (job / "working.phil").read_text()
            assert "MOLECULAR REPLACEMENT REFINEMENT AND PHASING" in [
                m.get("name") for m in xml.findall("Modules/Module")]
