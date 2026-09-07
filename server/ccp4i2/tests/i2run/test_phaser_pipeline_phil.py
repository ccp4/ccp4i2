from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import demoData, i2run


def _beta_blip_args():
    args = ["phaser_pipeline_phil"]
    args += ["--F_SIGF", f"fullPath={demoData('beta_blip', 'beta_blip_P3221.mtz')}",
             "columnLabels=/*/*/[Fobs,Sigma]"]
    for model in ("beta", "blip"):
        args += ["--ENSEMBLES", f"label={model}", "use=True", "number=1",
                 "pdbItemList/identity_to_target=0.9",
                 f"pdbItemList/structure={demoData('beta_blip', model + '.pdb')}"]
    return args


@pytest.mark.order("first")
def test_beta_blip_asu_no_refinement():
    args = _beta_blip_args() + ["--COMP_BY", "ASU",
                                "--ASUFILE", f"seqFile={demoData('beta_blip', 'beta.seq')}",
                                f"seqFile={demoData('beta_blip', 'blip.seq')}",
                                "--RUNSHEETBEND", "False", "--RUNREFMAC", "False"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "PHASER.1.pdb"))
        for mtz in ("DIFMAPOUT_1", "MAPOUT_1", "PHASEOUT_1"):
            gemmi.read_mtz_file(str(job / f"{mtz}.mtz"))
        xml = ET.parse(job / "program.xml")
        record = xml.find("PhaserMrResults")
        assert record is not None and record.findtext("Verdict") == "Single Solution"
        llgs = [float(e.text) for e in record.findall("Solutions/Solution/LLG")]
        assert max(llgs) > 1000
        names = [c.text for c in record.findall("Solutions/Solution/Components/Component/Name")]
        assert sorted(names) == ["beta", "blip"]
        assert record.find("Strategy").get("unparsed") == "0"
        # the pipeline hosted the tree and handed it down: the same working phil
        phil = (job / "CCP4_JOBS" / "job_1" / "working.phil") if (job / "CCP4_JOBS").exists() else None
        assert "sheetbend" not in "".join(f.name for f in job.iterdir()).lower()


@pytest.mark.order("first")
def test_beta_blip_with_sheetbend_and_refmac():
    args = _beta_blip_args() + ["--COMP_BY", "DEFAULT", "--RUNSHEETBEND", "True", "--RUNREFMAC", "True"]
    with i2run(args) as job:
        gemmi.read_pdb(str(job / "XYZOUT_SHEETBEND.pdb"))
        gemmi.read_pdb(str(job / "XYZOUT_REFMAC.pdb"))
        xml = ET.parse(job / "program.xml")
        assert xml.find("SheetbendResult") is not None and xml.find("REFMAC") is not None
        rworks = [float(e.text) for e in xml.iter("r_factor")]
        assert rworks and min(rworks) < 0.3
