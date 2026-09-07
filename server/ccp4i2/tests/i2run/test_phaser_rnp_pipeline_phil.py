"""phaser_rnp_pipeline_phil: a parent model cut into rigid bodies, refined
with Phaser over its PHIL."""
import xml.etree.ElementTree as ET

import gemmi
import pytest

from .utils import i2run


@pytest.mark.order("first")
def test_chain_selections(pdb1h1s, mtz1h1s):
    args = ["phaser_rnp_pipeline_phil"]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--XYZIN_PARENT", pdb1h1s]
    for chain in "ABCD":
        args += ["--SELECTIONS", f"text={chain}/"]
    args += ["--RUNREFMAC", "False", "--RUNSHEETBEND", "False"]
    with i2run(args) as job:
        for index in range(1, 5):
            assert (job / f"Fragment_{index}.pdb").is_file()
        gemmi.read_pdb(str(job / "PHASER.1.pdb"))
        gemmi.read_mtz_file(str(job / "MAPOUT_1.mtz"))
        gemmi.read_mtz_file(str(job / "F_SIGF_OUT.mtz"))
        xml = ET.parse(job / "program.xml")
        assert xml.find("POINTLESS") is not None
        record = xml.find("PhaserMrResults")
        assert record is not None
        llgs = [float(e.text) for e in record.findall("Solutions/Solution/LLG")]
        assert llgs and max(llgs) > 0
        names = [c.text for c in record.findall("Solutions/Solution/Components/Component/Name")]
        assert sorted(set(names)) == [f"Fragment_{i}" for i in range(1, 5)]
        # The reindexing is sub-job 1, Phaser sub-job 2
        phil = (job / "job_2" / "working.phil").read_text()
        assert "*MR_RNP" in phil and phil.count("solution_at_origin = True") == 4
