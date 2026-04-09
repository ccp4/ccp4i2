import xml.etree.ElementTree as ET
import gemmi
import pytest
from .utils import i2run


@pytest.mark.order("first")
def test_phaser_rnp_pipeline_chain_selections(pdb1h1s, mtz1h1s):
    """Test phaser_rnp_pipeline with per-chain selections from PDB 1h1s."""
    args = ["phaser_rnp_pipeline"]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--XYZIN_PARENT", pdb1h1s]
    args += ["--SELECTIONS", "text=A/"]
    args += ["--SELECTIONS", "text=B/"]
    args += ["--SELECTIONS", "text=C/"]
    args += ["--SELECTIONS", "text=D/"]
    args += ["--RUNREFMAC", "False"]
    args += ["--RUNSHEETBEND", "False"]
    with i2run(args) as job:
        # Phaser should produce at least one set of output coordinates
        pdb_files = list(job.glob("PHASER.*.pdb")) + list(job.glob("PHASER.*.cif"))
        assert len(pdb_files) > 0, "Phaser produced no output coordinate files"
        for pdb_file in pdb_files:
            gemmi.read_structure(str(pdb_file))

        # Check for map outputs
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) > 0, "No MTZ output files produced"

        # Parse program.xml for LLG scores
        xml_path = job / "program.xml"
        if xml_path.exists():
            xml = ET.parse(xml_path)
            llgs = [float(e.text) for e in xml.findall(".//Solution/LLG")]
            if llgs:
                assert max(llgs) > 0, f"LLG scores should be positive, got {llgs}"
