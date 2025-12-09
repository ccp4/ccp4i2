import xml.etree.ElementTree as ET
from pytest import approx
from .utils import i2run


def test_8xfm(cif8xfm, seq8xfm):
    args = ["modelASUCheck"]
    args += ["--XYZIN", cif8xfm]
    args += ["--ASUIN", f"seqFile={seq8xfm}"]
    with i2run(args) as job:
        tree = ET.parse(job / "program.xml")
        alignments = tree.findall(".//SequenceAlignment/Alignment")
        assert len(alignments) == 1
        assert int(alignments[0].find("match_count").text) == 276
        assert float(alignments[0].find("identity_1").text) == approx(87.898, abs=1e-3)
        assert float(alignments[0].find("identity_2").text) == approx(100.0)
        assert alignments[0].find("CIGAR").text == "1I159M19I54M5I63M13I"
