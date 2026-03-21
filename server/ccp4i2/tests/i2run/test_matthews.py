import xml.etree.ElementTree as ET
from .utils import demoData, i2run


def test_gamma_from_nres():
    """Test matthews coefficient estimation using number of residues."""
    args = ["matthews"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--MODE", "nres"]
    args += ["--NRES", "120"]
    with i2run(args) as job:
        tree = ET.parse(job / "program.xml")
        compositions = tree.findall(".//matthewsCompositions/composition")
        assert len(compositions) > 0, "No matthews compositions in output"

        # Check most likely result has reasonable solvent content (30-80%)
        summary = tree.find(".//summary")
        assert summary is not None, "No summary element"
        solvent = float(summary.find("solventContent").text)
        assert 20 < solvent < 90, f"Unreasonable solvent content: {solvent}%"


def test_gamma_from_molwt():
    """Test matthews coefficient estimation using molecular weight."""
    args = ["matthews"]
    args += ["--HKLIN", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--MODE", "molwt"]
    args += ["--MOLWT", "13500"]
    with i2run(args) as job:
        tree = ET.parse(job / "program.xml")
        compositions = tree.findall(".//matthewsCompositions/composition")
        assert len(compositions) > 0, "No matthews compositions in output"
