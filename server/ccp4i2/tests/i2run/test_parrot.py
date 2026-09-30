import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


def test_parrot():
    args = ["parrot"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--ABCD", demoData("gamma", "initial_phases.mtz")]
    args += ["--ASUIN", demoData("gamma", "gamma.asu.xml")]
    with i2run(args) as job:
        for name in ["ABCDOUT", "FPHIOUT"]:
            gemmi.read_mtz_file(str(job / f"{name}.mtz"))
        xml = ET.parse(job / "program.xml")
        foms = [float(e.text) for e in xml.findall(".//MeanFOM")]
        assert max(foms) > 0.794


def test_parrot_honours_asu_copies():
    """The ASU says 6 copies of AHIR: seqin must carry all 6, so parrot's
    solvent content is that of 6 copies (0.632), not a Matthews guess (it
    chose 9 copies, 0.448, when seqin held one)."""
    args = ["parrot"]
    args += ["--F_SIGF", demoData("ahir", "ahir_observed.mtz")]
    args += ["--ABCD", demoData("ahir", "ahir_phases_mr.mtz")]
    args += ["--ASUIN", demoData("ahir", "ahir.asu.xml")]
    args += ["--CYCLES", "1"]
    with i2run(args) as job:
        with open(job / "seqin.fasta", encoding="utf-8") as fh:
            assert fh.read().count(">") == 6
        xml = ET.parse(job / "program.xml")
        solvent = float(xml.find("./SolventContent/SolventContent").text)
        assert abs(solvent - 0.632) < 0.005
