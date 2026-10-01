"""Single-atom MR says what it found.

Its outputs were "Positioned coordinates for solution 1", "SingleMR.1.mtz"
(no annotation) and "H-L Co-efficients1", and its report opened with
Phaser's own text: the R-factor and the atoms placed (on 3njw, two sulfurs
completed to 132 atoms, R 22.3%) were nowhere a user would look.
"""
import xml.etree.ElementTree as ET

import pytest

gemmi = pytest.importorskip("gemmi")
smr = pytest.importorskip("ccp4i2.wrappers.phaser_singleMR.phaser_singleMR")

LOG = """
   Final Log-Likelihood = 6051.11
   Final R-factor =  25.9
   Final Log-Likelihood = 7194.38
   Final R-factor =  22.6
   Final Log-Likelihood = 7347.56
   Final R-factor =  22.3
"""


def write_atoms(path, elements):
    structure = gemmi.Structure()
    model, chain = gemmi.Model("1"), gemmi.Chain("A")
    for i, element in enumerate(elements, 1):
        residue = gemmi.Residue()
        residue.name, residue.seqid = "UNK", gemmi.SeqId(i, " ")
        atom = gemmi.Atom()
        atom.name, atom.element = element, gemmi.Element(element)
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    structure.write_pdb(str(path))
    return path


def test_summary(tmp_path):
    pdb = write_atoms(tmp_path / "SingleMR.1.pdb", ["S", "S"] + ["N"] * 130)
    summary = smr.summarise(LOG, pdb)
    assert [(c.get("llg"), c.get("r")) for c in summary.findall("Cycle")] == [
        ("6051.11", "25.9"), ("7194.38", "22.6"), ("7347.56", "22.3")]
    atoms = summary.find("Atoms")
    assert atoms.get("total") == "132"
    assert {e.get("name"): e.get("count") for e in atoms} == {"N": "130", "S": "2"}


def test_report_leads_with_the_result(tmp_path):
    from ccp4i2.wrappers.phaser_singleMR.phaser_singleMR_report import phaser_singleMR_report
    pdb = write_atoms(tmp_path / "SingleMR.1.pdb", ["S", "S"] + ["N"] * 130)
    root = ET.fromstring("<Job><SubJob_000><RunDate>today</RunDate>"
                         "<KeyText_003>There were 6 solutions</KeyText_003></SubJob_000></Job>")
    root.append(smr.summarise(LOG, pdb))
    report = phaser_singleMR_report(xmlnode=root, jobInfo={}, jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "Solution 1 has 132 atoms (130 N, 2 S); after 3 rounds" in text
    assert "R 22.3%, LLG 7347.56" in text


def test_acorn_report_states_the_correlation():
    """ACORN's report gave its correlation by cycle only as a plot."""
    from ccp4i2.wrappers.acorn.script.acorn_report import acorn_report
    cycles = "".join(f"<Cycle><NCycle>{n}</NCycle><CorrelationCoef>{cc}</CorrelationCoef></Cycle>"
                     for n, cc in enumerate([0.0, 0.55216, 0.65769, 0.68181, 0.69033, 0.68844, 0.68642]))
    root = ET.fromstring(f"<acorn><RunInfo>{cycles}</RunInfo></acorn>")
    report = acorn_report(xmlnode=root, jobInfo={}, jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "0.552 after cycle 1, best 0.690 (cycle 4), final 0.686 after 6 cycles" in text
