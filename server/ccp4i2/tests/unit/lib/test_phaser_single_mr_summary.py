"""Single-atom MR says what it found.

Its outputs were "Positioned coordinates for solution 1", "SingleMR.1.mtz"
(no annotation) and "H-L Co-efficients1", and its report opened with
Phaser's own text: the R-factor and the atoms placed (on 3njw, two sulfurs
completed to 132 atoms, R 22.3%) were nowhere a user would look.

Phaser completes each placement in turn, so the log's final LLG/R pairs are
one per placement (six on 3njw), not rounds of one; the kept solution is
the one with the highest LLG, wherever it falls in the log.
"""
import xml.etree.ElementTree as ET

import pytest

gemmi = pytest.importorskip("gemmi")
smr = pytest.importorskip("ccp4i2.wrappers.phaser_singleMR.phaser_singleMR")

LOG = "".join(f"   Final Log-Likelihood = {llg}\n   Final R-factor =  {r}\n" for llg, r in [
    ("4438.86", "31.0"), ("4425.49", "31.2"), ("6046.67", "26.0"),
    ("6051.17", "25.9"), ("7327.06", "22.3"), ("7347.56", "22.3")])


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
    assert len(summary.findall("Completion")) == 6
    assert (summary.find("Best").get("llg"), summary.find("Best").get("r")) == ("7347.56", "22.3")
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
    assert "Phaser completed 6 placements" in text
    assert "solution 1 (kept here), has 132 atoms (130 N, 2 S): LLG 7347.56, R 22.3%" in text


def test_best_is_the_highest_llg_not_the_last():
    log = ("Final Log-Likelihood = 900.0\nFinal R-factor = 21.0\n"
           "Final Log-Likelihood = 100.0\nFinal R-factor = 45.0\n")
    best = smr.summarise(log).find("Best")
    assert (best.get("llg"), best.get("r")) == ("900.0", "21.0")


def test_acorn_report_states_the_correlation():
    """ACORN's report gave its correlation by cycle only as a plot."""
    from ccp4i2.wrappers.acorn.script.acorn_report import acorn_report
    cycles = "".join(f"<Cycle><NCycle>{n}</NCycle><CorrelationCoef>{cc}</CorrelationCoef></Cycle>"
                     for n, cc in enumerate([0.0, 0.55216, 0.65769, 0.68181, 0.69033, 0.68844, 0.68642]))
    root = ET.fromstring(f"<acorn><RunInfo>{cycles}</RunInfo></acorn>")
    report = acorn_report(xmlnode=root, jobInfo={}, jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "0.552 after cycle 1, best 0.690 (cycle 4), final 0.686 after 6 cycles" in text
