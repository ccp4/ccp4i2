"""editbfac says what it kept and how it split the model.

Its outputs were annotated with their file names (converted_model_chainA1.pdb
and so on) and its report said "Edit B-factors finished." above the program's
log, so which file held which domain could be learnt only by opening them.
"""
import xml.etree.ElementTree as ET

import pytest

gemmi = pytest.importorskip("gemmi")
editbfac = pytest.importorskip("ccp4i2.wrappers.editbfac.script.editbfac")


def model(path, chain, numbers):
    structure = gemmi.Structure()
    m = gemmi.Model("1")
    c = gemmi.Chain(chain)
    for n in numbers:
        r = gemmi.Residue()
        r.name, r.seqid = "ALA", gemmi.SeqId(n, " ")
        a = gemmi.Atom()
        a.name, a.element = "CA", gemmi.Element("C")
        r.add_atom(a)
        c.add_residue(r)
    m.add_chain(c)
    structure.add_model(m)
    structure.setup_entities()
    structure.write_pdb(str(path))
    return path


def test_residue_ranges():
    assert editbfac.residue_ranges([28, 26, 27, 40, 41, 50]) == "26-28, 40-41, 50"
    assert editbfac.residue_ranges([]) == ""


@pytest.fixture
def summary(tmp_path):
    full = model(tmp_path / "in.pdb", "A", range(1, 492))
    kept = model(tmp_path / "converted_model.pdb", "A", [*range(26, 112), *range(435, 491)])
    a1 = model(tmp_path / "converted_model_chainA1.pdb", "A1", range(26, 112))
    a2 = model(tmp_path / "converted_model_chainA2.pdb", "A2", range(435, 491))
    return editbfac.summarise(full, kept, [a1, a2])


def test_summary_counts_and_ranges(summary):
    assert summary.findtext("InputResidues") == "491"
    assert summary.find("Model").get("residues") == "142"
    assert summary.find("Model").get("ranges") == "26-111, 435-490"
    domains = [(d.get("chain"), d.get("ranges"), d.get("residues"))
               for d in summary.findall("Domain")]
    assert domains == [("A1", "26-111", "86"), ("A2", "435-490", "56")]


def test_report_leads_with_what_was_kept(summary, tmp_path, monkeypatch):
    from lxml import etree
    from ccp4i2.wrappers.editbfac.script.editbfac_report import editbfac_report
    monkeypatch.setattr(editbfac_report, "getJobFolder", lambda self: str(tmp_path))
    (tmp_path / "log.txt").write_text("Maximum B-value to be included: 59.22 A**2 < 60\n")
    report = editbfac_report(xmlnode=ET.fromstring(etree.tostring(summary)),
                             jobInfo={}, jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "142 of 491 residues kept: 26-111, 435-490." in text
    assert "Split into 2 regions" in text
    assert "435-490" in text and "A2" in text
    assert "59.22" in text   # the log, escaped, in its fold
