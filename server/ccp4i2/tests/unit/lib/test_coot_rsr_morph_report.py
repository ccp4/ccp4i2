"""RSR morph says how far it moved the model.

Its report said "Full reporting is not yet available in this task" and its
output was "XYZOUT.pdb": nothing told a user whether the model had moved a
few tenths of an Angstrom (the local correction morphing is for) or several.
"""
import xml.etree.ElementTree as ET

import pytest

gemmi = pytest.importorskip("gemmi")
morph = pytest.importorskip("ccp4i2.wrappers.coot_rsr_morph.script.coot_rsr_morph")


def model(path, shift_last=0.0):
    structure = gemmi.Structure()
    m, chain = gemmi.Model("1"), gemmi.Chain("A")
    for i in range(1, 4):
        residue = gemmi.Residue()
        residue.name, residue.seqid = "ALA", gemmi.SeqId(i, " ")
        atom = gemmi.Atom()
        atom.name, atom.element = "CA", gemmi.Element("C")
        atom.pos = gemmi.Position(3.8 * i + (shift_last if i == 3 else 0.0), 0, 0)
        residue.add_atom(atom)
        chain.add_residue(residue)
    m.add_chain(chain)
    structure.add_model(m)
    structure.write_pdb(str(path))
    return path


def test_shifts_and_report(tmp_path):
    from lxml import etree
    from ccp4i2.wrappers.coot_rsr_morph.script.coot_rsr_morph_report import coot_rsr_morph_report
    moved = morph.shifts(model(tmp_path / "in.pdb"), model(tmp_path / "out.pdb", 0.9))
    assert (moved.get("atoms"), moved.get("max"), moved.get("mean")) == ("3", "0.90", "0.30")
    assert moved.get("rms") == "0.52"
    assert moved.find("Residue").get("name") == "A/ALA 3"
    root = etree.Element("coot_rsr_morph"); root.append(moved)
    report = coot_rsr_morph_report(xmlnode=ET.fromstring(etree.tostring(root)), jobInfo={},
                                   jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "moved 3 atoms: RMS shift 0.52" in text and "largest 0.90" in text
    assert "A/ALA 3" in text
