"""A search model's UNK chains are found, per chain (#680)."""

import gemmi

from ccp4i2.lib.utils.formats.unknown_residues import describe, unknown_residue_chains


def write_model(path, chains):
    structure = gemmi.Structure()
    structure.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    model = gemmi.Model("1")
    for name, residues in chains.items():
        chain = gemmi.Chain(name)
        for number, resname in enumerate(residues, start=1):
            residue = gemmi.Residue()
            residue.name = resname
            residue.seqid = gemmi.SeqId(number, " ")
            atom = gemmi.Atom()
            atom.name = "CA"
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(number * 3.8, 0, 0)
            residue.add_atom(atom)
            chain.add_residue(residue)
        model.add_chain(chain)
    structure.add_model(model)
    structure.setup_entities()
    structure.make_mmcif_document().write_file(str(path))
    return path


def test_unk_chains_are_reported_and_sequenced_chains_are_not(tmp_path):
    path = write_model(tmp_path / "built.cif", {
        "A": ["ALA", "GLY", "SER", "LEU", "VAL"],
        "B": ["UNK"] * 8 + ["ALA", "GLY"],
    })
    assert unknown_residue_chains(path) == [("B", 8, 10)]
    assert describe(unknown_residue_chains(path)) == "chain B has 8 UNK of 10 residues"


def test_a_fully_sequenced_model_has_none(tmp_path):
    path = write_model(tmp_path / "model.cif", {"A": ["ALA", "GLY", "SER", "LEU"]})
    assert unknown_residue_chains(path) == []


def test_molrep_warns_but_does_not_block_on_an_unk_search_model(tmp_path):
    import pytest
    pytest.importorskip("django")
    import os
    import django
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.test_settings")
    django.setup()
    from ccp4i2.core import CCP4ErrorHandling
    from ccp4i2.wrappers.molrep_mr.script.molrep_mr import molrep_mr

    path = write_model(tmp_path / "built.cif", {"A": ["UNK"] * 6})
    plugin = molrep_mr(parent=None, workDirectory=str(tmp_path))
    plugin.container.inputData.XYZIN.setFullPath(str(path))
    report = plugin.runTimeValidity()
    unk = [e for e in report._errors if e.get("code") == 180]
    assert len(unk) == 1
    assert "chain A has 6 UNK of 6 residues" in str(unk[0].get("details"))
    assert unk[0].get("severity") == CCP4ErrorHandling.SEVERITY_WARNING
    assert unk[0].get("name") == "molrep_mr.container.inputData.XYZIN"
