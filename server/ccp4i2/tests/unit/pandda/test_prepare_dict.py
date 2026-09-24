"""``prepare_dict_for_pandda`` (design note §4.7, assertions in §15.2).

A ``value_order`` dictionary comes out carrying a ``type`` column whose tokens
the old PanDDA reader's map accepts; a ``type``-only dictionary passes through
untouched; both round-trip through gemmi ``ChemComp`` with identical bond
orders. CCP4-free: gemmi only.
"""
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.wrappers.pandda_campaign.script.pandda_dict import (
    FALLBACK_TYPE, bond_spellings, needs_normalising, prepare_dict_for_pandda,
    type_token)
from .conftest import pdbx_spelling, source_files

#: The exact set the CCP4-bundled reader maps (autobuild/inbuilt.py:139).
OLD_READER_TOKENS = {"single", "double", "triple", "SINGLE", "DOUBLE", "TRIPLE",
                     "aromatic", "deloc"}


def _bond_orders(path):
    doc = gemmi.cif.read(str(path))
    block = next(b for b in doc if b.find_values("_chem_comp_bond.atom_id_1"))
    cc = gemmi.make_chemcomp_from_block(block)
    return sorted((b.id1.atom, b.id2.atom, str(b.type), b.aromatic) for b in cc.rt.bonds)


def _bond_block(path):
    doc = gemmi.cif.read(str(path))
    return next(b for b in doc if b.find_values("_chem_comp_bond.atom_id_1"))


@pytest.fixture
def type_dict():
    return source_files("BAZ2BA-x425")[2]


@pytest.fixture
def value_order_dict(tmp_path, type_dict):
    path = tmp_path / "ligand_pdbx.cif"
    path.write_text(pdbx_spelling(open(type_dict).read()))
    assert bond_spellings(_bond_block(path)) == {"value_order"}, "fixture must be value_order-only"
    return path


@pytest.fixture
def named_dict(tmp_path, type_dict):
    """acedrg names the block after the ligand: comp_MZ0, not comp_LIG."""
    path = tmp_path / "MZ0.cif"
    text = open(type_dict).read().replace("comp_LIG", "comp_MZ0").replace("LIG ", "MZ0 ")
    path.write_text(text)
    assert [b.name for b in gemmi.cif.read(str(path))] == ["comp_list", "comp_MZ0"]
    return path


def test_a_ligand_named_block_gains_a_comp_lig_alias(tmp_path, named_dict):
    assert needs_normalising(named_dict)
    out = prepare_dict_for_pandda(named_dict, tmp_path / "staging")
    names = [b.name for b in gemmi.cif.read(str(out))]
    assert names == ["comp_list", "comp_MZ0", "comp_LIG"], "original first, for the content-based reader"
    doc = gemmi.cif.read(str(out))
    original, alias = doc["comp_MZ0"], doc["comp_LIG"]
    assert list(alias.find_values("_chem_comp_atom.atom_id")) == list(original.find_values("_chem_comp_atom.atom_id"))
    assert list(alias.find_values("_chem_comp_bond.type")) == list(original.find_values("_chem_comp_bond.type"))
    assert _bond_orders(out) == _bond_orders(named_dict), "the first restraint block is still the true one"
    assert prepare_dict_for_pandda(out, tmp_path / "again") == out, "idempotent"


def test_type_only_dictionary_passes_through(tmp_path, type_dict):
    out = prepare_dict_for_pandda(type_dict, tmp_path / "staging")
    assert str(out) == str(type_dict)
    assert not (tmp_path / "staging").exists(), "nothing should be written"
    assert not needs_normalising(type_dict)


def test_value_order_dictionary_gains_type_column(tmp_path, value_order_dict):
    assert needs_normalising(value_order_dict)
    out = prepare_dict_for_pandda(value_order_dict, tmp_path / "staging")
    assert out != value_order_dict
    assert out.parent == tmp_path / "staging"
    block = _bond_block(out)
    assert bond_spellings(block) == {"type", "value_order"}, "both spellings, consistently"
    types = list(block.find_values("_chem_comp_bond.type"))
    assert types and set(types) <= OLD_READER_TOKENS, types
    # aromatic derived from the PDBx flag, in the CCP4 spelling
    aromatic = list(block.find_values("_chem_comp_bond.aromatic"))
    flags = list(block.find_values("_chem_comp_bond.pdbx_aromatic_flag"))
    assert aromatic == ["y" if f == "Y" else "n" for f in flags]
    assert all(t == "aromatic" for t, f in zip(types, flags) if f == "Y")


def test_normalised_dictionary_has_the_same_bond_orders(tmp_path, type_dict, value_order_dict):
    out = prepare_dict_for_pandda(value_order_dict, tmp_path / "staging")
    assert _bond_orders(out) == _bond_orders(value_order_dict) == _bond_orders(type_dict)


def test_both_spellings_already_present_passes_through(tmp_path, value_order_dict):
    once = prepare_dict_for_pandda(value_order_dict, tmp_path / "s1")
    again = prepare_dict_for_pandda(once, tmp_path / "s2")
    assert again == once
    assert not (tmp_path / "s2").exists()


def test_other_blocks_survive_the_rewrite(tmp_path, value_order_dict):
    out = prepare_dict_for_pandda(value_order_dict, tmp_path / "staging")
    before = [(b.name, sorted(b.get_mmcif_category_names())) for b in gemmi.cif.read(str(value_order_dict))]
    after = [(b.name, sorted(b.get_mmcif_category_names())) for b in gemmi.cif.read(str(out))]
    assert after == before


def test_dictionary_without_a_bond_table_passes_through(tmp_path):
    path = tmp_path / "nobonds.cif"
    path.write_text("data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.name\nLIG 'a ligand'\n")
    assert prepare_dict_for_pandda(path, tmp_path / "staging") == path


@pytest.mark.parametrize("value_order,flag,expected", [
    ("SING", "N", "single"), ("sing", "", "single"), ("DOUB", "N", "double"),
    ("TRIP", "N", "triple"), ("AROM", "Y", "aromatic"), ("DELO", "N", "deloc"),
    ("SING", "Y", "aromatic"),   # the flag wins
    ("QUAD", "N", FALLBACK_TYPE),  # unknown: connected beats collapsed
])
def test_type_token_mapping(value_order, flag, expected):
    assert type_token(value_order, flag) == expected
