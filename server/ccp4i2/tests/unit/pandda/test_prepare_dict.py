"""``prepare_dict_for_pandda`` (design note §4.7, assertions in §15.2).

A ``value_order`` dictionary comes out carrying a ``type`` column whose tokens
the old PanDDA reader's map accepts; a ``type``-only dictionary passes through
untouched; both round-trip through gemmi ``ChemComp`` with identical bond
orders. CCP4-free: gemmi only.
"""
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.wrappers.pandda_campaign.script.pandda_dict import (
    UnknownBondOrder, bond_spellings, needs_normalising, prepare_dict_for_pandda,
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
])
def test_type_token_mapping(value_order, flag, expected):
    assert type_token(value_order, flag) == expected


# ---------------------------------------------------------------------------
# Full-word value_order tokens (what acedrg writes), and unknown tokens.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value_order, flag, expected", [
    ("DOUB", "N", "double"), ("DOUBLE", "N", "double"), ("double", "", "double"),
    ("SING", "N", "single"), ("SINGLE", "N", "single"),
    ("TRIP", "N", "triple"), ("TRIPLE", "N", "triple"),
    ("AROM", "Y", "aromatic"), ("AROMATIC", "N", "aromatic"),
    ("DELO", "N", "deloc"), ("DELOC", "N", "deloc"),
    ("SINGLE", "Y", "aromatic"),          # the aromatic flag wins
])
def test_every_spelling_of_an_order_maps_to_the_readers_token(value_order, flag, expected):
    assert type_token(value_order, flag) == expected


def test_an_unknown_order_refuses_rather_than_guessing():
    """Reversed from "connected beats collapsed": a guessed order produced a
    wrong molecule that nothing downstream could notice."""
    with pytest.raises(UnknownBondOrder, match="QUAD"):
        type_token("QUAD")


ACETAMIDE_FULL_WORDS = """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
LIG LIG 'acetamide' non-polymer 9 4
data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
LIG C1 C CH3 0 0.0 0.0 0.0
LIG C2 C C 0 1.5 0.0 0.0
LIG O1 O O 0 2.1 1.1 0.0
LIG N1 N NH2 0 2.1 -1.2 0.0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LIG C1 C2 SINGLE N 1.500 0.016
LIG C2 O1 DOUBLE N 1.226 0.017
LIG C2 N1 SINGLE N 1.353 0.012
"""


def test_a_double_bond_spelled_in_full_stays_double_when_staged(tmp_path):
    """The bug: DOUBLE fell to the single fallback, and PanDDA, which reads
    only the added type column, built a tetrahedral amide."""
    import gemmi
    src = tmp_path / "amide.cif"
    src.write_text(ACETAMIDE_FULL_WORDS)
    staged = prepare_dict_for_pandda(src, tmp_path / "staging")
    assert staged != src, "a value_order-only dictionary must be rewritten with a type column"
    block = gemmi.cif.read(str(staged))["comp_LIG"]
    types = {(r[0], r[1]): r[2] for r in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "type"])}
    assert types[("C2", "O1")] == "double"
    assert types[("C2", "N1")] == "single" and types[("C1", "C2")] == "single"


def test_a_dictionary_with_an_unknown_order_is_refused_at_staging(tmp_path):
    src = tmp_path / "odd.cif"
    src.write_text(ACETAMIDE_FULL_WORDS.replace("C2 O1 DOUBLE", "C2 O1 QUAD"))
    with pytest.raises(UnknownBondOrder, match="QUAD"):
        prepare_dict_for_pandda(src, tmp_path / "staging")
