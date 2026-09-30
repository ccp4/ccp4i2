"""MakeLink writes AceDRG's link instructions from a declared description.

The words themselves are pinned against what AceDRG's parser (covLink.py,
getInstructionsForLinkFreeFormat2) accepts; each shape asserted here was run
through `acedrg -L` from ccp4-20260702 when this was written. The i2run test
test_make_link.py::test_several_edits_to_one_monomer runs one end to end.
"""
import pytest

from ccp4i2.pipelines.MakeLink.script.link_instruction import (
    MonomerEdits, edit_problems, edit_words, extra_instruction_problems,
)


# --- edit_words: one canonical form ------------------------------------------

def test_every_deletion_gets_its_own_delete():
    # "DELETE ATOM OE2 2 ATOM OXT 2" stops AceDRG with "Unknown keyword ATOM".
    words = edit_words(MonomerEdits(deletes=["OE2", "OXT"]), 2)
    assert words == ["DELETE", "ATOM", "OE2", "2", "DELETE", "ATOM", "OXT", "2"]


def test_bond_and_charge_changes_share_one_change_section():
    edits = MonomerEdits(bond_orders=[("CD", "OE1", "double")], charges=[("NZ", 1)])
    assert " ".join(edit_words(edits, 1)) == "CHANGE BOND CD OE1 DOUBLE 1 CHARGE 1 NZ 1"


def test_charge_puts_the_monomer_number_first():
    # Every other item ends with the monomer number; CHARGE starts with it.
    assert edit_words(MonomerEdits(charges=[("N1", -1)]), 2) == [
        "CHANGE", "CHARGE", "2", "N1", "-1"]


def test_deletions_come_before_changes_whatever_order_they_were_declared():
    edits = MonomerEdits(bond_orders=[("C4A", "O4A", "SINGLE")], deletes=["O4A"])
    # (Contradictory, but edit_words is only asked after edit_problems.)
    assert edit_words(edits, 2)[:4] == ["DELETE", "ATOM", "O4A", "2"]


def test_a_repeated_edit_is_written_once():
    edits = MonomerEdits(
        deletes=["OE2", "OE2"],
        bond_orders=[("CD", "OE1", "DOUBLE"), ("OE1", "CD", "DOUBLE")],
        charges=[("NZ", 1), ("NZ", 1)])
    assert " ".join(edit_words(edits, 1)) == (
        "DELETE ATOM OE2 1 CHANGE BOND CD OE1 DOUBLE 1 CHARGE 1 NZ 1")


def test_no_edits_no_words():
    assert edit_words(MonomerEdits(), 1) == []


def test_monomer_must_be_one_or_two():
    with pytest.raises(ValueError):
        edit_words(MonomerEdits(deletes=["O"]), 3)


# --- edit_problems: contradictions are found before AceDRG runs ------------

def test_a_consistent_description_has_no_problems():
    edits = MonomerEdits(deletes=["OE2"], bond_orders=[("CD", "OE1", "DOUBLE")])
    assert edit_problems(edits, link_atom="CD") == []


def test_the_linking_atom_cannot_be_deleted():
    problems = edit_problems(MonomerEdits(deletes=["NZ"]), link_atom="NZ")
    assert problems and "linking atom" in problems[0]


def test_a_bond_to_a_deleted_atom_cannot_change_order():
    edits = MonomerEdits(deletes=["O4A"], bond_orders=[("C4A", "O4A", "SINGLE")])
    assert any("O4A is deleted" in p for p in edit_problems(edits))


def test_a_deleted_atom_cannot_change_charge():
    edits = MonomerEdits(deletes=["N1"], charges=[("N1", 1)])
    assert any("N1 cannot change charge" in p for p in edit_problems(edits))


def test_one_bond_cannot_have_two_orders_whichever_way_round():
    edits = MonomerEdits(bond_orders=[("CD", "OE1", "DOUBLE"), ("OE1", "CD", "SINGLE")])
    assert any("two orders" in p for p in edit_problems(edits))


def test_one_atom_cannot_have_two_charges():
    edits = MonomerEdits(charges=[("NZ", 1), ("NZ", 0)])
    assert any("two charges" in p for p in edit_problems(edits))


@pytest.mark.parametrize("edits", [
    MonomerEdits(deletes=[""]),
    MonomerEdits(deletes=["O 4A"]),
    MonomerEdits(bond_orders=[("CD", "", "DOUBLE")]),
    MonomerEdits(bond_orders=[("CD", "CD", "DOUBLE")]),
    MonomerEdits(bond_orders=[("CD", "OE1", "QUADRUPLE")]),
    MonomerEdits(charges=[("", 1)]),
])
def test_a_malformed_edit_is_a_problem(edits):
    # A blank or spaced name would shift every later argument out of place.
    assert edit_problems(edits)


# --- extra_instruction_problems: the free-text box ----------------------------

@pytest.mark.parametrize("text", [
    "",
    "# only comments\n#   DELETE ATOM X 1\n",
    "DELETE ATOM C4A 1",
    "DELETE ATOM C4A 1\nDELETE ATOM O4A 1",
    "DELETE BOND C1 C2 2",
    "CHANGE BOND C O DOUBLE 2 CHARGE 2 SD 1",
    "change bond c o double 2",
    "ADD ATOM X1 C 0 1 BOND X1 C4 SINGLE 1",
    "LINK: DELETE ATOM C4A 1",
])
def test_well_formed_extra_instructions_pass(text):
    assert extra_instruction_problems(text) == []


@pytest.mark.parametrize("text, fragment", [
    # The shape that hangs AceDRG: an unknown word inside CHANGE.
    ("CHANGE OE1 DOUBLE 2", 'expected BOND or CHARGE, not "OE1"'),
    ("CHANGE BOND CD OE1 DOUBLE 2 OXT", 'not "OXT"'),
    ("ADD ATOM X1 C 0 1 FOO", 'not "FOO"'),
    # One item per DELETE.
    ("DELETE ATOM OE2 2 ATOM OXT 2", '"ATOM" is not a keyword here'),
    ("BOND CD OE1 DOUBLE 2", "must begin with DELETE, CHANGE or ADD"),
    ("DELETE ATOM C4A 3", "monomer number"),
    ("DELETE ATOM C4A", "needs 2 values"),
    ("CHANGE CHARGE SD 2 1", "monomer number"),
    ("CHANGE CHARGE 2 SD plus", "whole-number charge"),
    ("CHANGE BOND C O QUADRUPLE 2", "bond order"),
    ("DELETE", "must be followed"),
])
def test_malformed_extra_instructions_are_refused(text, fragment):
    problems = extra_instruction_problems(text)
    assert len(problems) == 1 and fragment in problems[0], problems
