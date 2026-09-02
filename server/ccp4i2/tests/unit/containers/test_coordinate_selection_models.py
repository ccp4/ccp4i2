"""Model-qualified selections, against structures gemmi actually built.

A CID names a model first: '/1/A/23(HOH)' is model 1, chain A, residue 23.
The evaluator matched that against ``Model.name``, a string gemmi carried
until 0.7 and does not carry now, so every model-qualified selection raised
TypeError on the fallback path --- reaching ``get_subchain(0)``, which wants a
subchain name, not an index. Nothing caught it, because the only structures
these tests had ever run against were mocks that still offered ``name``.

So this file uses real gemmi structures throughout. It also covers the second
half of the same mistake: the evaluator worked the model number out from the
model's *position* in the structure. A file numbered 3 and 7 has no model 1,
and calling its first model 1 selects the wrong atoms rather than none ---
the failure that says nothing.

The identifier is ``Model.num``, an int, and gemmi preserves what it read.
"""
import textwrap

import gemmi
import pytest

from ccp4i2.core.coordinate_selection.evaluator import evaluate_selection
from ccp4i2.core.coordinate_selection.parser import parse_selection


def _structure(numbers, chain="A"):
    """A structure with one CA atom per model, in models numbered *numbers*.

    Each model's atom sits at x = its model number, so which model was
    selected is readable off the coordinates.
    """
    text = ""
    for n in numbers:
        text += textwrap.dedent(f"""\
            MODEL     {n:>4}
            ATOM      1  CA  ALA {chain}   1     {float(n):7.3f}   0.000   0.000  1.00 20.00           C
            ENDMDL
            """)
    return gemmi.read_structure_string(text + "END\n", format=gemmi.CoorFormat.Pdb)


def _select(structure, selection):
    return evaluate_selection(parse_selection(selection), structure)


def _models_selected(structure, selection):
    """The model numbers a selection picked out."""
    return sorted({model.num for model, _chain, _res, _atom in _select(structure, selection)})


# --- the crash ---------------------------------------------------------------

def test_a_model_qualified_selection_does_not_raise():
    """It used to reach get_subchain(0) and raise TypeError."""
    structure = _structure([1])
    assert _select(structure, "/1/A") == _select(structure, "/1/A"), "evaluated twice without raising"


def test_a_model_qualified_selection_selects_the_model():
    structure = _structure([1])
    assert len(_select(structure, "/1/A")) == 1


def test_the_whole_of_one_model_can_be_asked_for():
    structure = _structure([1, 2])
    assert _models_selected(structure, "/1/*") == [1]


# --- numbers are what gemmi read, not positions ------------------------------

def test_model_numbers_are_preserved_as_read():
    """The premise the rest of this file rests on."""
    assert [m.num for m in _structure([3, 7])] == [3, 7]


def test_a_model_is_found_by_its_number():
    structure = _structure([3, 7])
    assert _models_selected(structure, "/7/A") == [7]


def test_position_is_not_the_number():
    """Numbering by position would call these 1 and 2 and select the wrong one."""
    structure = _structure([3, 7])
    assert _models_selected(structure, "/1/A") == [], "there is no model 1 in this file"
    assert _models_selected(structure, "/2/A") == []


def test_the_atoms_come_from_the_model_that_was_asked_for():
    """Selecting the wrong model returns atoms, so a count alone proves nothing."""
    structure = _structure([3, 7])
    (_model, _chain, _residue, atom), = _select(structure, "/7/A")
    assert atom.pos.x == pytest.approx(7.0)


# --- the pattern forms _matches_value offers, over model numbers -------------

def test_a_wildcard_model_takes_every_model():
    structure = _structure([3, 7])
    assert _models_selected(structure, "/*/A") == [3, 7]


def test_an_unqualified_selection_takes_every_model():
    structure = _structure([3, 7])
    assert _models_selected(structure, "A") == [3, 7]


def test_a_list_of_model_numbers():
    structure = _structure([1, 2, 3])
    assert _models_selected(structure, "/1,3/A") == [1, 3]


def test_a_range_of_model_numbers():
    structure = _structure([1, 2, 3, 4])
    assert _models_selected(structure, "/2-3/A") == [2, 3]


# --- the form the rest of the codebase emits ---------------------------------

def test_the_cid_coordinate_selector_writes_can_be_read_back():
    """coordinate_selector.py builds ids as '/{model.num}/{chain}/{seq}({name})'.

    It was writing CIDs the evaluator could not evaluate.
    """
    structure = _structure([1])
    model = structure[0]
    chain = model[0]
    residue = chain[0]
    cid = f"/{model.num}/{chain.name}/{residue.seqid.num}({residue.name})"

    assert len(_select(structure, cid)) == 1, f"{cid} selected nothing"


def test_a_single_model_file_still_answers_to_model_one():
    """The common case, which position-based numbering got right by accident."""
    assert _models_selected(_structure([1]), "/1/A") == [1]
