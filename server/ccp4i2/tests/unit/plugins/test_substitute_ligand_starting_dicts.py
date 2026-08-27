"""
A ligand already present in the starting model needs somewhere to declare
itself.

SubstituteLigand refines with servalcat before it fits the new ligand, so what
servalcat sees is the *starting* model --- including any ligand already in it.
That happens two ways: a model built against different data still carrying the
old ligand, and a genuine second soak on a crystal that already contained one.
Servalcat's pre-flight monomer check refuses a residue it has no dictionary
for, and SubstituteLigand passed it no dictionaries at all: `DICTIN` describes
the ligand being substituted and goes to coot, after refinement. So the job
stopped at refinement with no way past it, however well the user understood
what was wrong.

These tests cover the route, not the check --- that a dictionary offered to
SubstituteLigand arrives in the list servalcat actually reads.
"""
from pathlib import Path

import pytest

from ccp4i2.core.tasks import get_plugin_class


@pytest.fixture(name="old_ligand_dict")
def old_ligand_dict_fixture(tmp_path):
    path = tmp_path / "old_ligand.cif"
    path.write_text("data_comp_list\ndata_comp_LIG\n", encoding="utf-8")
    return str(path)


def _plugin(task, tmp_path):
    plugin_class = get_plugin_class(task)
    assert plugin_class is not None, f"{task} is not in the registry"
    directory = tmp_path / task
    directory.mkdir(exist_ok=True)
    return plugin_class(workDirectory=str(directory), name=task)


def test_the_input_exists_and_starts_empty(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert len(plugin.container.inputData.STARTING_DICT_LIST) == 0


def test_unset_entries_are_not_offered(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert plugin._startingLigandDictionaries() == []


def test_a_supplied_dictionary_is_offered_by_full_path(tmp_path, old_ligand_dict):
    plugin = _plugin("SubstituteLigand", tmp_path)
    plugin.container.inputData.STARTING_DICT_LIST.append(old_ligand_dict)
    assert plugin._startingLigandDictionaries() == [old_ligand_dict]


def test_the_dictionary_reaches_the_list_servalcat_reads(tmp_path, old_ligand_dict):
    """Both ends of the route: SubstituteLigand's input, servalcat's DICT_LIST.

    servalcat_pipe.runTimeValidity() builds its dictionary paths from
    inputData.DICT_LIST and hands them to checkMonomeCoverage; this is the
    list that has to be non-empty for an already-present ligand to be
    describable at all.
    """
    substitute = _plugin("SubstituteLigand", tmp_path)
    substitute.container.inputData.STARTING_DICT_LIST.append(old_ligand_dict)
    servalcat = _plugin("servalcat_pipe", tmp_path)

    for dictionary in substitute._startingLigandDictionaries():
        servalcat.container.inputData.DICT_LIST.append(dictionary)

    seen_by_the_check = [str(d.fullPath) for d in servalcat.container.inputData.DICT_LIST
                         if d.isSet()]
    assert seen_by_the_check == [old_ligand_dict]


def test_supplying_a_dictionary_does_not_take_it_out_of_the_pipeline(tmp_path, old_ligand_dict):
    """Appending the object itself would re-parent it; paths are passed instead."""
    substitute = _plugin("SubstituteLigand", tmp_path)
    substitute.container.inputData.STARTING_DICT_LIST.append(old_ligand_dict)
    servalcat = _plugin("servalcat_pipe", tmp_path)

    for dictionary in substitute._startingLigandDictionaries():
        servalcat.container.inputData.DICT_LIST.append(dictionary)

    assert substitute._startingLigandDictionaries() == [old_ligand_dict], (
        "the pipeline lost its own input when handing it on")
