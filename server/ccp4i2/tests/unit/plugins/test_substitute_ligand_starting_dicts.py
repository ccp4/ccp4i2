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


# --- the same check, made at the Run dialog instead of four steps in ---------

class _RecordedCheck:
    """Stands in for checkMonomeCoverage, recording how it was called."""

    def __init__(self, plugin, complaint="No dictionary at all for: ZZZ"):
        self.plugin = plugin
        self.complaint = complaint
        self.model = None
        self.dict_paths = None

    def __call__(self, xyzin_path, dict_paths=None):
        self.model = str(xyzin_path)
        self.dict_paths = list(dict_paths or [])
        if self.complaint:
            self.plugin.appendErrorReport(350, f"Monomer coverage check failed:\n{self.complaint}")
            return self.plugin.FAILED
        return self.plugin.SUCCEEDED


def _model(tmp_path, name="start.pdb", residue=None):
    path = tmp_path / name
    if residue:
        path.write_text(
            f"HETATM    1  C1  {residue} A 401       0.000   0.000   0.000  1.00 20.00           C  \n"
            "END\n",
            encoding="utf-8",
        )
    else:
        path.write_text("END\n", encoding="utf-8")
    return path


def _prepared(tmp_path, ligand_as="SMILES", xyzin=None, complaint="No dictionary at all for: ZZZ"):
    plugin = _plugin("SubstituteLigand", tmp_path)
    plugin.container.inputData.XYZIN.setFullPath(str(xyzin or _model(tmp_path)))
    plugin.container.controlParameters.LIGANDAS.set(ligand_as)
    check = _RecordedCheck(plugin, complaint=complaint)
    plugin.checkMonomeCoverage = check
    return plugin, check


def _codes(error):
    return [e["code"] for e in error.entries()]


# --- coverage: every ligand already in the model needs a dictionary ---------

def test_an_undescribed_ligand_is_refused_before_the_job_runs(tmp_path):
    plugin, _check = _prepared(tmp_path)
    assert 213 in _codes(plugin.runTimeValidity()), \
        "the pipeline would have run four steps to learn this"


def test_the_refusal_says_what_to_do_about_it(tmp_path):
    plugin, _check = _prepared(tmp_path)
    details = " ".join(e["details"] for e in plugin.runTimeValidity().entries())
    assert "ZZZ" in details, "the residue is named"
    assert "atom selection" in details, "remedy one: take it out of the model"
    assert "supply its dictionary" in details, "remedy two: describe it"


def test_the_refusal_points_at_the_coordinates(tmp_path):
    plugin, _check = _prepared(tmp_path)
    names = [e["name"] for e in plugin.runTimeValidity().entries() if e["code"] == 213]
    assert names == ["SubstituteLigand.container.inputData.XYZIN"]


def test_the_code_the_pipeline_will_mint_excuses_nothing(tmp_path):
    """A DRG in the starting model is a different molecule from the DRG we make.

    An earlier draft waved this through on the grounds that a DRG dictionary
    was coming. It describes the ligand being fitted, and covers nothing that
    is already there.
    """
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="No dictionary at all for: DRG")
    assert 213 in _codes(plugin.runTimeValidity())


def test_only_dictionaries_for_the_starting_model_count(tmp_path, old_ligand_dict):
    """DICTIN describes the ligand being fitted, so it covers nothing here."""
    plugin, check = _prepared(tmp_path, ligand_as="DICT")
    plugin.container.inputData.STARTING_DICT_LIST.append(old_ligand_dict)
    substituted = tmp_path / "substituted.cif"
    substituted.write_text("data_comp_LIG\n", encoding="utf-8")
    plugin.container.inputData.DICTIN.setFullPath(str(substituted))

    plugin.runTimeValidity()

    assert check.dict_paths == [old_ligand_dict]


def test_a_covered_model_passes(tmp_path):
    plugin, _check = _prepared(tmp_path, complaint="")
    assert 213 not in _codes(plugin.runTimeValidity())


def test_the_check_reads_the_selected_atoms_not_the_file_as_supplied(tmp_path):
    """Excluding a ligand with the atom selection is a real answer."""
    plugin, check = _prepared(tmp_path, complaint="")
    supplied = str(plugin.container.inputData.XYZIN.fullPath)
    plugin.runTimeValidity()
    assert check.model != supplied, \
        "checked the file as supplied, so the atom selection could never clear it"


def test_an_unusable_selection_falls_back_to_the_file(tmp_path):
    plugin, check = _prepared(tmp_path, complaint="")

    def refuse(*args, **kwargs):
        raise RuntimeError("selection cannot be applied")

    plugin.container.inputData.XYZIN.getSelectedAtomsFile = refuse
    plugin.runTimeValidity()
    assert check.model == str(plugin.container.inputData.XYZIN.fullPath), \
        "a check that cannot apply the selection must still be made"


# --- identity: two ligands cannot share a code ------------------------------

def test_a_clash_with_the_fitted_ligands_code_is_refused(tmp_path):
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="")
    assert 214 in _codes(plugin.runTimeValidity())


def test_the_two_problems_are_reported_separately(tmp_path):
    """An undescribed DRG in the model is both: it needs a dictionary, and it
    clashes. Fixing either one alone does not make the job runnable."""
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="No dictionary at all for: DRG")
    codes = _codes(plugin.runTimeValidity())
    assert 213 in codes and 214 in codes


def test_the_clash_message_offers_both_ways_out(tmp_path):
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="")
    details = " ".join(e["details"] for e in plugin.runTimeValidity().entries())
    assert "LIG" in details, "say what to call it instead"
    assert "atom selection" in details, "or take the old copy out"
    assert "repositioned against new data" in details, "the same-ligand case, named"


def test_the_clash_points_at_the_code_not_the_coordinates(tmp_path):
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="")
    names = [e["name"] for e in plugin.runTimeValidity().entries() if e["code"] == 214]
    assert names == ["SubstituteLigand.container.controlParameters.LIGAND_CODE"]


def test_choosing_another_code_resolves_it(tmp_path):
    plugin, _check = _prepared(tmp_path, xyzin=_model(tmp_path, residue="DRG"),
                               complaint="")
    plugin.container.controlParameters.LIGAND_CODE.set("LIG")
    assert 214 not in _codes(plugin.runTimeValidity())


def test_a_supplied_dictionary_clashes_by_the_code_it_declares(tmp_path):
    """In DICT mode the pipeline does not choose the code; the file does."""
    plugin, _check = _prepared(tmp_path, ligand_as="DICT",
                               xyzin=_model(tmp_path, residue="LIG"), complaint="")
    dictin = tmp_path / "supplied.cif"
    dictin.write_text(
        "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.name\n"
        "LIG  ligand\ndata_comp_LIG\n", encoding="utf-8")
    plugin.container.inputData.DICTIN.setFullPath(str(dictin))

    assert 214 in _codes(plugin.runTimeValidity())


def test_no_clash_when_no_ligand_is_being_fitted(tmp_path):
    plugin, _check = _prepared(tmp_path, ligand_as="NONE",
                               xyzin=_model(tmp_path, residue="DRG"), complaint="")
    assert 214 not in _codes(plugin.runTimeValidity())


def test_the_chosen_code_is_what_acedrg_and_coot_are_told(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert plugin._generatedLigandCode() == "DRG"
    plugin.container.controlParameters.LIGAND_CODE.set("lig")
    assert plugin._generatedLigandCode() == "LIG", "codes are upper case"


def test_an_unreadable_model_is_not_read_as_containing_nothing(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert plugin._residueCodesIn(tmp_path / "no_such_file.pdb") is None
