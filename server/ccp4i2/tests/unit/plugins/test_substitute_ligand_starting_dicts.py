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
        self.dict_paths = None
        self.ignore = None

    def __call__(self, xyzin_path, dict_paths=None, ignoreCodes=()):
        self.dict_paths = list(dict_paths or [])
        self.ignore = tuple(ignoreCodes)
        if self.complaint:
            self.plugin.appendErrorReport(350, f"Monomer coverage check failed:\n{self.complaint}")
            return self.plugin.FAILED
        return self.plugin.SUCCEEDED


def _prepared(tmp_path, ligand_as="SMILES", xyzin=None):
    plugin = _plugin("SubstituteLigand", tmp_path)
    model = xyzin or (tmp_path / "start.pdb")
    if xyzin is None:
        model.write_text("END\n", encoding="utf-8")
    plugin.container.inputData.XYZIN.setFullPath(str(model))
    plugin.container.controlParameters.LIGANDAS.set(ligand_as)
    check = _RecordedCheck(plugin)
    plugin.checkMonomeCoverage = check
    return plugin, check


def test_an_undescribed_ligand_is_refused_before_the_job_runs(tmp_path):
    plugin, _check = _prepared(tmp_path)
    error = plugin.runTimeValidity()
    codes = [e["code"] for e in error.entries()]
    assert 213 in codes, "the pipeline would have run four steps to learn this"


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


def test_the_code_the_pipeline_will_generate_is_not_held_against_it(tmp_path):
    """A model from an earlier run contains DRG, and DRG is what this run makes."""
    plugin, check = _prepared(tmp_path, ligand_as="SMILES")
    plugin.runTimeValidity()
    assert check.ignore == ("DRG",)


def test_nothing_is_excused_when_no_dictionary_will_be_generated(tmp_path):
    for ligand_as in ("NONE", "DICT"):
        plugin, check = _prepared(tmp_path, ligand_as=ligand_as)
        plugin.runTimeValidity()
        assert check.ignore == (), f"{ligand_as} generates no dictionary"


def test_every_dictionary_that_will_be_there_is_counted(tmp_path, old_ligand_dict):
    plugin, check = _prepared(tmp_path, ligand_as="DICT")
    plugin.container.inputData.STARTING_DICT_LIST.append(old_ligand_dict)
    substituted = tmp_path / "substituted.cif"
    substituted.write_text("data_comp_LIG\n", encoding="utf-8")
    plugin.container.inputData.DICTIN.setFullPath(str(substituted))

    plugin.runTimeValidity()

    assert check.dict_paths == [old_ligand_dict, str(substituted)]


def test_a_covered_model_passes(tmp_path):
    plugin, check = _prepared(tmp_path)
    check.complaint = ""
    codes = [e["code"] for e in plugin.runTimeValidity().entries()]
    assert 213 not in codes


# --- two ligands cannot share a code ----------------------------------------

@pytest.fixture(name="model_with_drg")
def model_with_drg_fixture(tmp_path):
    path = tmp_path / "already_has_drg.pdb"
    path.write_text(
        "HETATM    1  C1  DRG A 401       0.000   0.000   0.000  1.00 20.00           C  \n"
        "END\n",
        encoding="utf-8",
    )
    return path


def test_a_clash_with_the_fitted_ligands_code_is_refused(tmp_path, model_with_drg):
    """A model from an earlier run is full of DRG, and DRG is what we would mint."""
    plugin, check = _prepared(tmp_path, xyzin=model_with_drg)
    check.complaint = ""
    codes = [e["code"] for e in plugin.runTimeValidity().entries()]
    assert 214 in codes


def test_the_clash_message_offers_the_house_convention(tmp_path, model_with_drg):
    plugin, check = _prepared(tmp_path, xyzin=model_with_drg)
    check.complaint = ""
    details = " ".join(e["details"] for e in plugin.runTimeValidity().entries())
    assert "LIG" in details, "say what to call it instead"
    assert "could not be told apart" in details


def test_the_clash_points_at_the_code_not_the_coordinates(tmp_path, model_with_drg):
    plugin, check = _prepared(tmp_path, xyzin=model_with_drg)
    check.complaint = ""
    names = [e["name"] for e in plugin.runTimeValidity().entries() if e["code"] == 214]
    assert names == ["SubstituteLigand.container.controlParameters.LIGAND_CODE"]


def test_choosing_another_code_resolves_it(tmp_path, model_with_drg):
    plugin, check = _prepared(tmp_path, xyzin=model_with_drg)
    check.complaint = ""
    plugin.container.controlParameters.LIGAND_CODE.set("LIG")

    codes = [e["code"] for e in plugin.runTimeValidity().entries()]

    assert 214 not in codes
    assert check.ignore == ("LIG",), "and it is LIG that a dictionary will describe"


def test_no_clash_when_the_pipeline_mints_nothing(tmp_path, model_with_drg):
    """With a dictionary supplied or no ligand at all, no code is being minted."""
    for ligand_as in ("DICT", "NONE"):
        plugin, check = _prepared(tmp_path, ligand_as=ligand_as, xyzin=model_with_drg)
        check.complaint = ""
        codes = [e["code"] for e in plugin.runTimeValidity().entries()]
        assert 214 not in codes, ligand_as


def test_the_chosen_code_is_what_acedrg_is_told(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert plugin._generatedLigandCode() == "DRG"
    plugin.container.controlParameters.LIGAND_CODE.set("lig")
    assert plugin._generatedLigandCode() == "LIG", "codes are upper case"


def test_an_unreadable_model_is_not_read_as_containing_nothing(tmp_path):
    plugin = _plugin("SubstituteLigand", tmp_path)
    assert plugin._residueCodesIn(tmp_path / "no_such_file.pdb") is None
