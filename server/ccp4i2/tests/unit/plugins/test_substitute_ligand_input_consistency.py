"""The mode menus and the inputs they name have to agree.

SubstituteLigand is driven by two menus. 'Observed data provided as' says
whether the reflections arrive merged or unmerged; 'Chemistry of the ligand'
says whether the molecule to be placed arrives as a Mol file, a dictionary or
a SMILES string. Neither menu carries the input it names, and --- because any
one of them is optional depending on the other menu --- no input is required
on its own. So nothing stopped the two disagreeing: 'unmerged' with the
unmerged file list empty, 'SMILES String' with no string in the box.

What happened then depended on which way they disagreed, and none of it was
legible. Unmerged-with-nothing reached aimless_pipe with an empty file list.
SMILES-with-nothing reached acedrg with nothing to build. A batch fragment
campaign submitted these fifty at a time and got fifty failures, each one
apparently about a different program.

These tests cover the refusal, at the Run dialog, where the menu and the empty
field are both on screen. Alongside it: coordinates that hold no atoms --- the
shape a truncated download leaves behind --- which are refused for the same
reason, that the pipeline is built on placing them.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class


def _plugin(tmp_path, obs_as=None, ligand_as=None):
    plugin_class = get_plugin_class("SubstituteLigand")
    assert plugin_class is not None, "SubstituteLigand is not in the registry"
    directory = tmp_path / "SubstituteLigand"
    directory.mkdir(exist_ok=True)
    plugin = plugin_class(workDirectory=str(directory), name="SubstituteLigand")
    if obs_as is not None:
        plugin.container.controlParameters.OBSAS.set(obs_as)
    if ligand_as is not None:
        plugin.container.controlParameters.LIGANDAS.set(ligand_as)
    return plugin


def _model(tmp_path, name="start.pdb", atoms=True):
    path = tmp_path / name
    path.write_text(
        ("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
         if atoms else "HEADER    NOTHING HERE                            01-JAN-70   XXXX\n")
        + "END\n",
        encoding="utf-8",
    )
    return path


def _file(tmp_path, name):
    path = tmp_path / name
    path.write_text("not really, but it exists\n", encoding="utf-8")
    return path


def _sweep(plugin, tmp_path, name="sweep.mtz"):
    """One unmerged sweep on the pipeline's list."""
    unmerged = plugin.container.inputData.UNMERGEDFILES
    unmerged.addItem()
    unmerged[-1].file.setFullPath(str(_file(tmp_path, name)))
    return unmerged[-1]


def _entries(plugin, code):
    return [e for e in plugin.runTimeValidity().entries() if e["code"] == code]


def _codes(plugin):
    return [e["code"] for e in plugin.runTimeValidity().entries()]


# --- observed data: the form chosen must be the form supplied ---------------

def test_unmerged_with_no_unmerged_file_is_refused(tmp_path):
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="NONE")
    assert 216 in _codes(plugin), "aimless would have been handed an empty list"


def test_the_unmerged_refusal_points_at_the_empty_list(tmp_path):
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="NONE")
    entry, = _entries(plugin, 216)
    assert entry["name"] == "SubstituteLigand.container.inputData.UNMERGEDFILES"
    assert entry["severity"] >= 4, "there is nothing to merge, so nothing to run"


def test_the_unmerged_refusal_offers_both_ways_out(tmp_path):
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="NONE")
    entry, = _entries(plugin, 216)
    assert "unmerged" in entry["details"] and "merged" in entry["details"], (
        "supplying the file and changing the menu are both answers")


def test_merged_with_no_merged_file_is_refused(tmp_path):
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="NONE")
    entry, = _entries(plugin, 216)
    assert entry["name"] == "SubstituteLigand.container.inputData.F_SIGF_IN"


def test_unmerged_files_satisfy_the_unmerged_mode(tmp_path):
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="NONE")
    _sweep(plugin, tmp_path)
    assert 216 not in _codes(plugin)


def test_a_merged_file_satisfies_the_merged_mode(tmp_path):
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="NONE")
    plugin.container.inputData.F_SIGF_IN.setFullPath(str(_file(tmp_path, "obs.mtz")))
    assert 216 not in _codes(plugin)


def test_an_unmerged_file_does_not_satisfy_the_merged_mode(tmp_path):
    """The wrong menu with the right-looking file is the whole failure mode."""
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="NONE")
    _sweep(plugin, tmp_path)
    assert 216 in _codes(plugin)


def test_merged_data_without_a_free_r_set_is_advised_not_refused(tmp_path):
    """Nothing makes one on this route --- aimless only runs on unmerged data."""
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="NONE")
    plugin.container.inputData.F_SIGF_IN.setFullPath(str(_file(tmp_path, "obs.mtz")))
    entry, = _entries(plugin, 219)
    assert entry["severity"] == 2, "the job runs, it just cannot report R-free"
    assert entry["name"] == "SubstituteLigand.container.inputData.FREERFLAG_IN"


def test_no_free_r_advisory_once_one_is_given(tmp_path):
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="NONE")
    plugin.container.inputData.F_SIGF_IN.setFullPath(str(_file(tmp_path, "obs.mtz")))
    plugin.container.inputData.FREERFLAG_IN.setFullPath(str(_file(tmp_path, "free.mtz")))
    assert 219 not in _codes(plugin)


def test_no_free_r_advisory_on_the_unmerged_route(tmp_path):
    """aimless makes the set there, so its absence at input means nothing."""
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="NONE")
    _sweep(plugin, tmp_path)
    assert 219 not in _codes(plugin)


# --- ligand chemistry: the mechanism chosen must be the one supplied --------

@pytest.mark.parametrize("mode,field", [
    ("MOL", "MOLIN"),
    ("DICT", "DICTIN"),
    ("SMILES", "SMILESIN"),
])
def test_a_chemistry_mode_with_nothing_in_it_is_refused(tmp_path, mode, field):
    plugin = _plugin(tmp_path, ligand_as=mode)
    entry, = _entries(plugin, 217)
    assert entry["name"] == f"SubstituteLigand.container.inputData.{field}"
    assert entry["severity"] >= 4


def test_the_chemistry_refusal_names_the_menu_entry_the_user_chose(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="MOL")
    entry, = _entries(plugin, 217)
    assert "MDL Mol file" in entry["details"], (
        "the message has to name what the user actually sees in the menu")


def test_the_chemistry_refusal_offers_none_as_a_way_out(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="SMILES")
    entry, = _entries(plugin, 217)
    assert "NONE" in entry["details"], (
        "running with no ligand at all is a legitimate answer and is easy to miss")


@pytest.mark.parametrize("mode,field,value", [
    ("MOL", "MOLIN", "ligand.mol"),
    ("DICT", "DICTIN", "ligand.cif"),
])
def test_supplying_what_the_mode_reads_clears_it(tmp_path, mode, field, value):
    plugin = _plugin(tmp_path, ligand_as=mode)
    getattr(plugin.container.inputData, field).setFullPath(str(_file(tmp_path, value)))
    assert 217 not in _codes(plugin)


def test_supplying_the_smiles_string_clears_it(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="SMILES")
    plugin.container.inputData.SMILESIN.set("CCO")
    assert 217 not in _codes(plugin)


def test_no_ligand_needs_no_chemistry(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    assert 217 not in _codes(plugin)


def test_a_sketch_needs_no_chemistry_supplied(tmp_path):
    """SKETCH draws the molecule; there is no file whose absence means anything."""
    plugin = _plugin(tmp_path, ligand_as="SKETCH")
    assert 217 not in _codes(plugin)


def test_a_mol_file_in_smiles_mode_is_flagged_as_ignored(tmp_path):
    """Silently ignored input reads exactly like input that was used."""
    plugin = _plugin(tmp_path, ligand_as="SMILES")
    plugin.container.inputData.SMILESIN.set("CCO")
    plugin.container.inputData.MOLIN.setFullPath(str(_file(tmp_path, "ligand.mol")))
    entry, = _entries(plugin, 218)
    assert entry["severity"] == 2, "the job is runnable, it just ignores the file"
    assert entry["name"] == "SubstituteLigand.container.inputData.MOLIN"


def test_the_ignored_input_message_says_which_menu_entry_would_use_it(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="SMILES")
    plugin.container.inputData.SMILESIN.set("CCO")
    plugin.container.inputData.DICTIN.setFullPath(str(_file(tmp_path, "ligand.cif")))
    entry, = _entries(plugin, 218)
    assert "REFMAC Dict" in entry["details"]
    assert "REFMAC dictionary" in entry["details"], (
        "acronyms must survive the sentence they are put in")


def test_everything_is_ignored_when_no_ligand_is_being_fitted(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    plugin.container.inputData.MOLIN.setFullPath(str(_file(tmp_path, "ligand.mol")))
    entry, = _entries(plugin, 218)
    assert "no ligand is to be fitted" in entry["details"]


def test_what_the_mode_does_read_is_not_called_ignored(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="MOL")
    plugin.container.inputData.MOLIN.setFullPath(str(_file(tmp_path, "ligand.mol")))
    assert 218 not in _codes(plugin)


# --- coordinates: there has to be something to place ------------------------

def test_coordinates_that_are_not_set_are_refused(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    entries = [e for e in plugin.runTimeValidity().entries()
               if e["name"].endswith(".XYZIN") and e["severity"] >= 4]
    assert entries, "the pipeline has nothing to place and must not start"


def test_coordinates_naming_a_file_that_is_not_there_are_refused(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    plugin.container.inputData.XYZIN.setFullPath(str(tmp_path / "never_arrived.pdb"))
    entries = [e for e in plugin.runTimeValidity().entries()
               if e["name"].endswith(".XYZIN") and e["severity"] >= 4]
    assert entries


def test_coordinates_holding_no_atoms_are_refused(tmp_path):
    """A truncated download parses cleanly and contains nothing."""
    plugin = _plugin(tmp_path, ligand_as="NONE")
    plugin.container.inputData.XYZIN.setFullPath(str(_model(tmp_path, atoms=False)))
    entry, = _entries(plugin, 215)
    assert entry["severity"] >= 4
    assert entry["name"] == "SubstituteLigand.container.inputData.XYZIN"


def test_the_empty_model_message_does_not_blame_the_file_format(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    plugin.container.inputData.XYZIN.setFullPath(str(_model(tmp_path, atoms=False)))
    entry, = _entries(plugin, 215)
    assert "no coordinates" in entry["details"], (
        "the file was read, so 'unreadable' would send the user looking in "
        "the wrong place")


def test_coordinates_with_atoms_are_accepted(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    plugin.container.inputData.XYZIN.setFullPath(str(_model(tmp_path)))
    assert 215 not in _codes(plugin)


def test_an_unreadable_model_is_refused_rather_than_waved_through(tmp_path):
    plugin = _plugin(tmp_path, ligand_as="NONE")
    assert plugin._atomCountIn(tmp_path / "no_such_file.pdb") is None, (
        "0 would read as 'checked, and empty' rather than 'never checked'")


# --- the checks are independent of one another ------------------------------

def test_a_missing_model_does_not_hide_a_disagreeing_menu(tmp_path):
    """One problem per submission is how a batch of fifty took fifty attempts."""
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="SMILES")
    codes = _codes(plugin)
    assert 216 in codes, "the merged file is missing and should be said so now"
    assert 217 in codes, "so is the SMILES string"


def test_both_menus_are_reported_at_once(tmp_path):
    plugin = _plugin(tmp_path, obs_as="UNMERGED", ligand_as="MOL")
    plugin.container.inputData.XYZIN.setFullPath(str(_model(tmp_path)))
    codes = _codes(plugin)
    assert 216 in codes and 217 in codes


def test_a_complete_job_is_not_refused(tmp_path):
    plugin = _plugin(tmp_path, obs_as="MERGED", ligand_as="SMILES")
    plugin.container.inputData.XYZIN.setFullPath(str(_model(tmp_path)))
    plugin.container.inputData.F_SIGF_IN.setFullPath(str(_file(tmp_path, "obs.mtz")))
    plugin.container.inputData.FREERFLAG_IN.setFullPath(str(_file(tmp_path, "free.mtz")))
    plugin.container.inputData.SMILESIN.set("CCO")
    assert plugin.runTimeValidity().maxSeverity() < 4
