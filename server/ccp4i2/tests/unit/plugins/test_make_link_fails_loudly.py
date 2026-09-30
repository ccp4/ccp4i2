"""MakeLink must not finish successfully holding no model.

AceDRG writes the link dictionary; applying that link to the user's model is
our own gemmi code, and every way it could fail used to be a bare `return` or
a swallowed exception, leaving a job that said Finished and produced nothing
(or, worse, an empty two-line PDB annotated "Model with links applied").

Reported by Martin Maly against CCP4 9.0.015; the same holes were here.
"""
import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.pipelines.MakeLink.script.MakeLink import _has_atoms

gemmi = pytest.importorskip("gemmi", reason="_has_atoms works on gemmi structures")


@pytest.fixture
def plugin():
    # Kept alive by the fixture: a destroyed plugin empties the container.
    return get_plugin_class("MakeLink")(parent=None, name="mk")


def codes(error, severity=None):
    """The error codes in a CErrorReport, optionally only at one severity."""
    return [item["code"] for item in error.entries()
            if severity is None or item["severity"] == severity]


# --- _has_atoms: the empty-structure trap -------------------------------

def test_has_atoms_is_false_for_a_structure_gemmi_could_not_fill(tmp_path):
    # gemmi returns an EMPTY structure for junk rather than raising, which is
    # how a garbage input file became a valid-looking "model with links".
    junk = tmp_path / "junk.pdb"
    junk.write_text("this is not a model file at all\nnonsense\n")
    assert _has_atoms(gemmi.read_structure(str(junk))) is False


def test_has_atoms_is_true_for_a_real_model():
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "LYS"
    atom = gemmi.Atom()
    atom.name = "NZ"
    residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    assert _has_atoms(structure) is True


# --- validity(): the two ways to ask for a model and not get one --------

def test_a_model_without_the_toggle_warns_but_does_not_block(plugin, tmp_path):
    model = tmp_path / "in.pdb"
    model.write_text("END\n")
    plugin.container.inputData.XYZIN.setFullPath(str(model))
    plugin.container.controlParameters.TOGGLE_LINK = False
    error = plugin.validity()
    assert 306 in codes(error, CCP4ErrorHandling.SEVERITY_WARNING)
    # Advisory only: the dictionary is still a perfectly good result.
    assert 306 not in codes(error, CCP4ErrorHandling.SEVERITY_ERROR)


def test_the_toggle_without_a_model_is_an_error(plugin):
    plugin.container.controlParameters.TOGGLE_LINK = True
    error = plugin.validity()
    assert 307 in codes(error, CCP4ErrorHandling.SEVERITY_ERROR)


def test_toggle_and_model_together_raise_neither(plugin, tmp_path):
    model = tmp_path / "in.pdb"
    model.write_text("END\n")
    plugin.container.inputData.XYZIN.setFullPath(str(model))
    plugin.container.controlParameters.TOGGLE_LINK = True
    error = plugin.validity()
    assert 306 not in codes(error)
    assert 307 not in codes(error)


# --- applyLinksToModel(): the guard paths return a status ---------------

def test_no_toggle_is_success_not_failure(plugin):
    # Nothing was asked for, so the job is still a good job.
    plugin.container.controlParameters.TOGGLE_LINK = False
    assert plugin.applyLinksToModel(1.5) == CPluginScript.SUCCEEDED


def test_toggle_without_a_model_fails(plugin):
    plugin.container.controlParameters.TOGGLE_LINK = True
    assert plugin.applyLinksToModel(1.5) == CPluginScript.FAILED


def test_no_link_in_the_dictionary_fails(plugin, tmp_path):
    # get_link_bond_value() returned None: it could not find the link it just
    # asked AceDRG to make. That used to be a silent skip.
    model = tmp_path / "in.pdb"
    model.write_text("END\n")
    plugin.container.inputData.XYZIN.setFullPath(str(model))
    plugin.container.controlParameters.TOGGLE_LINK = True
    assert plugin.applyLinksToModel(None) == CPluginScript.FAILED


def test_an_unreadable_model_fails_rather_than_writing_an_empty_one(plugin, tmp_path):
    junk = tmp_path / "junk.pdb"
    junk.write_text("this is not a model file at all\nnonsense\n")
    plugin.container.inputData.XYZIN.setFullPath(str(junk))
    plugin.container.controlParameters.TOGGLE_LINK = True
    plugin.container.inputData.RES_NAME_1_TLC.set("LYS")
    plugin.container.inputData.RES_NAME_2_TLC.set("PLP")
    plugin.container.inputData.ATOM_NAME_1.set("NZ")
    plugin.container.inputData.ATOM_NAME_2.set("C4A")
    assert plugin.applyLinksToModel(1.5) == CPluginScript.FAILED
    assert not plugin.container.outputData.XYZOUT.isSet()


# --- residue codes: "Lys" is LYS -----------------------------------------
# AceDRG finds LYS.cif for "Lys", then looks inside it for a comp "Lys" and
# stops. The task interface upper-cases its own lookup, so the atom dropdown
# filled and nothing warned before the job failed (BAZ2B campaign, 2026-09-29).

def test_library_codes_are_upper_cased_before_the_instruction(plugin):
    inp = plugin.container.inputData
    inp.RES_NAME_1_TLC.set("Lys")
    inp.RES_NAME_2_TLC.set(" glu ")
    inp.ATOM_NAME_1.set("NZ")
    inp.ATOM_NAME_2.set("CD")
    plugin.normaliseResidueCodes()
    assert str(inp.RES_NAME_1_TLC) == "LYS"
    assert str(inp.RES_NAME_2_TLC) == "GLU"
    instruct = plugin.createLinkInstruction()
    assert "RES-NAME-1 LYS " in instruct
    assert "RES-NAME-2 GLU " in instruct


def test_a_dictionary_code_is_left_as_the_dictionary_spells_it(plugin):
    # CIF mode names a comp in the user's own file, which may be lower case.
    inp = plugin.container.inputData
    inp.MON_1_TYPE.set("CIF")
    inp.RES_NAME_1_CIF.set("Lig")
    plugin.normaliseResidueCodes()
    assert str(inp.RES_NAME_1_CIF) == "Lig"


# --- the edits as a description ------------------------------------------
# DELETE_ATOMS_n / BOND_ORDERS_n / CHARGES_n hold the modified monomer; the
# plugin writes AceDRG's words (test_make_link_instruction.py pins those).

def test_the_edit_lists_hold_their_own_types(plugin):
    # An unresolved subItem silently becomes CString (docs/cdata.md).
    inp = plugin.container.inputData
    assert type(inp.DELETE_ATOMS_1.makeItem()).__name__ == "CString"
    assert type(inp.BOND_ORDERS_1.makeItem()).__name__ == "CMakeLinkBondOrder"
    assert type(inp.CHARGES_2.makeItem()).__name__ == "CMakeLinkCharge"


def _set_link(plugin):
    inp = plugin.container.inputData
    inp.RES_NAME_1_TLC.set("LYS")
    inp.RES_NAME_2_TLC.set("GLU")
    inp.ATOM_NAME_1.set("NZ")
    inp.ATOM_NAME_2.set("CD")


def test_the_instruction_is_written_from_the_lists(plugin):
    _set_link(plugin)
    c = plugin.container
    c.set_parameter("inputData.DELETE_ATOMS_2", ["OE2", "OXT"], skip_first=True)
    c.set_parameter("inputData.BOND_ORDERS_2",
                    [{"ATOM_1": "CD", "ATOM_2": "OE1", "ORDER": "DOUBLE"}], skip_first=True)
    c.set_parameter("inputData.CHARGES_1", [{"ATOM": "NZ", "CHARGE": 0}], skip_first=True)
    instruct = plugin.createLinkInstruction()
    assert instruct.endswith(
        "CHANGE CHARGE 1 NZ 0"
        " DELETE ATOM OE2 2 DELETE ATOM OXT 2 CHANGE BOND CD OE1 DOUBLE 2")


def test_the_single_edit_fields_of_older_jobs_are_folded_in(plugin):
    # An older job, a clone of one, or an i2run command written for them.
    inp = plugin.container.inputData
    inp.TOGGLE_DELETE_2 = True
    inp.DELETE_2.set("OE2")
    inp.TOGGLE_CHANGE_2 = True
    inp.CHANGE_BOND_2.set("CD -- OE1")
    inp.CHANGE_2_TYPE.set("DOUBLE")
    inp.TOGGLE_CHARGE_1 = True
    inp.CHARGE_1.set("NZ")
    inp.CHARGE_1_VALUE.set(0)
    assert plugin.monomerEdits(2).deletes == ["OE2"]
    assert plugin.monomerEdits(2).bond_orders == [("CD", "OE1", "DOUBLE")]
    assert plugin.monomerEdits(1).charges == [("NZ", 0)]


def test_an_unticked_single_edit_is_not_folded_in(plugin):
    inp = plugin.container.inputData
    inp.TOGGLE_DELETE_2 = False
    inp.DELETE_2.set("OE2")
    assert plugin.monomerEdits(2).is_empty()


def test_a_contradictory_description_blocks_the_job(plugin):
    _set_link(plugin)
    plugin.container.set_parameter("inputData.DELETE_ATOMS_1", ["NZ"], skip_first=True)
    error = plugin.validity()
    assert 308 in codes(error, CCP4ErrorHandling.SEVERITY_ERROR)
    assert plugin.createLinkInstruction() == CPluginScript.FAILED


def test_extra_instructions_that_would_hang_acedrg_block_the_job(plugin):
    _set_link(plugin)
    plugin.container.controlParameters.EXTRA_ACEDRG_INSTRUCTIONS.set("CHANGE OE1 DOUBLE 2")
    assert 309 in codes(plugin.validity(), CCP4ErrorHandling.SEVERITY_ERROR)
    assert plugin.createLinkInstruction() == CPluginScript.FAILED


def test_the_default_extra_instructions_are_only_comments(plugin):
    # The box ships with a syntax crib; it must not itself be an error.
    _set_link(plugin)
    assert 309 not in codes(plugin.validity())
