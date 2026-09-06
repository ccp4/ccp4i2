"""A follow-on Phaser step takes the search models and composition from
the step before it; the user's own entries are never overwritten."""
from pathlib import Path

import pytest

pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

import ccp4i2
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.lib.utils.parameters.set_input_by_context import inherit_inputs_from_context

GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma"


def make(task, tmp_path):
    return get_plugin_class(task)(workDirectory=str(tmp_path / task), name=task)


def frf_with_models(tmp_path):
    frf = make("phaser_mr_frf_phil", tmp_path)
    inp = frf.container.inputData
    inp.ENSEMBLES.append(inp.ENSEMBLES.makeItem())
    inp.ENSEMBLES[-1].label.set("beta")
    inp.ENSEMBLES[-1].number.set(1)
    item = inp.ENSEMBLES[-1].pdbItemList.makeItem()
    inp.ENSEMBLES[-1].pdbItemList.append(item)
    item.structure.setFullPath(str(GAMMA / "gamma_model.pdb"))
    item.identity_to_target.set(0.9)
    inp.COMP_BY.set("ASU")
    inp.ASUFILE.setFullPath(str(GAMMA / "gamma.asu.xml"))
    return frf


def test_translation_function_inherits_the_rotation_functions_models(tmp_path):
    frf = frf_with_models(tmp_path)
    ftf = make("phaser_mr_ftf_phil", tmp_path)
    copied = inherit_inputs_from_context(ftf, frf)
    inp = ftf.container.inputData
    assert "ENSEMBLES" in copied and "ASUFILE" in copied
    assert [str(e.label) for e in inp.ENSEMBLES] == ["beta"]
    assert str(inp.ENSEMBLES[0].pdbItemList[0].structure.fullPath) == str(GAMMA / "gamma_model.pdb")
    assert str(inp.COMP_BY) == "ASU" and str(inp.ASUFILE.fullPath) == str(GAMMA / "gamma.asu.xml")


def test_what_the_user_set_is_kept(tmp_path):
    frf = frf_with_models(tmp_path)
    ftf = make("phaser_mr_ftf_phil", tmp_path)
    inp = ftf.container.inputData
    inp.ENSEMBLES.append(inp.ENSEMBLES.makeItem())
    inp.ENSEMBLES[-1].label.set("mine")
    copied = inherit_inputs_from_context(ftf, frf)
    assert "ENSEMBLES" not in copied
    assert [str(e.label) for e in inp.ENSEMBLES] == ["mine"]


def test_a_task_that_declares_nothing_inherits_nothing(tmp_path):
    frf = frf_with_models(tmp_path)
    other = get_plugin_class("pointless_reindexToMatch")(workDirectory=str(tmp_path / "p"))
    assert inherit_inputs_from_context(other, frf) == []


def test_the_pipeline_shares_the_declaration(tmp_path):
    frf = frf_with_models(tmp_path)
    pipeline = make("phaser_pipeline_phil", tmp_path)
    assert "ENSEMBLES" in inherit_inputs_from_context(pipeline, frf)
