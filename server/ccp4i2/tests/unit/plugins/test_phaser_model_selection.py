"""Phaser searches with the atoms the user selected, not the whole file.

phaser_simple_phil's XYZIN, and each ensemble model of phaser_pipeline_phil,
take an atom selection (one chain of a downloaded file); until this was
fixed every Phaser job searched with the whole file whatever the selection
said (found while reviewing the sculptor judgement).
"""
import os

import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.CCP4Utils import getCCP4I2Dir  # noqa: E402
from ccp4i2.core.tasks import get_plugin_class  # noqa: E402

MODEL = os.path.join(getCCP4I2Dir(), "demo_data", "CDK1CyclinBCKS2", "1buh.pdb")  # chains A, B


@pytest.fixture
def plugin(tmp_path):
    return get_plugin_class("phaser_mr_auto_phil")(workDirectory=str(tmp_path), name="p")


def _add_model(plugin, selection=None):
    ensembles = plugin.container.inputData.ENSEMBLES
    ensembles.append(ensembles.makeItem())
    item = ensembles[-1].pdbItemList.makeItem()
    ensembles[-1].pdbItemList.append(item)
    item = ensembles[-1].pdbItemList[-1]
    item.structure.setFullPath(MODEL)
    if selection:
        item.structure.selection.text.set(selection)
    return item


def test_a_selected_model_reaches_phaser_as_the_selection(plugin):
    _add_model(plugin, "A/")
    assert plugin._prepare_models() == plugin.SUCCEEDED
    written = plugin._model_paths[(0, 0)]
    chains = {ch.name for ch in gemmi.read_structure(written)[0]}
    assert chains == {"A"}


def test_the_same_file_twice_with_different_selections(plugin):
    _add_model(plugin, "A/")
    _add_model(plugin, "B/")
    assert plugin._prepare_models() == plugin.SUCCEEDED
    assert {ch.name for ch in gemmi.read_structure(plugin._model_paths[(0, 0)])[0]} == {"A"}
    assert {ch.name for ch in gemmi.read_structure(plugin._model_paths[(1, 0)])[0]} == {"B"}


def test_a_model_without_a_selection_is_given_as_it_is(plugin):
    _add_model(plugin)
    assert plugin._prepare_models() == plugin.SUCCEEDED
    assert (0, 0) not in plugin._model_paths


def test_the_shim_reads_the_selected_file(plugin):
    from ccp4i2.wrappers.phaser_phil.script.phaser_shims import EnsembleListShim
    _add_model(plugin, "A/")
    plugin._prepare_models()
    shim = EnsembleListShim("ENSEMBLES", "phaser.ensemble", "phaser.search",
                            path_map=plugin._model_paths)
    text = repr(shim.convert(plugin.container, str(plugin.getWorkDirectory())))
    assert plugin._model_paths[(0, 0)] in text and MODEL not in text


def test_phaser_simple_phil_carries_its_xyzin_selection_into_the_ensemble(tmp_path):
    simple = get_plugin_class("phaser_simple_phil")(workDirectory=str(tmp_path), name="s")
    simple.container.inputData.XYZIN.setFullPath(MODEL)
    simple.container.inputData.XYZIN.selection.text.set("A/")
    simple.createEnsembleElements()
    structure = simple.container.inputData.ENSEMBLES[0].pdbItemList[0].structure
    assert structure.isSelectionSet() and str(structure.selection.text) == "A/"
