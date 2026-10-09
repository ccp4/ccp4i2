"""Crank2 asks for the anomalous scatterer before it runs.

Without one it defines no substructure, runs its setup, and stops at FA
estimation with "No target specified for process of FA estimation" (an
agent trial, 2026-10-03, could not tell what was missing).
"""
from ccp4i2.core.tasks import get_plugin_class


def _codes(plugin):
    return [r["code"] for r in plugin.validity()._reports]


def test_no_atom_type_is_an_error(tmp_path):
    plugin = get_plugin_class("crank2")(workDirectory=str(tmp_path), name="c")
    assert 202 in _codes(plugin)
    plugin.container.inputData.ATOM_TYPE.set("Xe")
    assert 202 not in _codes(plugin)


def test_not_needed_when_starting_after_phasing(tmp_path):
    plugin = get_plugin_class("crank2")(workDirectory=str(tmp_path), name="c")
    plugin.container.inputData.START_PIPELINE.set("dmfull")
    assert 202 not in _codes(plugin)


# The default FREE "new" makes Crank2 choose its own free set, silently: an
# agent trial caught it only because the guide stresses free sets


def test_a_new_free_set_is_warned_of(tmp_path):
    plugin = get_plugin_class("crank2")(workDirectory=str(tmp_path), name="c")
    assert 203 in _codes(plugin)
    plugin.container.inputData.FREERFLAG.setFullPath(str(tmp_path / "free.mtz"))
    codes = _codes(plugin)
    assert 204 in codes and 203 not in codes  # given, but ignored
    plugin.container.inputData.FREE.set("existing")
    codes = _codes(plugin)
    assert 203 not in codes and 204 not in codes


def test_shelx_warns_too(tmp_path):
    plugin = get_plugin_class("shelx")(workDirectory=str(tmp_path), name="s")
    assert 203 in _codes(plugin)

