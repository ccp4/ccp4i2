"""Make Ligand's "a sketch" says before Run when it cannot work.

The sketch option runs Coot's Lidia as a desktop program. Current CCP4
(Coot 1) ships no lidia executable, so the option could only fail after Run.
Validity now reports it on the field while the job is being set up.

CCP4-free: nothing is run; whether Lidia is installed is pinned both ways.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class


def _codes(error):
    return [e.get("code") for e in error.getErrors()]


def _plugin(tmp_path, mode):
    plugin = get_plugin_class("LidiaAcedrgNew")(
        workDirectory=str(tmp_path), name="LidiaAcedrgNew")
    plugin.container.inputData.MOLSMILESORSKETCH.set(mode)
    plugin.container.inputData.TLC.set("LIG")
    return plugin


def test_sketch_without_lidia_is_an_error_on_the_field(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, "SKETCH")
    monkeypatch.setattr(type(plugin), "_lidiaAvailable", staticmethod(lambda: False))
    error = plugin.validity()
    ours = [e for e in error.getErrors() if e.get("code") == 203]
    assert ours, _codes(error)
    assert ours[0]["name"] == "LidiaAcedrgNew.container.inputData.MOLSMILESORSKETCH"


def test_sketch_with_lidia_is_not_flagged(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, "SKETCH")
    monkeypatch.setattr(type(plugin), "_lidiaAvailable", staticmethod(lambda: True))
    assert 203 not in _codes(plugin.validity())


def test_smiles_is_never_flagged(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, "SMILES")
    monkeypatch.setattr(type(plugin), "_lidiaAvailable", staticmethod(lambda: False))
    assert 203 not in _codes(plugin.validity())
