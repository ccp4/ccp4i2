"""Crank2 and shelx refuse FREE "existing" without a free-R set.

With FREE "existing" and no FREERFLAG, SetFree() passes nothing and REFMAC
refines with no free reflections ("Switching off use of free R"), so the
build has no R-free and its judgement could not cross-check it. An agent set
FREE "existing" on warning 203's advice and gave no set (HypF, 2026-10-13).
"""
from ccp4i2.core.tasks import get_plugin_class


def _codes(plugin):
    return [r["code"] for r in plugin.validity()._reports]


def test_existing_without_a_set_is_an_error(tmp_path):
    for task in ("crank2", "shelx"):
        plugin = get_plugin_class(task)(workDirectory=str(tmp_path), name="c")
        plugin.container.inputData.FREE.set("existing")
        assert 205 in _codes(plugin), task
        plugin.container.inputData.FREERFLAG.setFullPath(str(tmp_path / "free.mtz"))
        assert 205 not in _codes(plugin), task


def test_new_and_no_are_unaffected(tmp_path):
    plugin = get_plugin_class("crank2")(workDirectory=str(tmp_path), name="c")
    for free in ("new", "no"):
        plugin.container.inputData.FREE.set(free)
        assert 205 not in _codes(plugin)
