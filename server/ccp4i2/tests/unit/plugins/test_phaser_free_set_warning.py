"""MR followed by refinement, with no free-R set: advice, not silence.

The refinement after MR then reports no R-free, which is what judges whether
the placement is right; an agent driving the task was told nothing.
"""
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class


def _codes(error, severity):
    return [r["code"] for r in error._reports if r["severity"] == severity]


def test_no_free_set_is_a_warning(tmp_path):
    plugin = get_plugin_class("phaser_simple_phil")(workDirectory=str(tmp_path), name="mr")
    assert bool(plugin.container.inputData.RUNREFMAC)  # refinement is on by default
    error = plugin.validity()
    assert 220 in _codes(error, CCP4ErrorHandling.SEVERITY_WARNING)


def test_no_warning_without_refinement(tmp_path):
    plugin = get_plugin_class("phaser_simple_phil")(workDirectory=str(tmp_path), name="mr")
    plugin.container.inputData.RUNREFMAC.set(False)
    assert 220 not in [r["code"] for r in plugin.validity()._reports]
