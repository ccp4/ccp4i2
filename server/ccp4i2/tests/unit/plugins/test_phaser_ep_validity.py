"""phaser_EP wrappers must validate without an F_OR_I parameter.

They subclass phaser_MR, whose validity() reads inputData.F_OR_I to advise on
intensities versus amplitudes. EP declares no such parameter -- it always
phases on F+/F- pairs -- so the inherited check raised AttributeError at
runTimeValidity and every phaser_EP job died before running.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class


@pytest.mark.parametrize("task", ["phaser_EP_AUTO", "phaser_EP_LLG"])
def test_phaser_ep_validity_does_not_need_f_or_i(task, tmp_path):
    plugin = get_plugin_class(task)(workDirectory=str(tmp_path), name=task)
    assert not hasattr(plugin.container.inputData, "F_OR_I")
    plugin.validity()   # raised AttributeError before the override
