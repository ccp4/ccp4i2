"""The SAD pipeline says before it runs what its later steps need.

Parrot and ModelCraft take the AU contents only from ASUFILE, and ModelCraft
needs the free set; without them the job failed after Phaser and both hands'
density modification had run (found reviewing the agent judgement)."""
from ccp4i2.core.tasks import get_plugin_class


def _codes(plugin):
    return [r["code"] for r in plugin.validity()._reports]


def test_density_modification_needs_the_au_contents_file(tmp_path):
    plugin = get_plugin_class("phaser_ep_phil")(workDirectory=str(tmp_path), name="ep")
    assert bool(plugin.container.controlParameters.RUNPARROT)  # on by default
    assert 115 in _codes(plugin)
    plugin.container.controlParameters.RUNPARROT.set(False)
    assert 115 not in _codes(plugin)


def test_building_needs_the_free_set(tmp_path):
    plugin = get_plugin_class("phaser_ep_phil")(workDirectory=str(tmp_path), name="ep")
    assert 116 not in _codes(plugin)  # building is off by default
    plugin.container.controlParameters.RUNMODELCRAFT.set(True)
    assert 116 in _codes(plugin)
