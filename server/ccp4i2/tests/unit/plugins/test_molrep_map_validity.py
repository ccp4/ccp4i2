"""molrep_map offers a phaser engine it cannot yet run: say so while editing.

engines.place() raises NotImplementedError for phaser, which used to surface
only once the job had started. validity() now flags the choice on the field.
"""
from ccp4i2.core.tasks import get_plugin_class


def _codes(engine, tmp_path):
    plugin = get_plugin_class("molrep_map")(workDirectory=str(tmp_path), name="molrep_map")
    plugin.container.controlParameters.ENGINE.set(engine)
    return [e.get("code") for e in plugin.validity().getErrors()]


def test_phaser_engine_is_flagged(tmp_path):
    assert 207 in _codes("phaser", tmp_path)


def test_molrep_engine_is_not(tmp_path):
    assert 207 not in _codes("molrep", tmp_path)
