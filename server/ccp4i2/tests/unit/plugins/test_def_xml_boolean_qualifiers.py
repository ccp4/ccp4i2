"""<onlyEnumerators>1</onlyEnumerators> means true, to the interface as well.

37 def.xml qualifiers spell a boolean as 1. Kept as the string "1" it was
true in Python but not to the interface, which tests === true and so gave
Sculptor's pruning and B-factor choices, and Chainsaw's mode, free text
boxes instead of menus.
"""
from ccp4i2.core.tasks import get_plugin_class


def test_one_is_true_for_boolean_qualifiers(tmp_path):
    plugin = get_plugin_class("sculptor")(workDirectory=str(tmp_path))
    ctrl = plugin.container.controlParameters
    for name in ("PRUNING", "BFACTOR", "DELETION"):
        assert getattr(ctrl, name).get_qualifier("onlyEnumerators") is True, name


def test_target_sequence_defaults_to_the_first(tmp_path):
    """A numeric default is still a number, and the target menu starts on
    the first sequence, as the wrappers assume when it is unset."""
    for task in ("chainsaw", "sculptor"):
        plugin = get_plugin_class(task)(workDirectory=str(tmp_path))
        assert plugin.container.controlParameters.TARGETINDEX.value == 0, task
