"""find() searches past string parameters.

A CString's own find() is str.find, and its -1 ("not in this string") was
taken as the search's answer: a name searched for from the root stopped at
the first string parameter on the way. ProvideSequence's SEQIN, behind
SEQUENCETEXT, could not be found, so an upload to it failed (an agent
trial, 2026-10-03).
"""
from ccp4i2.core.tasks import get_plugin_class


def test_a_parameter_behind_a_string_is_found(tmp_path):
    plugin = get_plugin_class("ProvideSequence")(workDirectory=str(tmp_path), name="seq")
    container = plugin.container
    found = container.find("SEQIN")
    assert found is container.inputData.SEQIN
    assert container.find_by_path("inputData.SEQIN") is container.inputData.SEQIN
    assert container.find("NO_SUCH_PARAMETER") is None
