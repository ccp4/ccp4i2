"""A value of the wrong shape is refused when set, not later at validation.

An agent set an ASU_CONTENT item's sequence to {"text": "..."}; it was
stored as the dict's repr and failed only at validation.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.lib.utils.parameters.set_parameter import shape_errors


@pytest.fixture
def asu_content(tmp_path):
    plugin = get_plugin_class("ProvideAsuContents")(workDirectory=str(tmp_path), name="p")
    return plugin.container.inputData.ASU_CONTENT


def test_an_object_for_a_text_field_is_refused(asu_content):
    errors = shape_errors(asu_content, [{"name": "A", "sequence": {"text": "MKV"}, "nCopies": 1}],
                          "inputData.ASU_CONTENT")
    assert len(errors) == 1 and errors[0].startswith(
        "inputData.ASU_CONTENT[0].sequence takes a single value")


def test_a_list_takes_a_list(asu_content):
    assert shape_errors(asu_content, {"name": "A"}, "inputData.ASU_CONTENT") == [
        "inputData.ASU_CONTENT is a list: give a list of items"]


def test_the_right_shape_passes(asu_content):
    assert shape_errors(asu_content, [{"name": "A", "sequence": "MKV", "nCopies": 2}],
                        "inputData.ASU_CONTENT") == []
