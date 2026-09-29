"""A file inside a CList is recorded under the list's name, not a bare index.

`CList` names its elements by index alone -- `item.rename(f"[{i}]")` -- so a
file inside one answers `"[0]"` to `objectName()`. Recording that bare loses
which list the file went into, and two file lists on one job then record
`"[0]"` twice, indistinguishably. A live database had 162 such FileUse rows.

The one place that tried to put the name back together could never fire:

    item.objectName() if item.objectName() else f"{child.objectName()}[{i}]"

the fallback is dead, because CList guarantees the element has a name and it is
truthy.
"""
import pytest

from ccp4i2.core.base_object.fundamental_types import (
    CList,
    qualified_parameter_name,
)
from ccp4i2.core.CCP4Container import CContainer


class _Named:
    """Minimal stand-in: a name, and a parent that may be a CList."""

    def __init__(self, name, parent=None):
        self._name = name
        self._parent = parent

    def objectName(self):
        return self._name

    def parent(self):
        return self._parent


class TestTheRule:
    def test_a_plain_parameter_is_its_own_name(self):
        assert qualified_parameter_name(_Named("XYZIN")) == "XYZIN"

    def test_an_unnamed_object_stays_unnamed(self):
        assert qualified_parameter_name(_Named("")) == ""

    def test_a_list_element_takes_the_list_name(self):
        container = CContainer(name="inputData")
        dict_list = CList(name="DICT_LIST", parent=container)
        item = _Named("[0]", parent=dict_list)
        assert qualified_parameter_name(item) == "DICT_LIST[0]"

    def test_two_lists_are_told_apart(self):
        """The whole point: bare '[0]' from two lists is indistinguishable."""
        container = CContainer(name="inputData")
        first = CList(name="DICT_LIST", parent=container)
        second = CList(name="REFERENCE_MODELS", parent=container)

        assert qualified_parameter_name(_Named("[0]", parent=first)) == "DICT_LIST[0]"
        assert (
            qualified_parameter_name(_Named("[0]", parent=second))
            == "REFERENCE_MODELS[0]"
        )

    def test_an_index_outside_a_list_is_left_alone(self):
        """Only a CList parent supplies a name; anything else is not guessed at."""
        container = CContainer(name="inputData")
        assert qualified_parameter_name(_Named("[0]", parent=container)) == "[0]"

    def test_no_parent_is_survivable(self):
        assert qualified_parameter_name(_Named("[0]")) == "[0]"


class TestAgainstARealTask:
    @pytest.fixture(autouse=True)
    def _needs_plugins(self):
        pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

    def test_servalcat_dict_list_elements_are_qualified(self, tmp_path):
        from ccp4i2.core.tasks import get_plugin_class

        plugin = get_plugin_class("servalcat_pipe")(
            workDirectory=str(tmp_path), parent=None
        )
        dict_list = plugin.container.inputData.DICT_LIST
        dict_list.append(dict_list.makeItem())
        dict_list.append(dict_list.makeItem())

        names = [qualified_parameter_name(item) for item in dict_list]
        assert names == ["DICT_LIST[0]", "DICT_LIST[1]"]

        # And the plugin's own descendant walk agrees, which is what feeds the
        # error messages and used to emit a bare index.
        found = [
            name
            for name, _ in plugin._find_datafile_descendants(
                plugin.container.inputData
            )
        ]
        assert "DICT_LIST[0]" in found, found
        assert not [n for n in found if n.startswith("[")], found
