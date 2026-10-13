"""The first listMinLength items of a list are required, however empty.

An unset optional object downgrades its children's errors to warnings, and
an empty list item looked like one: aimless_pipe's UNMERGEDFILES[0] reported
its missing file as a warning, run_job let the job start, and it failed at
once, seven times over for an agent (2026-10-13). Items beyond the minimum
keep the lenient rule. Saving is unaffected: CList.getEtree writes every
item whatever its state, and nothing that saves consults validity.
"""
import tempfile

import pytest

from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core.base_object.error_reporting import SEVERITY_ERROR, SEVERITY_WARNING
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.tasks import TASKS, get_plugin_class


def _plugin(task):
    return get_plugin_class(task)(workDirectory=tempfile.mkdtemp(), name=task)


def _child(obj, name):
    return next(c for c in obj.children() if c.objectName() == name)


def _unmerged_files():
    """The plugin and its UNMERGEDFILES: the plugin is kept, since one
    collected empties the container it made (destroy() on GC)."""
    plugin = _plugin("aimless_pipe")
    return plugin, _child(_child(plugin.container, "inputData"), "UNMERGEDFILES")


def _severities(report, suffix):
    return [int(e["severity"]) for e in report.getErrors() if str(e.get("name", "")).endswith(suffix)]


def test_a_required_item_missing_its_file_blocks():
    plugin, files = _unmerged_files()
    assert _severities(files[0].validity(), "UNMERGEDFILES[0].file") == [SEVERITY_ERROR]


def test_an_item_beyond_the_minimum_stays_lenient():
    plugin, files = _unmerged_files()
    files.addItem()
    assert _severities(files[1].validity(), "UNMERGEDFILES[1].file") == [SEVERITY_WARNING]


def _required_items(container):
    """(list path, item) for each of the first listMinLength items of every
    list with a minimum, anywhere in the container."""
    for child in container.children():
        name = child.objectName()
        if not name or name.startswith("_"):
            continue
        if isinstance(child, CList):
            try:
                minimum = int(child.get_qualifier("listMinLength") or 0)
            except (TypeError, ValueError):
                minimum = 0
            for item in list(child)[:minimum]:
                yield name, item
        elif hasattr(child, "children") and not isinstance(child, CDataFile):
            yield from _required_items(child)


@pytest.mark.parametrize("task", sorted(TASKS))
def test_no_required_item_reports_a_missing_required_field_as_a_warning(task):
    """The guard: anywhere a def.xml gives a list a minimum and starts it
    with items, those items' unset required fields block."""
    try:
        plugin = _plugin(task)
    except Exception as err:  # a task needing what this environment lacks
        pytest.skip(f"{task} cannot be built here: {type(err).__name__}")
    for list_name, item in _required_items(plugin.container):
        if isinstance(item, CDataFile) or not hasattr(item, "children"):
            continue  # a bare file item never had the downgrade
        report = item.validity()
        for field in item.children():
            if (field.objectName() and field.get_qualifier("allowUndefined") is False
                    and not field.isSet()):
                found = _severities(report, f".{field.objectName()}")
                assert SEVERITY_WARNING not in found, (
                    f"{task}: {list_name} item's required {field.objectName()} "
                    f"is only a warning")
