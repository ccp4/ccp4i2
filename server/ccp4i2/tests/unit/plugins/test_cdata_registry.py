"""
A task registers its own CData classes; core/ knows nothing about them.

core/cdata_registry.py is the one place a def.xml class name is resolved, for
the def.xml handler and the file digest alike. Its sources are the core
implementation modules and every task's ``Task.dataTypes``. The PanDDA types
are the precedent: they live in their task directories and core/ has no
PanDDA module. CCP4-free.
"""
import importlib
import logging

import pytest

from ccp4i2.core import cdata_registry
from ccp4i2.core.base_object.fundamental_types import CString, CList
from ccp4i2.core.tasks import TASKS, Task, locate_def_xml
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser


@pytest.fixture(autouse=True)
def fresh_registry():
    cdata_registry.reset_cache()
    yield
    cdata_registry.reset_cache()


def test_core_has_no_pandda_module_and_the_registry_still_knows_the_types():
    for gone in ("ccp4i2.core.CPanddaEvent", "ccp4i2.core.CPanddaDataset"):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(gone)
    classes = cdata_registry.cdata_classes()
    assert classes["CPanddaEvent"].__module__ == "ccp4i2.wrappers.pandda_events.script.pandda_events_types"
    assert classes["CPanddaDataset"].__module__ == "ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_types"
    for name in ("CPanddaReceiptPerformance", "CPanddaRunPerformance",
                 "CPanddaFanoutPerformance", "CPanddaManifestDataFile"):
        assert name in classes, name
    # And the core types are all still there.
    for name in ("CInt", "CString", "CList", "CPdbDataFile", "CObsDataFile", "CCell", "CDmDomain",
                 "CMoorhenSceneDataFile"):
        assert name in classes, name


def test_the_tasks_declare_the_modules_they_need():
    assert TASKS["pandda_events"].dataTypes == ("ccp4i2.wrappers.pandda_events.script.pandda_events_types",)
    assert TASKS["pandda_campaign"].dataTypes == ("ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_types",)
    assert TASKS["pandda_fanout"].dataTypes == TASKS["pandda_campaign"].dataTypes
    assert "ccp4i2.wrappers.pandda_events.script.pandda_events_types" in list(cdata_registry.task_data_type_modules())


@pytest.mark.parametrize("task, path, item_class", [
    ("pandda_events", "outputData.EVENTS", "CPanddaEvent"),
    ("pandda_campaign", "inputData.DATASETS", "CPanddaDataset"),
])
def test_the_def_xml_resolves_the_list_item_class(task, path, item_class):
    """The silent failure: an unresolved subItem makes CString items with no
    warning. Assert the item class, not just that the def.xml parses."""
    container = DefXmlParser().parse_def_xml(locate_def_xml(task))
    obj = container
    for part in path.split("."):
        obj = getattr(obj, part)
    assert isinstance(obj, CList)
    sub = obj.get_qualifier("subItem")
    assert sub and sub["class"].__name__ == item_class, sub
    assert sub["class"] is not CString


def test_a_broken_data_types_module_is_reported_and_skipped(monkeypatch, caplog):
    monkeypatch.setitem(TASKS, "zz_broken", Task(
        title="broken", pluginPath="x:y", defXmlPath="x",
        dataTypes=("ccp4i2.no_such_module_zz", "ccp4i2.wrappers.pandda_events.script.pandda_events_types")))
    with caplog.at_level(logging.WARNING, logger="ccp4i2:ccp4i2.core.cdata_registry"):
        classes = cdata_registry.cdata_classes()
    assert "ccp4i2.no_such_module_zz" in caplog.text
    assert "CPanddaEvent" in classes and "CPdbDataFile" in classes


def test_the_file_digest_resolves_through_the_same_registry():
    from ccp4i2.lib.utils.files import digest
    assert digest._class_named("CPdbDataFile") is cdata_registry.cdata_classes()["CPdbDataFile"]
    assert digest._class_named("CNoSuchFile") is None
