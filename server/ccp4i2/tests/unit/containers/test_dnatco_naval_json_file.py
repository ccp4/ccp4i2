"""CDnatcoNavalJsonFile: the file type for DNATCO's NAVAL validation JSON.

A new CDataFile subclass is only usable once four registries agree on it:
the class itself, the mimetype tables the database is seeded from, the
def.xml class resolver, and the digest module's copy of those tables (which
is now the same object). These tests pin all four.
"""
import pytest

from ccp4i2.core.base_object.base_classes import CDataFile
from ccp4i2.core.CCP4ModelData import CDnatcoNavalJsonFile
from ccp4i2.db.ccp4i2_static_data import FILETYPELIST, FILETYPES_CLASS, FILETYPES_TEXT

MIMETYPE = "application/dnatco-naval-json"


def test_class_declares_its_mimetype():
    assert issubclass(CDnatcoNavalJsonFile, CDataFile)
    obj = CDnatcoNavalJsonFile()
    assert obj.qualifiers("mimeTypeName") == MIMETYPE
    assert obj.qualifiers("fileExtensions") == ["json"]


def test_mimetype_tables_agree():
    # FILETYPES_TEXT[i] <-> FILETYPES_CLASS[i] <-> FILETYPELIST id i
    assert len(FILETYPES_TEXT) == len(FILETYPES_CLASS) == len(FILETYPELIST)
    assert [item[0] for item in FILETYPELIST] == list(range(len(FILETYPELIST)))
    index = FILETYPES_TEXT.index(MIMETYPE)
    assert FILETYPES_CLASS[index] == "DnatcoNavalJsonFile"
    assert FILETYPELIST[index][1] == MIMETYPE


def test_class_maps_to_its_mimetype():
    from ccp4i2.lib.cdata_utils import get_file_type_from_class
    assert get_file_type_from_class(CDnatcoNavalJsonFile()) == MIMETYPE


def test_digest_uses_the_seeded_tables():
    pytest.importorskip("gemmi", reason="digest imports gemmi")
    pytest.importorskip("django", reason="digest imports the db models")
    try:
        from ccp4i2.lib.utils.files import digest
    except Exception as exc:  # unconfigured Django settings, etc.
        pytest.skip(f"digest not importable here: {exc}")
    assert digest.FILETYPES_TEXT is FILETYPES_TEXT
    assert digest.FILETYPES_CLASS is FILETYPES_CLASS
    assert "DnatcoNavalJsonFile" in digest.CLASS_REGISTRY or "CDnatcoNavalJsonFile" in digest.CLASS_REGISTRY


def test_def_xml_resolves_the_class_name():
    """A def.xml naming the class must get it, not the silent CString fallback."""
    from ccp4i2.core.tasks import locate_def_xml
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
    for task, item in (("dnatco", "JSONOUT"), ("dnatco_pipe", "JSONOUT1"), ("dnatco_pipe", "JSONOUT2")):
        container = DefXmlParser().parse_def_xml(locate_def_xml(task))
        obj = getattr(container.outputData, item)
        assert isinstance(obj, CDnatcoNavalJsonFile), f"{task}.{item} resolved to {type(obj).__name__}"
