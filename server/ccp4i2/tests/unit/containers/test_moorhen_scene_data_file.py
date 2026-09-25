"""CMoorhenSceneDataFile: a Moorhen scene a job authored about its outputs.

A new CDataFile subclass is only usable once the class, the mimetype tables
the database is seeded from, and the def.xml class resolver agree on it.
These tests pin all three; the last one parses a def.xml that names the
class, because a name the resolver does not know falls back to CString
without a word.
"""
import textwrap

from ccp4i2.core.base_object.base_classes import CDataFile
from ccp4i2.core.CMoorhenSceneDataFile import CMoorhenSceneDataFile
from ccp4i2.db.ccp4i2_static_data import FILETYPELIST, FILETYPES_CLASS, FILETYPES_TEXT

MIMETYPE = "application/moorhen-scene"


def test_class_declares_its_mimetype():
    assert issubclass(CMoorhenSceneDataFile, CDataFile)
    obj = CMoorhenSceneDataFile()
    assert obj.qualifiers("mimeTypeName") == MIMETYPE
    assert "scene.yaml" in obj.qualifiers("fileExtensions")


def test_mimetype_tables_agree():
    assert len(FILETYPES_TEXT) == len(FILETYPES_CLASS) == len(FILETYPELIST)
    assert [item[0] for item in FILETYPELIST] == list(range(len(FILETYPELIST)))
    index = FILETYPES_TEXT.index(MIMETYPE)
    assert FILETYPES_CLASS[index] == "MoorhenSceneDataFile"
    assert FILETYPELIST[index][1] == MIMETYPE


def test_class_maps_to_its_mimetype():
    from ccp4i2.lib.cdata_utils import get_file_type_from_class
    assert get_file_type_from_class(CMoorhenSceneDataFile()) == MIMETYPE


def test_def_xml_resolves_the_class_name(tmp_path):
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

    def_xml = tmp_path / "scene_task.def.xml"
    def_xml.write_text(textwrap.dedent("""\
        <?xml version='1.0' encoding='ASCII'?>
        <ccp4i2>
            <ccp4i2_header>
                <function>DEF</function>
                <pluginName>scene_task</pluginName>
            </ccp4i2_header>
            <ccp4i2_body id="scene_task">
                <container id="outputData">
                    <content id="SCENE">
                        <className>CMoorhenSceneDataFile</className>
                        <qualifiers>
                            <guiLabel>Scene</guiLabel>
                        </qualifiers>
                    </content>
                </container>
            </ccp4i2_body>
        </ccp4i2>
        """))
    container = DefXmlParser().parse_def_xml(str(def_xml))
    obj = container.outputData.SCENE
    assert isinstance(obj, CMoorhenSceneDataFile), f"SCENE resolved to {type(obj).__name__}"
