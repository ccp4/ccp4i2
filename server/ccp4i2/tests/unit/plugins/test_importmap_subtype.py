"""ImportMap lets the user tag the kind of map being imported.

A CCP4/MRC map carries no self-description of its role, so ImportMap offers a
``MAP_SUBTYPE`` choice (normal / difference / anomalous / mask / half map) and
stamps it onto the imported file's ``subType``. Half maps in particular need this
so they are recognised by servalcat / the cryo-EM placement task as one of a pair.

CCP4-free: the def.xml parser and the plugin's finalize_output hook, no binaries.
"""

from ccp4i2.core.CCP4XtalData import CMapDataFile
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
from ccp4i2.core.tasks import locate_def_xml
from ccp4i2.wrappers.ImportMap.script.ImportMap import ImportMap


def _container():
    return DefXmlParser().parse_def_xml(locate_def_xml("ImportMap"))


def test_map_subtype_enum_offers_every_kind_including_half_map():
    c = _container()
    p = c.controlParameters.MAP_SUBTYPE
    assert p.qualifiers("enumerators") == [1, 2, 3, 4, 5]
    assert int(p) == CMapDataFile.SUBTYPE_NORMAL          # sensible default
    menu = p.qualifiers("menuText")
    assert len(menu) == 5
    assert any("half map" in m.lower() for m in menu)


def test_finalize_output_stamps_chosen_subtype():
    c = _container()
    c.controlParameters.MAP_SUBTYPE.set(CMapDataFile.SUBTYPE_HALFMAP)
    # Exercise the real hook without spinning up the plugin machinery.
    plugin = object.__new__(ImportMap)
    plugin.container = c
    plugin.finalize_output(c.outputData.MAPOUT)
    assert int(c.outputData.MAPOUT.subType) == CMapDataFile.SUBTYPE_HALFMAP


def test_finalize_output_respects_mask_choice():
    c = _container()
    c.controlParameters.MAP_SUBTYPE.set(CMapDataFile.SUBTYPE_MASK)
    plugin = object.__new__(ImportMap)
    plugin.container = c
    plugin.finalize_output(c.outputData.MAPOUT)
    assert int(c.outputData.MAPOUT.subType) == CMapDataFile.SUBTYPE_MASK
