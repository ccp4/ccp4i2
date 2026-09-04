"""A class's declared qualifiers are readable without an instance.

``Meta.qualifiers`` replaced the QUALIFIERS class attribute in the declarative
rewrite. Code that still read the attribute -- the follow-on population of
file lists, acedrg's element table -- got nothing back and failed quietly.
``class_qualifier`` is the sanctioned way to read the declaration.
"""
from ccp4i2.core import CCP4ModelData, CCP4XtalData
from ccp4i2.core.base_object.cdata import CData
from ccp4i2.core.base_object.cdata_file import CDataFile


def test_a_file_class_declares_its_mime_type():
    assert CData.class_qualifier("mimeTypeName") is None
    assert CCP4ModelData.CPdbDataFile.class_qualifier("mimeTypeName") == "chemical/x-pdb"
    assert (
        CCP4XtalData.CMapCoeffsDataFile.class_qualifier("mimeTypeName")
        == "application/CCP4-mtz-map"
    )


def test_it_agrees_with_what_an_instance_sees():
    instance = CCP4ModelData.CPdbDataFile()
    assert instance.get_qualifier("mimeTypeName") == CCP4ModelData.CPdbDataFile.class_qualifier(
        "mimeTypeName"
    )


def test_the_default_is_returned_for_an_undeclared_key():
    assert CCP4ModelData.CPdbDataFile.class_qualifier("noSuchQualifier", "fallback") == "fallback"


def test_every_from_previous_job_list_item_class_can_be_looked_up():
    """The 14 tasks whose file lists auto-populate all use these item classes."""
    for cls in (
        CCP4ModelData.CPdbDataFile,
        CCP4XtalData.CMapCoeffsDataFile,
        CCP4ModelData.CDictDataFile,
        CCP4ModelData.CSeqDataFile,
    ):
        assert issubclass(cls, CDataFile)
        assert cls.class_qualifier("mimeTypeName"), cls.__name__


def test_the_element_table_acedrg_reads_is_declared():
    elements = CCP4ModelData.CElement.class_qualifier("enumerators")
    assert elements[:3] == ["H", "He", "Li"]
    assert "Zn" in elements
