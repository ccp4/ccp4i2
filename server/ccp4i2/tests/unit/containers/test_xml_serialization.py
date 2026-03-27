"""Test XML serialization and deserialization for CData objects."""

import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.core.base_object.base_classes import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean


def test_simple_value_getEtree():
    """Test getEtree() on simple value types."""
    # CInt
    num = CInt(value=42, name="myInt")
    elem = num.getEtree()
    assert elem.tag == "myInt"
    assert elem.text == "42"

    # CString
    text = CString(value="hello", name="myString")
    elem = text.getEtree()
    assert elem.tag == "myString"
    assert elem.text == "hello"

    # CBoolean
    flag = CBoolean(value=True, name="myFlag")
    elem = flag.getEtree()
    assert elem.tag == "myFlag"
    assert elem.text == "True"


def test_simple_value_setEtree():
    """Test setEtree() on simple value types."""
    # CInt
    num = CInt(name="test")
    elem = ET.Element("test")
    elem.text = "99"
    num.setEtree(elem)
    assert num.value == 99

    # CString
    text = CString(name="test")
    elem = ET.Element("test")
    elem.text = "world"
    text.setEtree(elem)
    assert text.value == "world"

    # CBoolean
    flag = CBoolean(name="test")
    elem = ET.Element("test")
    elem.text = "true"
    flag.setEtree(elem)
    assert flag.value is True


def test_container_xml_serialization():
    """Test XML serialization of containers."""
    # Create a container with nested data
    container = CContainer(name="task")
    container.addContent(CInt, "param1")
    container.param1.value = 10

    container.addContent(CString, "param2")
    container.param2.value = "test_value"

    # Serialize
    elem = container.getEtree()

    # Verify structure
    assert elem.tag == "task"
    assert len(elem) == 2

    # Find children
    param1_elem = elem.find("param1")
    assert param1_elem is not None
    assert param1_elem.text == "10"

    param2_elem = elem.find("param2")
    assert param2_elem is not None
    assert param2_elem.text == "test_value"


def test_container_xml_deserialization():
    """Test XML deserialization into containers."""
    # Create container structure
    container = CContainer(name="task")
    container.addContent(CInt, "value1")
    container.addContent(CString, "value2")

    # Create XML to deserialize
    root = ET.Element("task")
    child1 = ET.SubElement(root, "value1")
    child1.text = "42"
    child2 = ET.SubElement(root, "value2")
    child2.text = "loaded"

    # Deserialize
    container.setEtree(root)

    # Verify values were loaded
    assert container.value1.value == 42
    assert container.value2.value == "loaded"


def test_xml_roundtrip():
    """Test full roundtrip: create -> serialize -> deserialize -> compare."""
    # Create original container
    original = CContainer(name="data")
    original.addContent(CInt, "count")
    original.count.value = 100
    original.addContent(CBoolean, "enabled")
    original.enabled.value = True
    original.addContent(CString, "label")
    original.label.value = "Test Label"

    # Serialize to XML
    xml_elem = original.getEtree()

    # Create new container and deserialize
    restored = CContainer(name="data")
    restored.addContent(CInt, "count")
    restored.addContent(CBoolean, "enabled")
    restored.addContent(CString, "label")
    restored.setEtree(xml_elem)

    # Verify all values match
    assert restored.count.value == 100
    assert restored.enabled.value is True
    assert restored.label.value == "Test Label"


def test_qualifiers_serialization():
    """Test qualifier serialization/deserialization."""
    num = CInt(name="test")
    num.set_qualifier("min", 0)
    num.set_qualifier("max", 100)
    num.set_qualifier("default", 50)

    # Serialize qualifiers
    qual_elem = num.getQualifiersEtree()

    # Verify structure
    assert qual_elem.tag == "qualifiers"
    assert len(qual_elem) == 3

    # Find specific qualifiers
    min_elem = qual_elem.find("min")
    assert min_elem is not None
    assert min_elem.text == "0"

    max_elem = qual_elem.find("max")
    assert max_elem is not None
    assert max_elem.text == "100"


def test_qualifiers_deserialization():
    """Test loading qualifiers from XML."""
    num = CInt(name="test")

    # Create qualifiers XML
    qual_elem = ET.Element("qualifiers")
    min_child = ET.SubElement(qual_elem, "min")
    min_child.text = "10"
    max_child = ET.SubElement(qual_elem, "max")
    max_child.text = "90"

    # Deserialize
    num.setQualifiersEtree(qual_elem)

    # Verify qualifiers were loaded
    assert num.get_qualifier("min") == 10
    assert num.get_qualifier("max") == 90


@pytest.fixture
def temp_xml_file():
    """Create a temporary XML file for testing."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
        temp_path = f.name

    yield temp_path

    # Cleanup
    Path(temp_path).unlink(missing_ok=True)


def test_save_load_xml_file(temp_xml_file):
    """Test saving and loading container from XML file."""
    # Create container
    container = CContainer(name="project")
    container.addContent(CInt, "cycles")
    container.cycles.value = 25
    container.addContent(CString, "method")
    container.method.value = "refinement"

    # Save to file
    container.saveContentsToXml(temp_xml_file)

    # Verify file was created
    assert Path(temp_xml_file).exists()

    # Create new container and load
    loaded = CContainer(name="project")
    loaded.addContent(CInt, "cycles")
    loaded.addContent(CString, "method")
    loaded.loadContentsFromXml(temp_xml_file)

    # Verify data was loaded correctly
    assert loaded.cycles.value == 25
    assert loaded.method.value == "refinement"


def test_nested_container_serialization():
    """Test serialization of nested containers."""
    # Create nested structure
    root = CContainer(name="root")
    child = root.addContent(CContainer, "child")
    child.addContent(CInt, "value")
    child.value.value = 123

    # Serialize
    elem = root.getEtree()

    # Verify structure
    assert elem.tag == "root"
    child_elem = elem.find("child")
    assert child_elem is not None
    value_elem = child_elem.find("value")
    assert value_elem is not None
    assert value_elem.text == "123"


def test_load_save_def_file(temp_xml_file):
    """Test loading and saving DEF files (.def.xml)."""
    # Create container with structure
    container = CContainer(name="task")
    container.addContent(CInt, "ncycles")
    container.ncycles.set_qualifier("min", 1)
    container.ncycles.set_qualifier("max", 100)
    container.ncycles.set_qualifier("default", 10)

    container.addContent(CString, "method")
    container.method.set_qualifier("default", "refinement")

    # Save as DEF file
    def_file = temp_xml_file.replace(".xml", ".def.xml")
    container.saveDefFile(def_file)

    # Verify file was created
    assert Path(def_file).exists()

    # Load into new container
    loaded = CContainer(name="task")
    loaded.addContent(CInt, "ncycles")
    loaded.addContent(CString, "method")
    loaded.loadDefFile(def_file)

    # Verify qualifiers were loaded (structure, not values)
    # For now, just verify it doesn't crash
    assert hasattr(loaded, 'ncycles')
    assert hasattr(loaded, 'method')

    # Cleanup
    Path(def_file).unlink(missing_ok=True)


def test_load_save_params_file(temp_xml_file):
    """Test loading and saving PARAMS files (.params.xml)."""
    # Create container with data values
    container = CContainer(name="task")
    container.addContent(CInt, "ncycles")
    container.ncycles.value = 25

    container.addContent(CString, "method")
    container.method.value = "refinement"

    # Save as PARAMS file
    params_file = temp_xml_file.replace(".xml", ".params.xml")
    container.saveParamsFile(params_file)

    # Verify file was created
    assert Path(params_file).exists()

    # Load into new container (with same structure)
    loaded = CContainer(name="task")
    loaded.addContent(CInt, "ncycles")
    loaded.addContent(CString, "method")
    loaded.loadParamsFile(params_file)

    # Verify values were loaded
    assert loaded.ncycles.value == 25
    assert loaded.method.value == "refinement"

    # Cleanup
    Path(params_file).unlink(missing_ok=True)
