"""
Test prosmart_refmac plugin inheritance from refmac_i2.

This test verifies that prosmart_refmac properly inherits parameters
like NCYCLES from its parent plugin refmac_i2 via the <file> mechanism.
"""

import sys
import os
from pathlib import Path

import pytest

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Set CCP4I2_ROOT for plugin discovery
if "CCP4I2_ROOT" not in os.environ:
    os.environ["CCP4I2_ROOT"] = str(Path(__file__).parent.parent)

from ccp4i2.core.task_manager.def_xml_handler import parse_def_xml_file


def test_refmac_has_ncycles():
    """Test that refmac_i2 has NCYCLES parameter."""
    refmac_def = Path(__file__).parent.parent / "wrappers/refmac_i2/script/refmac.def.xml"

    if not refmac_def.exists():
        pytest.skip(f"refmac.def.xml not found at {refmac_def}")

    task = parse_def_xml_file(refmac_def)

    # Check that controlParameters exists
    assert hasattr(task, 'controlParameters'), "refmac should have controlParameters"

    # Check that NCYCLES exists
    assert hasattr(task.controlParameters, 'NCYCLES'), "refmac.controlParameters should have NCYCLES"

    # Check default value
    assert task.controlParameters.NCYCLES.value == 10, "NCYCLES default should be 10"

    # Check qualifiers
    assert task.controlParameters.NCYCLES.get_qualifier('min') == 0, "NCYCLES min should be 0"

    print(f"✓ refmac has NCYCLES with default={task.controlParameters.NCYCLES.value}")


def test_prosmart_refmac_inherits_ncycles():
    """Test that prosmart_refmac inherits NCYCLES from refmac_i2."""
    prosmart_def = Path(__file__).parent.parent / "pipelines/prosmart_refmac/script/prosmart_refmac.def.xml"

    if not prosmart_def.exists():
        pytest.skip(f"prosmart_refmac.def.xml not found at {prosmart_def}")

    task = parse_def_xml_file(prosmart_def)

    # Check that controlParameters exists
    assert hasattr(task, 'controlParameters'), \
        "prosmart_refmac should have controlParameters (inherited from refmac)"

    # Check that NCYCLES exists
    assert hasattr(task.controlParameters, 'NCYCLES'), \
        "prosmart_refmac.controlParameters should have NCYCLES (inherited from refmac)"

    # Check default value
    assert task.controlParameters.NCYCLES.value == 10, \
        "Inherited NCYCLES default should be 10"

    # Check qualifiers
    assert task.controlParameters.NCYCLES.get_qualifier('min') == 0, \
        "Inherited NCYCLES min should be 0"

    print(f"✓ prosmart_refmac inherited NCYCLES with default={task.controlParameters.NCYCLES.value}")


def test_prosmart_refmac_file_reference():
    """Test that prosmart_refmac has the correct <file> reference to refmac."""
    import xml.etree.ElementTree as ET

    prosmart_def = Path(__file__).parent.parent / "pipelines/prosmart_refmac/script/prosmart_refmac.def.xml"

    if not prosmart_def.exists():
        pytest.skip(f"prosmart_refmac.def.xml not found at {prosmart_def}")

    tree = ET.parse(prosmart_def)
    root = tree.getroot()

    # Find <file> element
    file_elem = root.find(".//{http://www.ccp4.ac.uk/ccp4ns}file")
    if file_elem is None:
        file_elem = root.find(".//file")

    assert file_elem is not None, "prosmart_refmac should have a <file> element"

    # Check it references refmac
    base_name_elem = file_elem.find(".//baseName")
    assert base_name_elem is not None, "<file> should have baseName element"
    assert base_name_elem.text == "refmac.def.xml", \
        "prosmart_refmac should reference refmac.def.xml"

    print(f"✓ prosmart_refmac correctly references {base_name_elem.text}")


def test_load_nested_xml_expansion():
    """Test that load_nested_xml properly expands <file> references."""
    import xml.etree.ElementTree as ET

    # Add server path for importing
    server_path = Path(__file__).parent.parent / 'server'
    if str(server_path) not in sys.path:
        sys.path.insert(0, str(server_path))

    from ccp4x.lib.utils.parameters.load_xml import load_nested_xml

    prosmart_def = Path(__file__).parent.parent / "pipelines/prosmart_refmac/script/prosmart_refmac.def.xml"

    if not prosmart_def.exists():
        pytest.skip(f"prosmart_refmac.def.xml not found at {prosmart_def}")

    # Parse the original XML
    tree = ET.parse(prosmart_def)
    root = tree.getroot()

    # Check that <file> exists before expansion
    file_elem = root.find(".//file")
    assert file_elem is not None, "Original XML should have <file> element"

    # Expand with load_nested_xml
    expanded_root = load_nested_xml(root)

    # After expansion, should have NCYCLES somewhere in the tree
    ncycles_found = False
    for elem in expanded_root.iter():
        if elem.get('id') == 'NCYCLES':
            ncycles_found = True
            break

    assert ncycles_found, "Expanded XML should contain NCYCLES (inherited from refmac)"

    print("✓ load_nested_xml successfully expanded file references")


if __name__ == "__main__":
    # Run tests
    test_refmac_has_ncycles()
    test_prosmart_refmac_file_reference()
    test_load_nested_xml_expansion()
    test_prosmart_refmac_inherits_ncycles()
    print("\n✅ All inheritance tests passed!")
