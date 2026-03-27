#!/usr/bin/env python3
"""
Test if direct __dict__ assignment works for CData serialization.
"""
import xml.etree.ElementTree as ET

from ccp4i2.core.CCP4ModelData import CAsuContentSeq
from ccp4i2.core.base_object.fundamental_types import CString, CInt


def test_dict_assignment_serialization():
    """Test that __dict__ assignment produces valid XML serialization."""
    item = CAsuContentSeq(name="test_item")

    # Try direct __dict__ assignment for 'name' and 'polymerType'
    item.__dict__['name'] = CString("BETA", name="name")
    item.__dict__['polymerType'] = CString("PROTEIN", name="polymerType")

    # Also set some regular attributes
    item.sequence = "HPETLVKVK"
    item.nCopies = 1
    item.description = "Test protein"

    # Serialize to XML
    etree = item.getEtree()
    xml_str = ET.tostring(etree, encoding='unicode')

    # Verify XML has content
    assert len(xml_str) > 0, "XML output should not be empty"

    # Check that child elements exist in the XML
    children = list(etree)
    assert len(children) > 0, "XML should have child elements after serialization"
