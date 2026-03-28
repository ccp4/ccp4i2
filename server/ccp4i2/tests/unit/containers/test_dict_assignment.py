#!/usr/bin/env python3
# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
