#!/usr/bin/env python3
"""
Test if direct __dict__ assignment works for CData serialization.
"""
import xml.etree.ElementTree as ET
from xml.dom import minidom

from ccp4i2.core.CCP4ModelData import CAsuContentSeq
from ccp4i2.core.base_object.fundamental_types import CString, CInt

print("\n" + "="*80)
print("Testing __dict__ Assignment for CData Serialization")
print("="*80)

# Create a CAsuContentSeq item
item = CAsuContentSeq(name="test_item")
print(f"\n1. Created CAsuContentSeq")

# Try direct __dict__ assignment for 'name' and 'polymerType'
print(f"\n2. Attempting __dict__ assignment...")
item.__dict__['name'] = CString("BETA", name="name")
item.__dict__['polymerType'] = CString("PROTEIN", name="polymerType")
print(f"   __dict__['name'] = CString('BETA')")
print(f"   __dict__['polymerType'] = CString('PROTEIN')")

# Also set some regular attributes
item.sequence = "HPETLVKVK"
item.nCopies = 1
item.description = "Test protein"

# Check what's in __dict__
print(f"\n3. Checking __dict__ contents:")
for key, val in item.__dict__.items():
    if not key.startswith('_') and key not in ['parent', 'children', 'signals']:
        print(f"   {key}: {type(val).__name__} = {val if not hasattr(val, 'value') else val.value}")

# Try to serialize
print(f"\n4. Serializing to XML...")
etree = item.getEtree()
xml_str = ET.tostring(etree, encoding='unicode')

# Pretty print
dom = minidom.parseString(xml_str)
pretty_xml = dom.toprettyxml(indent="  ")
print(pretty_xml)

# Check what made it into XML
print(f"\n5. Analyzing XML output...")
for child in etree:
    print(f"   {child.tag}: {child.text}")

print("="*80)
