#!/usr/bin/env python3
"""
Test CAsuContentSeqList serialization using i2run-style construction.
This mimics how i2run creates items and uses CList.set() to copy them.
"""
import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList
from ccp4i2.core.base_object.fundamental_types import CString, CInt

print("\n" + "="*80)
print("Testing CAsuContentSeqList with i2run-style Construction")
print("="*80)

# Create source list (like inputData.ASU_CONTENT in i2run)
source_list = CAsuContentSeqList(name="ASU_CONTENT")
print(f"\n1. Created source list: {type(source_list).__name__}")

# Create first item using setattr (like i2run does)
item1 = CAsuContentSeq()
setattr(item1, 'sequence', 'HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW')
setattr(item1, 'nCopies', 1)
setattr(item1, 'name', 'BETA')
setattr(item1, 'polymerType', 'PROTEIN')
setattr(item1, 'description', 'Beta-lactamase')
print(f"\n2. Created item1 using setattr")
print(f"   item1.name = {item1.name.value if hasattr(item1.name, 'value') else item1.name}")
print(f"   item1.nCopies = {item1.nCopies.value if hasattr(item1.nCopies, 'value') else item1.nCopies}")

# Append to source list
source_list.append(item1)
print(f"\n3. Appended item1 to source_list. Length: {len(source_list)}")

# Create second item
item2 = CAsuContentSeq()
setattr(item2, 'name', 'BLIP')
setattr(item2, 'sequence', 'AGVMTGAKFTQIQFGMTRQQVLDIAGAENCETGGSFGDSIHCRGHAAGDYYAYATFGFTSAAADAKVDSKSQEKLLAPSAPTLTLAKFNQVTVGMTRAQVLATVGQGSCTTWSEYYPAYPSTAGVTLSLSCFDVDGYSSTGFYRGSAHLWFTDGVLQGKRQWDLV')
setattr(item2, 'polymerType', 'PROTEIN')
setattr(item2, 'nCopies', 1)
setattr(item2, 'description', 'Beta-lactamase inhibitory protein')
print(f"\n4. Created item2 using setattr")
print(f"   item2.name = {item2.name.value if hasattr(item2.name, 'value') else item2.name}")

# Append to source list
source_list.append(item2)
print(f"\n5. Appended item2 to source_list. Length: {len(source_list)}")

# Now create target list and use .set() to copy (like ProvideAsuContents does)
target_list = CAsuContentSeqList(name="seqList")
print(f"\n6. Created target list: {type(target_list).__name__}")
print(f"   Calling target_list.set(source_list)...")

target_list.set(source_list)
print(f"   ✓ set() completed. Target list length: {len(target_list)}")

# Verify items in target list
print(f"\n7. Verifying items in target_list:")
for i, item in enumerate(target_list):
    print(f"   Item {i}: {type(item).__name__}")
    if hasattr(item, 'name'):
        name_val = item.name.value if hasattr(item.name, 'value') else item.name
        print(f"     name: {name_val}")
    if hasattr(item, 'sequence'):
        seq_val = item.sequence.value if hasattr(item.sequence, 'value') else item.sequence
        print(f"     sequence length: {len(seq_val) if seq_val else 0}")

# Serialize target list to XML
print(f"\n8. Serializing target_list to XML...")
etree = target_list.getEtree()
xml_str = ET.tostring(etree, encoding='unicode')

# Pretty print
from xml.dom import minidom
dom = minidom.parseString(xml_str)
pretty_xml = dom.toprettyxml(indent="  ")
print("-" * 80)
print(pretty_xml)
print("-" * 80)

# Check if items have content
print(f"\n9. Analyzing XML structure:")
items = list(etree)
print(f"   Found {len(items)} CAsuContentSeq elements")

all_populated = True
for i, item_elem in enumerate(items):
    children = list(item_elem)
    print(f"\n   Item {i}:")
    if len(children) == 0:
        print(f"     ❌ EMPTY - no child elements!")
        all_populated = False
    else:
        print(f"     ✓ Has {len(children)} child elements:")
        for child in children:
            text = child.text[:50] + "..." if child.text and len(child.text) > 50 else child.text
            print(f"       - {child.tag}: {text or '(empty)'}")

print("\n" + "="*80)
if all_populated:
    print("✅ SUCCESS: All items have populated attributes in XML")
else:
    print("❌ FAILURE: Some items are empty in XML")
print("="*80 + "\n")

# Return success/failure
sys.exit(0 if all_populated else 1)
