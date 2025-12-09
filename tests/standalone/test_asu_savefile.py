#!/usr/bin/env python3
"""
Test CAsuDataFile.saveFile() after using CList.set()
"""
import sys
from pathlib import Path
import tempfile

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList, CAsuDataFile

print("\n" + "="*80)
print("Testing CAsuDataFile.saveFile() with CList.set()")
print("="*80)

# Create source list with items (like inputData.ASU_CONTENT)
source_list = CAsuContentSeqList(name="ASU_CONTENT")
item1 = CAsuContentSeq()
item1.sequence = 'HPETLVKVK'
item1.nCopies = 1
item1.name = 'BETA'
item1.polymerType = 'PROTEIN'
item1.description = 'Beta-lactamase'
source_list.append(item1)

item2 = CAsuContentSeq()
item2.name = 'BLIP'
item2.sequence = 'AGVMTGAKFT'
item2.polymerType = 'PROTEIN'
item2.nCopies = 1
item2.description = 'BLIP protein'
source_list.append(item2)

print(f"\n1. Created source_list with 2 items")
print(f"   Item 0 name: {source_list[0].name.value if hasattr(source_list[0].name, 'value') else source_list[0].name}")
print(f"   Item 0 isSet('name'): {source_list[0].isSet('name')}")
print(f"   Item 1 name: {source_list[1].name.value if hasattr(source_list[1].name, 'value') else source_list[1].name}")
print(f"   Item 1 isSet('name'): {source_list[1].isSet('name')}")

# Create ASU data file
asu_file = CAsuDataFile()
print(f"\n2. Created CAsuDataFile")

# Set file path
with tempfile.TemporaryDirectory() as tmpdir:
    test_path = Path(tmpdir) / "test.asu.xml"
    asu_file.relPath = str(tmpdir)
    asu_file.baseName = "test.asu.xml"
    print(f"   Set path to: {test_path}")
    
    # Copy items using set() (like ProvideAsuContents does)
    print(f"\n3. Calling fileContent.seqList.set(source_list)...")
    asu_file.fileContent.seqList.set(source_list)
    print(f"   ✓ set() completed. Target list length: {len(asu_file.fileContent.seqList)}")
    
    # Check if items are populated
    print(f"\n4. Checking target list items:")
    for i, item in enumerate(asu_file.fileContent.seqList):
        print(f"   Item {i}:")
        print(f"     name: {item.name.value if hasattr(item.name, 'value') else item.name}")
        print(f"     isSet('name'): {item.isSet('name')}")
        print(f"     isSet('sequence'): {item.isSet('sequence')}")
        print(f"     isSet('nCopies'): {item.isSet('nCopies')}")
    
    # Save file
    print(f"\n5. Calling asu_file.saveFile()...")
    asu_file.saveFile()
    print(f"   ✓ saveFile() completed")
    
    # Read back and check
    print(f"\n6. Reading back saved XML:")
    with open(test_path, 'r') as f:
        content = f.read()
        print(content)
    
    # Check if XML has content
    if '<name>BETA</name>' in content and '<name>BLIP</name>' in content:
        print("\n✅ SUCCESS: XML contains name elements")
        sys.exit(0)
    else:
        print("\n❌ FAILURE: XML missing name elements")
        sys.exit(1)
