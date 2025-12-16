#!/usr/bin/env python3
"""
Test CAsuContentSeqList serialization to verify XML export works correctly.

This tests whether populated CAsuContentSeq items in a list properly serialize to XML.
"""

import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList


def test_clist_serialization():
    """Test that CAsuContentSeqList properly serializes populated items."""

    print("\n" + "="*80)
    print("Testing CAsuContentSeqList Serialization")
    print("="*80)

    # Create a CAsuContentSeqList
    seq_list = CAsuContentSeqList(name="ASU_CONTENT")
    print(f"\n1. Created CAsuContentSeqList: {type(seq_list).__name__}")
    print(f"   Initial length: {len(seq_list)}")

    # Create first CAsuContentSeq item (BETA)
    item1 = seq_list.makeItem()
    print(f"\n2. Created first item: {type(item1).__name__}")

    # Populate first item
    item1.sequence = "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
    item1.nCopies = 1
    item1.name = "BETA"
    item1.polymerType = "PROTEIN"
    item1.description = "Beta-lactamase"

    print(f"   Set attributes on item1:")
    print(f"     sequence: {item1.sequence.value if hasattr(item1.sequence, 'value') else item1.sequence}")
    print(f"     nCopies: {item1.nCopies.value if hasattr(item1.nCopies, 'value') else item1.nCopies}")
    print(f"     name: {item1.name.value if hasattr(item1.name, 'value') else item1.name}")
    print(f"     polymerType: {item1.polymerType.value if hasattr(item1.polymerType, 'value') else item1.polymerType}")
    print(f"     description: {item1.description.value if hasattr(item1.description, 'value') else item1.description}")

    # Add first item to list
    seq_list.append(item1)
    print(f"\n3. Appended item1 to list. List length: {len(seq_list)}")

    # Create second CAsuContentSeq item (BLIP)
    item2 = seq_list.makeItem()
    print(f"\n4. Created second item: {type(item2).__name__}")

    # Populate second item
    item2.sequence = "AGVMTGAKFTQIQFGMTRQQVLDIAGAENCETGGSFGDSIHCRGHAAGDYYAYATFGFTSAAADAKVDSKSQEKLLAPSAPTLTLAKFNQVTVGMTRAQVLATVGQGSCTTWSEYYPAYPSTAGVTLSLSCFDVDGYSSTGFYRGSAHLWFTDGVLQGKRQWDLV"
    item2.nCopies = 1
    item2.name = "BLIP"
    item2.polymerType = "PROTEIN"
    item2.description = "Beta-lactamase inhibitory protein"

    print(f"   Set attributes on item2:")
    print(f"     name: {item2.name.value if hasattr(item2.name, 'value') else item2.name}")
    print(f"     polymerType: {item2.polymerType.value if hasattr(item2.polymerType, 'value') else item2.polymerType}")
    print(f"     nCopies: {item2.nCopies.value if hasattr(item2.nCopies, 'value') else item2.nCopies}")

    # Add second item to list
    seq_list.append(item2)
    print(f"\n5. Appended item2 to list. List length: {len(seq_list)}")

    # Verify items are in the list
    print(f"\n6. Verifying items in list:")
    for i, item in enumerate(seq_list):
        print(f"   Item {i}: {type(item).__name__}")
        if hasattr(item, 'name'):
            name_val = item.name.value if hasattr(item.name, 'value') else item.name
            print(f"     name: {name_val}")

    # Now serialize to XML
    print(f"\n7. Calling getEtree() to serialize to XML...")
    try:
        etree = seq_list.getEtree()
        print(f"   ✓ getEtree() succeeded, returned: {type(etree)}")

        # Convert to string for display
        xml_str = ET.tostring(etree, encoding='unicode')
        print(f"\n8. XML Output:")
        print("-" * 80)
        # Pretty print
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ")
        # Skip XML declaration line
        pretty_lines = pretty_xml.split('\n')[1:]
        print('\n'.join(pretty_lines))
        print("-" * 80)

        # Check if items have content
        print(f"\n9. Analyzing XML structure:")
        items = etree.findall('.//CAsuContentSeq')
        print(f"   Found {len(items)} CAsuContentSeq elements")

        for i, item_elem in enumerate(items):
            print(f"\n   Item {i}:")
            children = list(item_elem)
            if len(children) == 0:
                print(f"     ❌ EMPTY - No child elements!")
            else:
                print(f"     ✓ Has {len(children)} child elements:")
                for child in children:
                    text = (child.text or '').strip()
                    if text:
                        display_text = text[:50] + "..." if len(text) > 50 else text
                        print(f"       - {child.tag}: {display_text}")
                    else:
                        print(f"       - {child.tag}: (empty)")

        # Final verdict
        print(f"\n" + "="*80)
        if all(len(list(item_elem)) > 0 for item_elem in items):
            print("✅ SUCCESS: All items have populated attributes in XML")
            return True
        else:
            print("❌ FAILURE: Some items are empty in XML")
            return False

    except Exception as e:
        print(f"   ❌ getEtree() failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_clist_serialization()
    sys.exit(0 if success else 1)
