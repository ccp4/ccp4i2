"""
Check if NCYCLES appears in the expanded prosmart_refmac XML.
"""

import sys
import os
import xml.etree.ElementTree as ET
from pathlib import Path

# Add project root and server to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "server"))

# Set CCP4I2_ROOT for plugin discovery
os.environ["CCP4I2_ROOT"] = str(project_root)

from ccp4x.lib.utils.parameters.load_xml import load_nested_xml


def test_ncycles_in_expanded_xml():
    """Test that NCYCLES appears in the expanded XML after load_nested_xml."""

    prosmart_def = project_root / "pipelines/prosmart_refmac/script/prosmart_refmac.def.xml"

    if not prosmart_def.exists():
        print(f"❌ prosmart_refmac.def.xml not found at {prosmart_def}")
        return

    # Parse the original XML
    tree = ET.parse(prosmart_def)
    root = tree.getroot()

    # Expand with load_nested_xml
    print("Expanding <file> references with load_nested_xml...")
    expanded_root = load_nested_xml(root)

    # Search for NCYCLES in the expanded XML
    ncycles_elements = []
    for elem in expanded_root.iter():
        if elem.get('id') == 'NCYCLES':
            ncycles_elements.append(elem)

    print(f"\nFound {len(ncycles_elements)} elements with id='NCYCLES'")

    if ncycles_elements:
        for i, elem in enumerate(ncycles_elements):
            print(f"\n=== NCYCLES element #{i+1} ===")
            print(f"Tag: {elem.tag}")
            print(f"Attributes: {elem.attrib}")

            # Find parent containers
            parent_chain = []
            current = elem
            while current is not None:
                if current.get('id'):
                    parent_chain.append(current.get('id'))
                current = find_parent(expanded_root, current)
            print(f"Parent chain: {' -> '.join(reversed(parent_chain))}")

            # Show className and qualifiers
            class_name = elem.find('.//className')
            if class_name is not None:
                print(f"className: {class_name.text}")

            qualifiers = elem.find('.//qualifiers')
            if qualifiers is not None:
                print(f"qualifiers:")
                for q in qualifiers:
                    print(f"  {q.tag}: {q.text}")

        print(f"\n✓ NCYCLES found in expanded XML")
    else:
        print(f"\n✗ NCYCLES not found in expanded XML!")

        # Show what parameters are in controlParameters
        print("\n=== Parameters in controlParameters ===")
        for elem in expanded_root.iter():
            if elem.get('id') == 'controlParameters':
                for content in elem.findall('.//content'):
                    content_id = content.get('id')
                    if content_id:
                        print(f"  {content_id}")
                break


def find_parent(root, target):
    """Find parent of target element in the tree."""
    for elem in root.iter():
        if target in list(elem):
            return elem
    return None


if __name__ == "__main__":
    test_ncycles_in_expanded_xml()
