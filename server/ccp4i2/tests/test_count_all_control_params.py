"""
Count ALL controlParameters containers in expanded XML.
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

from ccp4i2.lib.utils.parameters.load_xml import load_nested_xml


def test_count_all_control_params():
    """Find all controlParameters containers in expanded XML."""

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

    # Find ALL controlParameters
    count = 0
    for elem in expanded_root.iter():
        if elem.get('id') == 'controlParameters':
            count += 1
            # Count direct content children
            content_children = elem.findall('./content[@id]')
            print(f"\n=== controlParameters #{count} ===")
            print(f"Total <content> children: {len(content_children)}")

            # Show first 10 and check for NCYCLES
            param_names = [c.get('id') for c in content_children]
            has_ncycles = 'NCYCLES' in param_names

            print(f"First 10 parameters: {', '.join(param_names[:10])}")
            print(f"Has NCYCLES: {has_ncycles}")

    print(f"\n=== Total controlParameters containers: {count} ===")


if __name__ == "__main__":
    test_count_all_control_params()
