"""
Count ALL controlParameters containers in expanded XML.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

from ccp4i2.core import CCP4Utils
from ccp4i2.lib.utils.parameters.load_xml import load_nested_xml


def test_count_all_control_params():
    """Find all controlParameters containers in expanded XML."""

    prosmart_def = Path(CCP4Utils.getCCP4I2Dir()) / "pipelines/prosmart_refmac/script/prosmart_refmac.def.xml"

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
