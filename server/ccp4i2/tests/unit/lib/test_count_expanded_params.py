"""
Count how many parameters are in expanded controlParameters XML.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

from ccp4i2.core import CCP4Utils
from ccp4i2.lib.utils.parameters.load_xml import load_nested_xml


def test_count_params_in_expanded_xml():
    """Count parameters in expanded controlParameters."""

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

    # Find controlParameters
    for elem in expanded_root.iter():
        if elem.get('id') == 'controlParameters':
            # Count direct content children
            content_children = elem.findall('./content[@id]')
            print(f"\n=== controlParameters in expanded XML ===")
            print(f"Total <content> children: {len(content_children)}")

            # Show all parameter names
            param_names = [c.get('id') for c in content_children]
            param_names.sort()

            print(f"\nAll parameters:")
            for name in param_names:
                print(f"  {name}")

            # Check for NCYCLES
            if 'NCYCLES' in param_names:
                print(f"\n✅ NCYCLES is in the list")
            else:
                print(f"\n❌ NCYCLES is NOT in the list")

            break


if __name__ == "__main__":
    test_count_params_in_expanded_xml()
