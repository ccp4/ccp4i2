#!/usr/bin/env python3
"""
Quick test to verify that CUnmergedDataFileList properly specifies its subItem type.
"""

import sys
import os
from pathlib import Path

# Set up environment using relative path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
os.environ['CCP4I2_ROOT'] = str(PROJECT_ROOT)
os.environ['DJANGO_SETTINGS_MODULE'] = 'ccp4x.settings'

# Add project root to path
sys.path.insert(0, str(PROJECT_ROOT))

from core.CCP4XtalData import CUnmergedDataFileList, CUnmergedDataFile


def test_unmerged_list_subitem():
    """Test that CUnmergedDataFileList has correct subItem qualifier"""

    # Create a list
    unmerged_list = CUnmergedDataFileList(name="test_list")

    # Verify subItem qualifier is set
    sub_item = unmerged_list.get_qualifier('subItem')
    assert sub_item is not None, "subItem qualifier should be set"
    assert isinstance(sub_item, dict), "subItem should be a dict"
    assert 'class' in sub_item, "subItem should have 'class' key"
    assert sub_item['class'] == CUnmergedDataFile, f"subItem class should be CUnmergedDataFile, got {sub_item['class']}"

    print(f"✅ CUnmergedDataFileList.subItem qualifier correctly set to {sub_item['class'].__name__}")

    # Test makeItem() creates correct type
    new_item = unmerged_list.makeItem()
    assert isinstance(new_item, CUnmergedDataFile), f"makeItem() should create CUnmergedDataFile, got {type(new_item).__name__}"

    print(f"✅ CUnmergedDataFileList.makeItem() correctly creates {type(new_item).__name__} instances")

    # Test that items have expected attributes
    assert hasattr(new_item, 'fullPath'), "CUnmergedDataFile should have fullPath attribute"
    assert hasattr(new_item, 'loadFile'), "CUnmergedDataFile should have loadFile method"

    print("✅ Created items have expected CUnmergedDataFile attributes")

    return True


if __name__ == "__main__":
    try:
        test_unmerged_list_subitem()
        print("\n✅ SUCCESS: CUnmergedDataFileList subItem type is correctly configured")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
