#!/usr/bin/env python3
"""
Test that CList.append() performs smart type conversion for strings -> CDataFile.
This fixes the aimless wrapper bug where it appends string paths to MTZUNMERGEDOUT.
"""

import sys

from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.CCP4XtalData import CUnmergedMtzDataFile


def test_clist_smart_append_string_to_file():
    """Test that appending a string to a CDataFile list creates a file object."""

    # Create a list with CUnmergedMtzDataFile as subItem
    file_list = CList(name="MTZUNMERGEDOUT")
    file_list.set_qualifier('subItem', {'class': CUnmergedMtzDataFile, 'qualifiers': {}})

    # Simulate what aimless.py does: append a string path
    test_path = "/path/to/unmerged.mtz"
    file_list.append(test_path)

    # Verify smart conversion happened
    assert len(file_list) == 1, "List should have one item"

    item = file_list[0]
    assert isinstance(item, CUnmergedMtzDataFile), f"Item should be CUnmergedMtzDataFile, got {type(item).__name__}"
    assert hasattr(item, 'fullPath'), "Item should have fullPath attribute"

    # Verify the path was set correctly
    assert item.fullPath == test_path, f"Expected fullPath={test_path}, got {item.fullPath}"

    print(f"✅ Smart conversion: string '{test_path}' -> CUnmergedMtzDataFile with fullPath")
    print(f"✅ Item type: {type(item).__name__}")
    print(f"✅ Item has fullPath: {item.fullPath}")

    return True


def test_clist_append_file_object_directly():
    """Test that appending a file object directly still works."""

    file_list = CList(name="MTZUNMERGEDOUT")
    file_list.set_qualifier('subItem', {'class': CUnmergedMtzDataFile, 'qualifiers': {}})

    # Create file object manually
    file_obj = CUnmergedMtzDataFile()
    file_obj.setFullPath("/direct/path.mtz")

    file_list.append(file_obj)

    assert len(file_list) == 1
    assert file_list[0] is file_obj, "Direct file object append should preserve identity"
    assert file_list[0].fullPath == "/direct/path.mtz"

    print("✅ Direct file object append still works")

    return True


if __name__ == "__main__":
    try:
        test_clist_smart_append_string_to_file()
        test_clist_append_file_object_directly()
        print("\n✅ SUCCESS: CList.append() smart type conversion working correctly")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
