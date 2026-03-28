#!/usr/bin/env python3
# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
