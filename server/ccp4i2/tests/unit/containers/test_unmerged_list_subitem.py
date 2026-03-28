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
"""Test that CUnmergedDataFileList properly specifies its subItem type."""

from ccp4i2.core.CCP4XtalData import CUnmergedDataFileList, CUnmergedDataFile


def test_unmerged_list_subitem():
    """Test that CUnmergedDataFileList has correct subItem qualifier."""
    unmerged_list = CUnmergedDataFileList(name="test_list")

    sub_item = unmerged_list.get_qualifier('subItem')
    assert sub_item is not None, "subItem qualifier should be set"
    assert isinstance(sub_item, dict), "subItem should be a dict"
    assert sub_item['class'] == CUnmergedDataFile

    new_item = unmerged_list.makeItem()
    assert isinstance(new_item, CUnmergedDataFile)
    assert hasattr(new_item, 'fullPath')
    assert hasattr(new_item, 'loadFile')
