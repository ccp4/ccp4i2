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
from ccp4i2.core.base_object.base_classes import CDataFile

class TestCDataFile:
    @classmethod
    def setup_class(cls):
        # Class-level setup (runs once before all tests)
        cls.shared_resource = "initialized"

    def test_set_unsets_missing_attributes(self):
        f = CDataFile()
        f.baseName = "original"
        f.dbFileId = "Banana"
        assert f.baseName == "original"
        assert f.dbFileId == "Banana"
        f.set({"baseName": "changed"})
        assert f.baseName == "changed"
        # dbFileId should be unset (reset to empty/None/default)
        assert not f.dbFileId  # empty string, None, or falsy

    def test_update_only_changes_specified(self):
        f = CDataFile()
        f.baseName = "original"
        f.dbFileId = "Banana"
        f.update({"baseName": "changed"})
        assert f.baseName == "changed"
        assert f.dbFileId == "Banana"

    def test_set_and_update_methods_exist(self):
        f = CDataFile()
        assert hasattr(f, "set")
        assert hasattr(f, "update")