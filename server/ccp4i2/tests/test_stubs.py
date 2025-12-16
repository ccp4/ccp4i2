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
        # dbFileId should be unset (removed or None)
        assert not hasattr(f, "dbFileId") or f.dbFileId is None

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