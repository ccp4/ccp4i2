import pytest

# Try to import CPdbDataFile from the generated core module
from ccp4i2.core.CCP4ModelData import CPdbDataFile

def test_cpdbdatafile_instantiation():
    # Attempt to instantiate CPdbDataFile
    obj = CPdbDataFile()
    assert obj is not None
    # Optionally, check for expected attributes or behaviors
    # Example: assert hasattr(obj, 'some_expected_attribute')

def test_cpdbdatafile_set_from_dict():
    # Attempt to instantiate CPdbDataFile
    obj = CPdbDataFile()
    obj.set({"baseName": "test_file", "dbFileId": "12345"})

    # baseName and dbFileId are CData wrappers (CFilePath, CUUID), not plain strings
    # Access their values via .value attribute
    assert hasattr(obj.baseName, 'value')
    assert obj.baseName.value == "test_file"
    assert hasattr(obj.dbFileId, 'value')
    assert obj.dbFileId.value == "12345"

    # relPath should not be set (or None/empty)
    assert not hasattr(obj.relPath, 'value') or obj.relPath.value is None or obj.relPath.value == ""

    assert isinstance(obj, CPdbDataFile)
