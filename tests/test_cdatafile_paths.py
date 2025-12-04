"""Tests for CDataFile path handling behavior."""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.base_object.base_classes import CDataFile


class TestCDataFilePaths:
    """Tests for CDataFile path handling (non-database mode)."""

    def test_cdatafile_init_with_path(self):
        """Test that CDataFile can be initialized with a file path."""
        file_obj = CDataFile(file_path="/path/to/file.pdb")

        # Should set baseName via setFullPath
        assert file_obj.getFullPath() == "/path/to/file.pdb"
        assert str(file_obj) == "/path/to/file.pdb"

    def test_cdatafile_setfullpath(self):
        """Test setFullPath method."""
        file_obj = CDataFile()

        file_obj.setFullPath("/data/structure.cif")

        assert file_obj.getFullPath() == "/data/structure.cif"
        assert str(file_obj) == "/data/structure.cif"

    def test_cdatafile_getfullpath(self):
        """Test getFullPath method returns string."""
        file_obj = CDataFile(file_path="/my/file.mtz")

        result = file_obj.getFullPath()

        assert isinstance(result, str)
        assert result == "/my/file.mtz"

    def test_cdatafile_fullpath_property(self):
        """Test fullPath property."""
        file_obj = CDataFile()
        file_obj.setFullPath("/test/data.map")

        # Access via property
        path = file_obj.fullPath

        assert isinstance(path, str)
        assert path == "/test/data.map"

    def test_cdatafile_str_method(self):
        """Test __str__ returns the file path."""
        file_obj = CDataFile(file_path="/some/file.pdb")

        result = str(file_obj)

        assert result == "/some/file.pdb"

    def test_cdatafile_empty_path(self):
        """Test that CDataFile with no path returns empty string."""
        file_obj = CDataFile()

        assert file_obj.getFullPath() == ""
        assert file_obj.fullPath == ""
        assert str(file_obj) == ""

    def test_cdatafile_update_path(self):
        """Test that path can be updated."""
        file_obj = CDataFile(file_path="/old/path.pdb")

        assert file_obj.getFullPath() == "/old/path.pdb"

        # Update path
        file_obj.setFullPath("/new/path.pdb")

        assert file_obj.getFullPath() == "/new/path.pdb"
        assert str(file_obj) == "/new/path.pdb"

    def test_cdatafile_with_basename_attribute(self):
        """Test CDataFile when baseName attribute exists."""
        from core.CCP4ModelData import CPdbDataFile

        # CPdbDataFile should have baseName as CFilePath
        pdb_file = CPdbDataFile()

        # Set path
        pdb_file.setFullPath("/structures/model.pdb")

        # Should set baseName.value
        assert pdb_file.getFullPath() == "/structures/model.pdb"
        assert str(pdb_file) == "/structures/model.pdb"

    def test_cdatafile_string_conversion_in_expressions(self):
        """Test that CDataFile can be used in string expressions."""
        file_obj = CDataFile(file_path="/data/input.mtz")

        # Should work in f-strings (calls __str__)
        assert f"File: {file_obj}" == "File: /data/input.mtz"

        # Can concatenate by converting to string
        assert str(file_obj) + ".backup" == "/data/input.mtz.backup"

    def test_cdatafile_comparison_with_string(self):
        """Test that CDataFile can be compared to strings."""
        file_obj = CDataFile(file_path="/path/to/file.pdb")

        # String comparison should work
        assert str(file_obj) == "/path/to/file.pdb"
        assert file_obj.fullPath == "/path/to/file.pdb"

    def test_cdatafile_setfullpath_multiple_times(self):
        """Test that setFullPath can be called multiple times."""
        file_obj = CDataFile()

        file_obj.setFullPath("/first.pdb")
        assert file_obj.getFullPath() == "/first.pdb"

        file_obj.setFullPath("/second.mtz")
        assert file_obj.getFullPath() == "/second.mtz"

        file_obj.setFullPath("/third.cif")
        assert file_obj.getFullPath() == "/third.cif"

    def test_cdatafile_with_parent_and_name(self):
        """Test that CDataFile works with parent and name in hierarchy."""
        from core.base_object.base_classes import CContainer

        container = CContainer(name="files")
        file_obj = CDataFile(
            file_path="/data/model.pdb",
            parent=container,
            name="XYZIN"
        )

        assert file_obj.name == "XYZIN"
        assert file_obj.getFullPath() == "/data/model.pdb"
        assert str(file_obj) == "/data/model.pdb"
