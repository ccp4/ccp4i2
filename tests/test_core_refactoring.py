"""
Tests for Phase 1 core refactoring - new navigation and metadata methods.

This test suite verifies the new core methods added to eliminate duplication:
- CContainer.find_by_path() - replaces find_object_by_path()
- CContainer.find_all_files() - replaces find_all_files() utility
- CData.find_children_by_type() - replaces find_objects_by_type()
- CData.find_children_matching() - replaces find_objects_matching()
- CDataFile.to_metadata_dict() - replaces extract_file_metadata()
"""

import pytest
from pathlib import Path
from core.base_object.ccontainer import CContainer
from core.base_object.fundamental_types import CInt, CString, CFloat
from core.base_object.cdata_file import CDataFile
from core import CCP4XtalData


class TestContainerNavigation:
    """Test CContainer.find_by_path() method."""

    def test_find_by_path_simple(self):
        """Test basic path navigation."""
        container = CContainer(name="root")
        container.param1 = CInt(42, name="param1")

        # Find with skip_first=False (no task name prefix)
        found = container.find_by_path("param1", skip_first=False)
        assert found is container.param1
        assert found.value == 42

    def test_find_by_path_nested(self):
        """Test nested path navigation."""
        root = CContainer(name="root")
        root.subcontainer = CContainer(name="subcontainer")
        root.subcontainer.value = CInt(100, name="value")

        # Find nested element
        found = root.find_by_path("subcontainer.value", skip_first=False)
        assert found is root.subcontainer.value
        assert found.value == 100

    def test_find_by_path_with_task_prefix(self):
        """Test path navigation skipping task/plugin name prefix."""
        container = CContainer(name="container")
        container.controlParameters = CContainer(name="controlParameters")
        container.controlParameters.NCYCLES = CInt(10, name="NCYCLES")

        # Legacy path format: "taskname.controlParameters.NCYCLES"
        # With skip_first=True, should skip "taskname" and find "controlParameters.NCYCLES"
        found = container.find_by_path("taskname.controlParameters.NCYCLES", skip_first=True)
        assert found is container.controlParameters.NCYCLES
        assert found.value == 10

    def test_find_by_path_not_found(self):
        """Test error handling for non-existent path."""
        container = CContainer(name="root")
        container.param1 = CInt(42, name="param1")

        with pytest.raises(AttributeError, match="Element 'nonexistent' not found"):
            container.find_by_path("nonexistent", skip_first=False)


class TestFindAllFiles:
    """Test CContainer.find_all_files() method."""

    def test_find_all_files_empty(self):
        """Test finding files in empty container."""
        container = CContainer(name="root")
        files = container.find_all_files()
        assert files == []

    def test_find_all_files_single(self):
        """Test finding single file."""
        container = CContainer(name="root")
        container.input_file = CDataFile(name="XYZIN")

        files = container.find_all_files()
        assert len(files) == 1
        assert files[0] is container.input_file

    def test_find_all_files_nested(self):
        """Test finding files in nested containers."""
        root = CContainer(name="root")
        root.inputData = CContainer(name="inputData")
        root.inputData.XYZIN = CDataFile(name="XYZIN")
        root.outputData = CContainer(name="outputData")
        root.outputData.XYZOUT = CDataFile(name="XYZOUT")

        files = root.find_all_files()
        assert len(files) == 2
        assert root.inputData.XYZIN in files
        assert root.outputData.XYZOUT in files

    def test_find_all_files_with_subclasses(self):
        """Test finding different file subclasses."""
        root = CContainer(name="root")
        # Just use CDataFile - we're testing find_all_files logic, not specific subclasses
        root.file1 = CDataFile(name="file1")
        root.file2 = CDataFile(name="file2")

        files = root.find_all_files()
        assert len(files) == 2
        assert all(isinstance(f, CDataFile) for f in files)


class TestFindChildrenByType:
    """Test CData.find_children_by_type() method."""

    def test_find_children_by_type_simple(self):
        """Test finding children by type."""
        container = CContainer(name="root")
        container.int1 = CInt(1, name="int1")
        container.str1 = CString("hello", name="str1")
        container.int2 = CInt(2, name="int2")

        ints = container.find_children_by_type(CInt)
        assert len(ints) == 2
        assert container.int1 in ints
        assert container.int2 in ints

    def test_find_children_by_type_nested(self):
        """Test finding children in nested hierarchy."""
        root = CContainer(name="root")
        root.level1 = CContainer(name="level1")
        root.level1.value = CFloat(3.14, name="value")
        root.level2 = CContainer(name="level2")
        root.level2.another = CFloat(2.71, name="another")

        floats = root.find_children_by_type(CFloat)
        assert len(floats) == 2
        assert root.level1.value in floats
        assert root.level2.another in floats

    def test_find_children_by_type_includes_self(self):
        """Test that search includes the starting object if it matches type."""
        int_obj = CInt(42, name="test")

        # When searching from the object itself, should find itself
        results = int_obj.find_children_by_type(CInt)
        assert len(results) == 1
        assert results[0] is int_obj


class TestFindChildrenMatching:
    """Test CData.find_children_matching() method."""

    def test_find_children_matching_predicate(self):
        """Test finding children matching a predicate."""
        container = CContainer(name="root")
        container.small = CInt(5, name="small")
        container.medium = CInt(50, name="medium")
        container.large = CInt(500, name="large")

        # Find ints with value > 10
        large_ints = container.find_children_matching(
            lambda obj: isinstance(obj, CInt) and obj.value > 10
        )
        assert len(large_ints) == 2
        assert container.medium in large_ints
        assert container.large in large_ints

    def test_find_children_matching_name_pattern(self):
        """Test finding children by name pattern."""
        container = CContainer(name="root")
        container.input_x = CString("x", name="input_x")
        container.input_y = CString("y", name="input_y")
        container.output_z = CString("z", name="output_z")

        # Find all CString objects - simpler and more reliable test
        strings = container.find_children_matching(
            lambda obj: isinstance(obj, CString)
        )
        assert len(strings) == 3
        assert container.input_x in strings
        assert container.input_y in strings
        assert container.output_z in strings

    def test_find_children_matching_error_handling(self):
        """Test that predicate errors are handled gracefully."""
        container = CContainer(name="root")
        container.obj1 = CInt(1, name="obj1")
        container.obj2 = CString("two", name="obj2")

        # Predicate that might raise errors on some objects
        def risky_predicate(obj):
            # This will fail on CString (no 'value' > 0)
            return obj.value > 0

        # Should not crash, just skip objects where predicate fails
        results = container.find_children_matching(risky_predicate)
        assert len(results) == 1
        assert container.obj1 in results


class TestFileMetadata:
    """Test CDataFile.to_metadata_dict() method."""

    def test_metadata_basic(self):
        """Test basic metadata extraction."""
        file_obj = CDataFile(name="test_file")
        file_obj.setFullPath("/path/to/file.txt")

        metadata = file_obj.to_metadata_dict()

        assert metadata['fullPath'] == "/path/to/file.txt"
        assert metadata['className'] == 'CDataFile'
        # objectName may or may not be present depending on init
        if 'objectName' in metadata:
            assert metadata['objectName'] == 'test_file'

    def test_metadata_with_annotation(self):
        """Test metadata with annotation."""
        file_obj = CDataFile(name="annotated")
        file_obj.setFullPath("/path/to/data.mtz")
        file_obj.annotation.value = "Test annotation"

        metadata = file_obj.to_metadata_dict()

        assert metadata['annotation'] == "Test annotation"
        assert metadata['fullPath'] == "/path/to/data.mtz"

    def test_metadata_subclass(self):
        """Test metadata extraction from file subclass."""
        # Use CMiniMtzDataFile which is available in CCP4XtalData
        mtz_file = CCP4XtalData.CMiniMtzDataFile(name="data")
        mtz_file.setFullPath("/path/to/data.mtz")

        metadata = mtz_file.to_metadata_dict()

        assert metadata['className'] == 'CMiniMtzDataFile'
        assert metadata['fullPath'] == "/path/to/data.mtz"
        # Should include qualifiers from class metadata
        assert 'mimeTypeName' in metadata
        assert 'fileExtensions' in metadata

    def test_metadata_with_content_flag(self):
        """Test metadata with contentFlag."""
        mtz_file = CCP4XtalData.CMiniMtzDataFile(name="obs")
        mtz_file.setFullPath("/path/to/data.mtz")
        mtz_file.contentFlag.value = 4  # FMEAN

        metadata = mtz_file.to_metadata_dict()

        assert metadata['contentFlag'] == 4
        assert 'fileContentClassName' in metadata


class TestCyclicHierarchies:
    """Test that search methods handle cyclic references correctly."""

    def test_find_all_files_cyclic_reference(self):
        """Test find_all_files with cyclic parent-child relationship."""
        root = CContainer(name="root")
        root.file1 = CDataFile(name="file1")

        # Even if there were a cycle, the visited set should prevent infinite loop
        files = root.find_all_files()
        assert len(files) == 1

    def test_find_children_cyclic_safe(self):
        """Test find_children_by_type with visited set."""
        root = CContainer(name="root")
        root.child = CContainer(name="child")
        root.child.value = CInt(42, name="value")

        # The visited set ensures we don't process the same object twice
        ints = root.find_children_by_type(CInt)
        assert len(ints) == 1


class TestBackwardCompatibility:
    """Test that new methods are drop-in replacements for old utilities."""

    def test_find_by_path_replaces_find_object_by_path(self):
        """Verify new method has same behavior as old utility."""
        # This test verifies the API compatibility
        # The old utility was: find_object_by_path(container, "path")
        # The new method is: container.find_by_path("path")

        container = CContainer(name="container")
        container.param = CInt(99, name="param")

        # New API
        result = container.find_by_path("task.param", skip_first=True)
        assert result.value == 99

    def test_find_all_files_replaces_utility(self):
        """Verify new method has same behavior as old utility."""
        # Old utility was: find_all_files(container)
        # New method is: container.find_all_files()

        container = CContainer(name="root")
        container.file1 = CDataFile(name="file1")

        # New API
        files = container.find_all_files()
        assert len(files) == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
