"""
Tests to verify that deprecation warnings are emitted correctly.

This validates Phase 2 of the refactoring: ensuring all deprecated utility
functions emit proper warnings when called.
"""

import pytest
import warnings
from ccp4i2.core.base_object.ccontainer import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CString
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.core import CCP4XtalData


class TestDeprecationWarnings:
    """Test that deprecated utility functions emit warnings."""

    def test_find_object_by_path_warning(self):
        """Test find_object_by_path() emits deprecation warning."""
        from server.ccp4i2.lib.utils.containers.find_objects import find_object_by_path

        container = CContainer(name="root")
        container.param = CInt(42, name="param")

        # Check that deprecation warning is emitted
        with pytest.warns(DeprecationWarning, match="find_object_by_path.*deprecated"):
            result = find_object_by_path(container, "task.param")
            assert result.value == 42

    def test_find_all_files_warning(self):
        """Test find_all_files() emits deprecation warning."""
        from server.ccp4i2.lib.cdata_utils import find_all_files

        container = CContainer(name="root")
        container.file1 = CDataFile(name="file1")

        # Check that deprecation warning is emitted
        with pytest.warns(DeprecationWarning, match="find_all_files.*deprecated"):
            files = find_all_files(container)
            assert len(files) == 1

    def test_find_objects_by_type_warning(self):
        """Test find_objects_by_type() emits deprecation warning."""
        from server.ccp4i2.lib.cdata_utils import find_objects_by_type

        container = CContainer(name="root")
        container.int1 = CInt(1, name="int1")
        container.str1 = CString("hello", name="str1")

        # Check that deprecation warning is emitted
        with pytest.warns(DeprecationWarning, match="find_objects_by_type.*deprecated"):
            ints = find_objects_by_type(container, CInt)
            assert len(ints) == 1

    def test_find_objects_matching_warning(self):
        """Test find_objects_matching() emits deprecation warning."""
        from server.ccp4i2.lib.cdata_utils import find_objects_matching

        container = CContainer(name="root")
        container.int1 = CInt(5, name="int1")
        container.int2 = CInt(50, name="int2")

        # Check that deprecation warning is emitted
        with pytest.warns(DeprecationWarning, match="find_objects_matching.*deprecated"):
            large_ints = find_objects_matching(
                container,
                lambda obj: isinstance(obj, CInt) and obj.value > 10
            )
            assert len(large_ints) == 1

    def test_extract_file_metadata_warning(self):
        """Test extract_file_metadata() emits deprecation warning."""
        from server.ccp4i2.lib.cdata_utils import extract_file_metadata

        mtz_file = CCP4XtalData.CMiniMtzDataFile(name="data")
        mtz_file.setFullPath("/path/to/data.mtz")

        # Check that deprecation warning is emitted
        with pytest.warns(DeprecationWarning, match="extract_file_metadata.*deprecated"):
            metadata = extract_file_metadata(mtz_file)
            assert 'name' in metadata
            assert 'file_type' in metadata


class TestDeprecatedFunctionsStillWork:
    """Verify that deprecated functions still work correctly (backward compatibility)."""

    def test_find_object_by_path_still_works(self):
        """Verify find_object_by_path() still works despite deprecation."""
        from server.ccp4i2.lib.utils.containers.find_objects import find_object_by_path

        container = CContainer(name="root")
        container.subcontainer = CContainer(name="subcontainer")
        container.subcontainer.value = CInt(99, name="value")

        # Suppress the warning for this test
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            result = find_object_by_path(container, "task.subcontainer.value")
            assert result.value == 99

    def test_find_all_files_still_works(self):
        """Verify find_all_files() still works despite deprecation."""
        from server.ccp4i2.lib.cdata_utils import find_all_files

        container = CContainer(name="root")
        container.inputData = CContainer(name="inputData")
        container.inputData.XYZIN = CDataFile(name="XYZIN")
        container.outputData = CContainer(name="outputData")
        container.outputData.XYZOUT = CDataFile(name="XYZOUT")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            files = find_all_files(container)
            assert len(files) == 2

    def test_find_objects_by_type_still_works(self):
        """Verify find_objects_by_type() still works despite deprecation."""
        from server.ccp4i2.lib.cdata_utils import find_objects_by_type

        container = CContainer(name="root")
        container.int1 = CInt(1, name="int1")
        container.int2 = CInt(2, name="int2")
        container.str1 = CString("test", name="str1")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            ints = find_objects_by_type(container, CInt)
            assert len(ints) == 2

    def test_find_objects_matching_still_works(self):
        """Verify find_objects_matching() still works despite deprecation."""
        from server.ccp4i2.lib.cdata_utils import find_objects_matching

        container = CContainer(name="root")
        container.small = CInt(5, name="small")
        container.large = CInt(500, name="large")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            large_ints = find_objects_matching(
                container,
                lambda obj: isinstance(obj, CInt) and obj.value > 10
            )
            assert len(large_ints) == 1
            assert large_ints[0].value == 500

    def test_extract_file_metadata_still_works(self):
        """Verify extract_file_metadata() still works despite deprecation."""
        from server.ccp4i2.lib.cdata_utils import extract_file_metadata

        file_obj = CDataFile(name="test_file")
        file_obj.setFullPath("/path/to/file.txt")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            metadata = extract_file_metadata(file_obj)
            assert metadata['name'] == 'test_file'
            assert 'file_type' in metadata
            assert 'is_set' in metadata


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
