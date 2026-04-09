"""
Test smart list handling in upload_param API.

Tests that uploading files to list parameters works correctly:
1. Uploading without index appends a new item
2. Uploading with existing index uses that item
3. Uploading with index beyond list length creates items up to that index
"""

from unittest.mock import Mock
from ccp4i2.core.base_object.fundamental_types import CList
from ccp4i2.core.base_object.cdata_file import CDataFile
from ccp4i2.lib.utils.files.upload_param import resolve_list_upload_path


class TestSmartListResolution:
    """Test the resolve_list_upload_path helper function."""

    def test_append_to_list_without_index(self):
        """Test that uploading to a list without index appends a new item."""
        # Create a mock container with a list parameter
        mock_container = Mock()
        mock_list = Mock(spec=CList)
        mock_list.__len__ = Mock(return_value=2)  # List has 2 items

        # Create a new item that will be appended
        mock_new_item = Mock(spec=CDataFile)

        # Mock find_by_path to return the list
        def mock_find_by_path(path, skip_first=True):
            if path == "task.prosmartProtein.REFERENCE_MODELS":
                return mock_list
            elif path == "task.prosmartProtein.REFERENCE_MODELS[2]":
                return mock_new_item
            raise AttributeError(f"Path not found: {path}")

        mock_container.find_by_path = Mock(side_effect=mock_find_by_path)

        # Mock list.addItem() to increment length and return the new item
        def mock_add_item():
            mock_list.__len__ = Mock(return_value=3)
            return mock_new_item
        mock_list.addItem = Mock(side_effect=mock_add_item)

        # Call resolve_list_upload_path
        param_object, final_path = resolve_list_upload_path(
            mock_container,
            "task.prosmartProtein.REFERENCE_MODELS",
            skip_first=True
        )

        # Verify addItem was called
        mock_list.addItem.assert_called_once()

        # Verify the final path has [2] appended (new item at index 2)
        assert final_path == "task.prosmartProtein.REFERENCE_MODELS[2]"
        assert param_object == mock_new_item

    def test_use_existing_list_item(self):
        """Test that uploading with an existing index uses that item."""
        mock_container = Mock()
        mock_list = Mock(spec=CList)
        mock_list.__len__ = Mock(return_value=3)  # List has 3 items

        mock_existing_item = Mock(spec=CDataFile)

        def mock_find_by_path(path, skip_first=True):
            if path.endswith("REFERENCE_MODELS"):
                return mock_list
            elif path.endswith("REFERENCE_MODELS[1]"):
                return mock_existing_item
            raise AttributeError(f"Path not found: {path}")

        mock_container.find_by_path = Mock(side_effect=mock_find_by_path)
        mock_list.addItem = Mock()  # Should NOT be called

        # Call with existing index [1]
        param_object, final_path = resolve_list_upload_path(
            mock_container,
            "task.prosmartProtein.REFERENCE_MODELS[1]",
            skip_first=True
        )

        # Verify addItem was NOT called (item already exists)
        mock_list.addItem.assert_not_called()

        # Verify path unchanged
        assert final_path == "task.prosmartProtein.REFERENCE_MODELS[1]"
        assert param_object == mock_existing_item

    def test_expand_list_to_reach_index(self):
        """Test that uploading with index beyond length creates items."""
        mock_container = Mock()
        mock_list = Mock(spec=CList)
        current_length = [2]  # List has 2 items initially

        # Mock __len__ to reflect dynamic length
        mock_list.__len__ = Mock(side_effect=lambda: current_length[0])

        mock_target_item = Mock(spec=CDataFile)

        def mock_find_by_path(path, skip_first=True):
            if path.endswith("REFERENCE_MODELS"):
                return mock_list
            elif path.endswith("REFERENCE_MODELS[5]"):
                return mock_target_item
            raise AttributeError(f"Path not found: {path}")

        mock_container.find_by_path = Mock(side_effect=mock_find_by_path)

        # Mock addItem to increment length and return a mock item
        def mock_add_item():
            current_length[0] += 1
            return Mock(spec=CDataFile)
        mock_list.addItem = Mock(side_effect=mock_add_item)

        # Call with index [5] when list only has 2 items
        param_object, final_path = resolve_list_upload_path(
            mock_container,
            "task.prosmartProtein.REFERENCE_MODELS[5]",
            skip_first=True
        )

        # Verify addItem was called 4 times (to get from 2 items to 6 items)
        assert mock_list.addItem.call_count == 4

        # Verify path unchanged
        assert final_path == "task.prosmartProtein.REFERENCE_MODELS[5]"
        assert param_object == mock_target_item

    def test_non_list_parameter_unchanged(self):
        """Test that non-list parameters pass through unchanged."""
        mock_container = Mock()
        mock_file = Mock(spec=CDataFile)

        mock_container.find_by_path = Mock(return_value=mock_file)

        # Call with a non-list parameter
        param_object, final_path = resolve_list_upload_path(
            mock_container,
            "task.inputData.XYZIN",
            skip_first=True
        )

        # Verify path unchanged
        assert final_path == "task.inputData.XYZIN"
        assert param_object == mock_file
