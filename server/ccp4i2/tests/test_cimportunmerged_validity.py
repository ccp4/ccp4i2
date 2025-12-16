"""Test CImportUnmerged validity with content_qualifiers.

These tests verify that when a new CImportUnmerged item is added to a list
via addItem(), the per-field qualifiers (like allowUndefined=False) from
content_qualifiers are properly applied and generate validity errors.
"""

import pytest
from ccp4i2.core.cdata_stubs.CCP4XtalData import CImportUnmergedStub


class TestCImportUnmergedValidity:
    """Test validity checking for CImportUnmerged items."""

    def test_standalone_cimportunmerged_has_content_qualifiers(self):
        """Test that CImportUnmergedStub has content_qualifiers in metadata."""
        item = CImportUnmergedStub(name='test_item')

        # Check metadata has content_qualifiers
        assert hasattr(item, '_metadata')
        assert item._metadata is not None
        assert item._metadata.content_qualifiers is not None

        # Check expected fields have allowUndefined=False
        cq = item._metadata.content_qualifiers
        assert cq.get('file', {}).get('allowUndefined') is False
        assert cq.get('crystalName', {}).get('allowUndefined') is False
        assert cq.get('dataset', {}).get('allowUndefined') is False

    def test_content_qualifiers_applied_to_child_fields(self):
        """Test that content_qualifiers are applied to child field instances."""
        item = CImportUnmergedStub(name='test_item')

        # file should have allowUndefined=False from content_qualifiers
        assert item.file is not None
        assert item.file.get_qualifier('allowUndefined') is False

        # crystalName should have allowUndefined=False
        assert item.crystalName is not None
        assert item.crystalName.get_qualifier('allowUndefined') is False

        # dataset should have allowUndefined=False
        assert item.dataset is not None
        assert item.dataset.get_qualifier('allowUndefined') is False

    def test_validity_reports_required_fields_not_set(self):
        """Test that validity() reports errors for required but unset fields."""
        item = CImportUnmergedStub(name='test_item')

        report = item.validity()

        # Should have errors - minimum 3 for file, crystalName, dataset
        # Plus 3 for cell.a, cell.b, cell.c (from CCellLength's class qualifiers)
        assert len(report) >= 6

        # Extract error paths
        error_paths = [str(err) for err in report]

        # Check required fields are reported
        assert any('file' in path for path in error_paths), f"Missing 'file' error: {error_paths}"
        assert any('crystalName' in path for path in error_paths), f"Missing 'crystalName' error: {error_paths}"
        assert any('dataset' in path for path in error_paths), f"Missing 'dataset' error: {error_paths}"


class TestCImportUnmergedViaAimlessPipe:
    """Test CImportUnmerged validity when created via aimless_pipe plugin."""

    def test_additem_creates_item_with_content_qualifiers(self):
        """Test that addItem() creates items with content_qualifiers applied."""
        from ccp4i2.core.task_manager.plugin_registry import PluginRegistry

        registry = PluginRegistry()
        plugin_class = registry.get_plugin_class('aimless_pipe')
        plugin = plugin_class()

        # Get the UNMERGEDFILES list
        unmerged_list = plugin.container.inputData.UNMERGEDFILES

        # Add a new item via addItem()
        item = unmerged_list.addItem()

        # Check item class and metadata
        assert item.__class__.__name__ == 'CImportUnmerged'
        assert hasattr(item, '_metadata')
        assert item._metadata.content_qualifiers is not None

        # Check file has allowUndefined=False applied
        assert item.file is not None
        assert item.file.get_qualifier('allowUndefined') is False

    def test_list_validity_includes_item_errors(self):
        """Test that list validity() aggregates errors from items."""
        from ccp4i2.core.task_manager.plugin_registry import PluginRegistry

        registry = PluginRegistry()
        plugin_class = registry.get_plugin_class('aimless_pipe')
        plugin = plugin_class()

        # Get the UNMERGEDFILES list
        unmerged_list = plugin.container.inputData.UNMERGEDFILES

        # Add a new empty item
        item = unmerged_list.addItem()

        # Get validity report from the list
        list_report = unmerged_list.validity()

        # Should include errors from the item (file, crystalName, dataset, cell.*)
        assert len(list_report) >= 6

        # Check that item errors are in the list report
        error_paths = [str(err) for err in list_report]
        assert any('file' in path for path in error_paths)
        assert any('crystalName' in path for path in error_paths)
        assert any('dataset' in path for path in error_paths)

    def test_container_validity_includes_nested_errors(self):
        """Test that container validity includes nested list item errors."""
        from ccp4i2.core.task_manager.plugin_registry import PluginRegistry

        registry = PluginRegistry()
        plugin_class = registry.get_plugin_class('aimless_pipe')
        plugin = plugin_class()

        # Add an empty item to UNMERGEDFILES
        plugin.container.inputData.UNMERGEDFILES.addItem()

        # Get validity from container
        container_report = plugin.container.validity()

        # Should include errors from nested items
        error_paths = [str(err) for err in container_report]

        # Check for UNMERGEDFILES item errors
        assert any('UNMERGEDFILES' in path and 'file' in path for path in error_paths), \
            f"Missing UNMERGEDFILES.*.file error: {error_paths}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
