"""Test that CContainer allows attribute-style access to children."""

import pytest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ccp4i2.core.base_object.base_classes import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CString


def test_container_attribute_access_to_children():
    """Test that children can be accessed via attributes."""
    # Create a container
    container = CContainer(name="testContainer")

    # Add a child container named "inputData"
    input_data = CContainer(name="inputData")
    input_data.set_parent(container)

    # Test: Access via attribute should work
    assert container.inputData is input_data
    assert container.inputData.name == "inputData"


def test_container_nested_child_access():
    """Test that nested children can be accessed via chained attributes."""
    # Create hierarchy
    container = CContainer(name="task")
    input_data = CContainer(name="inputData")
    input_data.set_parent(container)

    file_param = CString(value="test.pdb", name="XYZIN")
    file_param.set_parent(input_data)

    # Test: Chained access should work
    assert container.inputData.XYZIN is file_param
    assert container.inputData.XYZIN.value == "test.pdb"


def test_container_nonexistent_child_raises_error():
    """Test that accessing non-existent children raises AttributeError."""
    container = CContainer(name="testContainer")

    # Test: Non-existent child should raise AttributeError
    with pytest.raises(AttributeError, match="no attribute 'nonExistent'"):
        _ = container.nonExistent


def test_container_child_access_with_addContent():
    """Test that children added via addContent are accessible."""
    container = CContainer(name="task")

    # Add using addContent (old API)
    input_data = container.addContent(CContainer, "inputData")
    ncycles = input_data.addContent(CInt, "NCYCLES")
    ncycles.value = 10

    # Test: Should be accessible via attributes
    assert container.inputData is input_data
    assert container.inputData.NCYCLES is ncycles
    assert container.inputData.NCYCLES.value == 10


def test_container_deleted_child_not_accessible():
    """Test that deleted children are no longer accessible."""
    container = CContainer(name="task")

    # Add child
    child = container.addContent(CContainer, "tempChild")
    assert container.tempChild is child

    # Delete child
    container.deleteObject("tempChild")

    # Test: Should raise AttributeError after deletion
    with pytest.raises(AttributeError, match="no attribute 'tempChild'"):
        _ = container.tempChild


def test_container_child_access_doesnt_interfere_with_methods():
    """Test that __getattr__ doesn't interfere with normal methods."""
    container = CContainer(name="task")
    child = CContainer(name="inputData")
    child.set_parent(container)

    # Test: Normal methods should still work
    assert container.children() == [child]
    assert container.dataOrder() == []
    assert len(container) == 0  # _container_items is empty
