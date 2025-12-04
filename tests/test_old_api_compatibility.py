"""Test old CCP4i2 API compatibility methods."""

import sys
import os

import pytest

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.base_object.base_classes import CData, CContainer
from core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean


def test_object_path():
    """Test objectPath() returns hierarchical path."""
    # Create a hierarchy: container -> subcontainer -> value
    root = CContainer(name="root")
    child = CContainer(parent=root, name="child")
    value = CInt(parent=child, name="value")

    # Test paths
    assert root.objectPath() == "root"
    assert child.objectPath() == "root.child"
    assert value.objectPath() == "root.child.value"


def test_object_name():
    """Test objectName() returns object name."""
    obj = CInt(name="myInt")
    assert obj.objectName() == "myInt"

    obj2 = CString(name="myString")
    assert obj2.objectName() == "myString"


def test_set_default():
    """Test setDefault() sets default value."""
    num = CInt()
    num.setDefault(42)

    # Value should be set
    assert num.value == 42
    # Should be marked as DEFAULT state
    from core.base_object.base_classes import ValueState
    assert num.getValueState('value') == ValueState.DEFAULT


def test_is_set_with_field_name():
    """Test isSet() with explicit field name."""
    obj = CInt(name="test")

    # Not set initially
    assert not obj.isSet('value')

    # Set the value
    obj.value = 10

    # Should be set now
    assert obj.isSet('value')


def test_unset():
    """Test unSet() returns field to not-set state."""
    obj = CInt(name="test")
    obj.value = 50

    # Should be set
    assert obj.isSet('value')

    # Unset it
    obj.unSet('value')

    # Should not be set anymore
    assert not obj.isSet('value')


def test_container_add_content():
    """Test addContent() adds new content item."""
    container = CContainer(name="test")

    # Add new CInt content
    new_int = container.addContent(CInt, "myInt")

    # Should be accessible
    assert hasattr(container, 'myInt')
    assert container.myInt is new_int
    assert new_int.name == "myInt"
    assert new_int.parent is container


def test_container_add_object():
    """Test addObject() adds existing object."""
    container = CContainer(name="test")

    # Create an object
    obj = CFloat(name="myFloat")
    obj.value = 3.14

    # Add it to container
    result = container.addObject(obj)

    # Should be accessible
    assert hasattr(container, 'myFloat')
    assert container.myFloat is obj
    assert result is obj
    assert obj.parent is container


def test_container_delete_object():
    """Test deleteObject() removes object from container."""
    container = CContainer(name="test")

    # Add content
    container.addContent(CInt, "value1")
    container.addContent(CString, "value2")

    # Verify it exists
    assert hasattr(container, 'value1')
    assert hasattr(container, 'value2')

    # Delete value1
    container.deleteObject('value1')

    # Should be gone
    assert not hasattr(container, 'value1')
    # value2 should still exist
    assert hasattr(container, 'value2')


def test_container_data_order():
    """Test dataOrder() returns content order."""
    container = CContainer(name="test")

    # Add content in specific order
    container.addContent(CInt, "third")
    container.addContent(CString, "first")
    container.addContent(CFloat, "second")

    # Get order
    order = container.dataOrder()

    # Should match insertion order
    assert order == ["third", "first", "second"]


def test_container_clear():
    """Test clear() removes all content."""
    container = CContainer(name="test")

    # Add multiple items
    container.addContent(CInt, "item1")
    container.addContent(CString, "item2")
    container.addContent(CBoolean, "item3")

    # Verify they exist
    assert hasattr(container, 'item1')
    assert hasattr(container, 'item2')
    assert hasattr(container, 'item3')

    # Clear the container
    container.clear()

    # All should be gone
    assert not hasattr(container, 'item1')
    assert not hasattr(container, 'item2')
    assert not hasattr(container, 'item3')

    # Data order should be empty
    assert container.dataOrder() == []


def test_container_add_content_by_string():
    """Test addContent() with string class name."""
    container = CContainer(name="test")

    # Add using string class name
    obj = container.addContent('CInt', 'myValue')

    # Should create correct type
    assert isinstance(obj, CInt)
    assert hasattr(container, 'myValue')
    assert container.myValue is obj


def test_full_hierarchy_workflow():
    """Test complete workflow using old API methods."""
    # Create task structure
    task = CContainer(name="servalcat_pipe")

    # Add input data container
    input_data = task.addContent(CContainer, "inputData")

    # Add control parameters container
    ctrl = task.addContent(CContainer, "controlParameters")

    # Add parameters to control
    ncycles = ctrl.addContent(CInt, "NCYCLES")
    ncycles.setDefault(10)

    add_waters = ctrl.addContent(CBoolean, "ADD_WATERS")
    add_waters.setDefault(False)

    # Test paths
    assert task.objectPath() == "servalcat_pipe"
    assert ctrl.objectPath() == "servalcat_pipe.controlParameters"
    assert ncycles.objectPath() == "servalcat_pipe.controlParameters.NCYCLES"

    # Test values
    assert ncycles.value == 10
    assert add_waters.value is False

    # Test data order
    assert ctrl.dataOrder() == ["NCYCLES", "ADD_WATERS"]

    # Modify and test isSet
    assert not ncycles.isSet('value')  # Default value
    ncycles.value = 25
    assert ncycles.isSet('value')  # Explicitly set

    # Delete an object
    ctrl.deleteObject("ADD_WATERS")
    assert not hasattr(ctrl, "ADD_WATERS")
    assert ctrl.dataOrder() == ["NCYCLES"]
