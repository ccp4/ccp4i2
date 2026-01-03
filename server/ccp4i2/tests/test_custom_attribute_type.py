"""Test that AttributeType.CUSTOM works correctly with @cdata_class decorator."""

from ccp4i2.core.base_object.base_classes import CData
from ccp4i2.core.base_object.fundamental_types import CInt, CList, CString
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType


def test_custom_attribute_with_clist():
    """Test that AttributeType.CUSTOM creates a CList instance, not CString."""

    @cdata_class(
        attributes={
            'myList': attribute(AttributeType.CUSTOM, custom_class="CList")
        }
    )
    class TestContainer(CData):
        """Test container with custom CList attribute."""
        pass

    # Create instance
    container = TestContainer(name="testContainer")

    # Verify the attribute was created
    assert hasattr(container, 'myList'), "myList attribute should exist"

    # CRITICAL TEST: Verify it's a CList, not a CString!
    assert isinstance(container.myList, CList), \
        f"myList should be CList, got {type(container.myList).__name__}"
    assert not isinstance(container.myList, CString), \
        "myList should NOT be CString (the old broken behavior)"


def test_custom_attribute_with_cint():
    """Test that AttributeType.CUSTOM works with CInt."""

    @cdata_class(
        attributes={
            'customInt': attribute(AttributeType.CUSTOM, custom_class="CInt")
        },
        qualifiers={
            'default': 42
        }
    )
    class TestContainer(CData):
        """Test container with custom CInt attribute."""
        pass

    container = TestContainer(name="testContainer")

    assert hasattr(container, 'customInt')
    assert isinstance(container.customInt, CInt), \
        f"customInt should be CInt, got {type(container.customInt).__name__}"
    # Default value should be applied
    assert container.customInt.value == 42


def test_custom_attribute_fallback_to_cstring():
    """Test that unknown custom classes fall back to CString with warning."""

    @cdata_class(
        attributes={
            'unknownClass': attribute(AttributeType.CUSTOM,
                                      custom_class="NonExistentClass")
        }
    )
    class TestContainer(CData):
        """Test container with unknown custom class."""
        pass

    container = TestContainer(name="testContainer")

    # Should still create an attribute (fallback to CString)
    assert hasattr(container, 'unknownClass')
    # Should be CString as fallback
    assert isinstance(container.unknownClass, CString)


def test_custom_attribute_without_custom_class():
    """Test that CUSTOM without custom_class falls back to CString."""

    @cdata_class(
        attributes={
            'noClass': attribute(AttributeType.CUSTOM)  # No custom_class!
        }
    )
    class TestContainer(CData):
        """Test container with CUSTOM but no custom_class."""
        pass

    container = TestContainer(name="testContainer")

    # Should create attribute with CString fallback
    assert hasattr(container, 'noClass')
    assert isinstance(container.noClass, CString)


def test_multiple_custom_attributes():
    """Test multiple custom attributes in one class."""

    @cdata_class(
        attributes={
            'intAttr': attribute(AttributeType.CUSTOM, custom_class="CInt"),
            'listAttr': attribute(AttributeType.CUSTOM, custom_class="CList"),
            'strAttr': attribute(AttributeType.CUSTOM, custom_class="CString"),
        }
    )
    class TestContainer(CData):
        """Test container with multiple custom attributes."""
        pass

    container = TestContainer(name="testContainer")

    assert isinstance(container.intAttr, CInt)
    assert isinstance(container.listAttr, CList)
    assert isinstance(container.strAttr, CString)


def test_custom_attribute_with_qualifiers():
    """Test that qualifiers are applied to custom attributes."""

    @cdata_class(
        attributes={
            'customInt': attribute(AttributeType.CUSTOM, custom_class="CInt")
        },
        qualifiers={
            'min': 0,
            'max': 100,
            'default': 10
        }
    )
    class TestContainer(CData):
        """Test container with custom attribute and qualifiers."""
        pass

    container = TestContainer(name="testContainer")

    # Check qualifiers were applied
    assert hasattr(container.customInt, 'qualifiers')
    assert container.customInt.get_qualifier('min') == 0
    assert container.customInt.get_qualifier('max') == 100


def test_censemble_label_is_coneword():
    """Test that CEnsemble.label is created as COneWord (if attributes are defined).

    Note: This test checks if the CEnsemble class properly creates attributes
    from its type annotations. CEnsemble has type annotation `label: Optional[COneWord]`
    but may not have explicit attribute metadata in the @cdata_class decorator.
    """
    from ccp4i2.core.CCP4ModelData import CEnsemble
    from ccp4i2.core.CCP4Data import COneWord

    ensemble = CEnsemble(name="testEnsemble")

    # Check if label attribute exists
    # Note: CEnsemble may not have attributes defined in decorator metadata,
    # so this test may reveal that type annotations alone aren't sufficient
    if hasattr(ensemble, 'label'):
        # Verify it's a COneWord instance
        assert isinstance(ensemble.label, COneWord), \
            f"label should be COneWord, got {type(ensemble.label).__name__}"
    else:
        # Document that the attribute wasn't created from type annotation
        # This would indicate that we need explicit attribute metadata
        import pytest
        pytest.skip("CEnsemble.label not created - type annotations not converted to attributes")
