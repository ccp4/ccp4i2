import pytest

from ccp4i2.core.base_object.base_classes import CContainer, CDataFile
from ccp4i2.core.base_object.fundamental_types import CInt, CList

class TestExample:
    @classmethod
    def setup_class(cls):
        # Class-level setup (runs once before all tests)
        cls.shared_resource = "initialized"

    def test_ccontainer_inheritance(self):
        c = CContainer()
        # Check that CContainer is indeed a subclass of the base CContainer
        assert isinstance(c, CContainer)

    def test_ccontainer_methods(self):
        c = CContainer(name="container")
        # Test adding items via attribute assignment (new API)
        b = CInt(10, name="test_int")
        d = CDataFile(name="test_file")
        d.update({"baseName": "changed"})

        # Add items as named attributes
        c.test_int = b
        c.test_file = d

        # Access via children() method
        items = c.children()
        assert len(items) == 2
        assert b.objectName() == "test_int"
        assert b.object_path() == "container.test_int"
        assert d.objectName() == "test_file"
        assert d.object_path() == "container.test_file"

        # Add another item
        l = CList(name="test_list")
        c.test_list = l
        items = c.children()
        assert len(items) == 3

        # Access by index
        the_list: CList = c[2]
        assert the_list.objectName() == "test_list"
        assert the_list.object_path() == "container.test_list"
        assert isinstance(the_list, CList)
        assert l is the_list

        # Test qualifiers
        print(b.get_qualifier("min"))
        b.set_qualifier("min", 2)
        b.set_qualifier("max", 100)
        with pytest.raises(ValueError):
            b.value = 120  # This should raise a ValueError due to the max qualifier
