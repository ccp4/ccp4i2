"""
Test that the hash collision bug is fixed.

This test verifies that multiple CInt/CFloat objects with the same value
can coexist as children of the same container.
"""

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat


def test_multiple_cints_with_same_value():
    """Test that multiple CInt objects with value=10 can be children of same container."""

    container = CContainer(name="test_container")

    # Create multiple CInt objects with value=10
    ncycles = CInt(10, name="NCYCLES")
    another_ten = CInt(10, name="ANOTHER_TEN")
    yet_another_ten = CInt(10, name="YET_ANOTHER_TEN")

    # Add them to container
    container.NCYCLES = ncycles
    container.ANOTHER_TEN = another_ten
    container.YET_ANOTHER_TEN = yet_another_ten

    # Verify all are in children list
    child_names = {c.objectName() for c in container.children()}

    print(f"Container has {len(child_names)} children")
    print(f"Children: {sorted(child_names)}")

    assert len(child_names) == 3, f"Expected 3 children, got {len(child_names)}"
    assert "NCYCLES" in child_names, "NCYCLES should be in children"
    assert "ANOTHER_TEN" in child_names, "ANOTHER_TEN should be in children"
    assert "YET_ANOTHER_TEN" in child_names, "YET_ANOTHER_TEN should be in children"

    print("✅ Multiple CInt objects with same value work correctly!")


def test_multiple_cfloats_with_same_value():
    """Test that multiple CFloat objects with value=1.5 can be children of same container."""

    container = CContainer(name="test_container")

    # Create multiple CFloat objects with value=1.5
    weight1 = CFloat(1.5, name="WEIGHT1")
    weight2 = CFloat(1.5, name="WEIGHT2")
    weight3 = CFloat(1.5, name="WEIGHT3")

    # Add them to container
    container.WEIGHT1 = weight1
    container.WEIGHT2 = weight2
    container.WEIGHT3 = weight3

    # Verify all are in children list
    child_names = {c.objectName() for c in container.children()}

    print(f"\nContainer has {len(child_names)} children")
    print(f"Children: {sorted(child_names)}")

    assert len(child_names) == 3, f"Expected 3 children, got {len(child_names)}"
    assert "WEIGHT1" in child_names, "WEIGHT1 should be in children"
    assert "WEIGHT2" in child_names, "WEIGHT2 should be in children"
    assert "WEIGHT3" in child_names, "WEIGHT3 should be in children"

    print("✅ Multiple CFloat objects with same value work correctly!")


def test_same_valued_objects_stay_distinct():
    """The property this file exists for, stated as a property.

    It used to assert `hash(obj1) != hash(obj2)` --- which pins the
    *mechanism*, not the guarantee. `_children` is a set of weak references,
    and a weakref hashes by its referent, so a value-based hash would make two
    same-valued children collide and one would silently vanish. That is why
    the identity hash exists.

    But if `_children` ever stops being hash-based, an assertion about hashes
    fails for a reason that has nothing to do with what it protects, and the
    obvious repair is to delete it --- taking the real guarantee with it. So
    the guarantee is asserted directly: same value, still two distinct
    children, both reachable.

    See also tests/unit/conformance/test_identity_and_lifetime.py.
    """
    container = CContainer(name="test_container")
    first = CInt(10, name="first")
    second = CInt(10, name="second")
    container.first = first
    container.second = second

    assert first == second, "equal by value, which is the interesting case"
    assert first is not second

    children = container.children()
    assert len({id(c) for c in children}) == len(children) == 2
    assert container.find_child("first") is first
    assert container.find_child("second") is second


if __name__ == "__main__":
    try:
        test_multiple_cints_with_same_value()
        test_multiple_cfloats_with_same_value()
        test_same_valued_objects_stay_distinct()
        print("\n🎉 All hash collision fix tests passed!")
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def test_a_cdata_value_is_not_a_dict_key():
    """Equal to a str, but not interchangeable with one as a key.

    The identity hash the tests above require means a CString cannot be used to
    look up a dict keyed by plain strings: it compares equal and hashes
    differently, so the lookup misses and returns the default. A wrapper author
    who writes ``ext_map.get(par.OUTPUT_TYPE, ".map")`` gets the default every
    time --- which is how every nucleofind RAW run came to look for a file
    nucleofind had not written.

    Convert at the call site: ``ext_map.get(str(par.OUTPUT_TYPE), ...)``.
    """
    from ccp4i2.core.base_object.fundamental_types import CString

    value = CString()
    value.set("RAW")
    table = {"RAW": ".raw.map"}

    assert value == "RAW"
    assert table.get(value) is None, "equal, but not the same key"
    assert table.get(str(value)) == ".raw.map"
