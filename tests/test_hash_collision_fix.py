"""
Test that the hash collision bug is fixed.

This test verifies that multiple CInt/CFloat objects with the same value
can coexist as children of the same container.
"""

import sys
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "server"))

# Set CCP4I2_ROOT for plugin discovery
os.environ["CCP4I2_ROOT"] = str(project_root)

from core.CCP4Container import CContainer
from core.base_object.fundamental_types import CInt, CFloat


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


def test_hash_is_identity_based():
    """Verify that hash is based on object identity, not value."""

    obj1 = CInt(10, name="obj1")
    obj2 = CInt(10, name="obj2")

    # Same value, different objects
    assert obj1.value == obj2.value, "Values should be equal"
    # Note: CInt has __eq__ that compares by value, which is OK
    # The important thing is that hash is identity-based
    assert hash(obj1) != hash(obj2), "Hashes should be different (identity-based)"
    assert id(obj1) != id(obj2), "Objects should have different identities"

    print("\n✅ Hash is identity-based (not value-based)!")


if __name__ == "__main__":
    try:
        test_multiple_cints_with_same_value()
        test_multiple_cfloats_with_same_value()
        test_hash_is_identity_based()
        print("\n🎉 All hash collision fix tests passed!")
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
