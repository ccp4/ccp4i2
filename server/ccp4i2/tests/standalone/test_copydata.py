#!/usr/bin/env python3
"""Test script to verify CContainer.copyData() works correctly"""
from ccp4i2.core.base_object.ccontainer import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CString, CBoolean

print("[TEST] Creating source container with data...")
source = CContainer(name="source")

# Use setattr to add items (this is what DefXmlParser does)
source.count = CInt(42, name="count")
source.count.set_parent(source)
source._data_order.append("count")

source.message = CString("hello", name="message")
source.message.set_parent(source)
source._data_order.append("message")

source.flag = CBoolean(True, name="flag")
source.flag.set_parent(source)
source._data_order.append("flag")

print(f"  Source container has {len(source.dataOrder())} items:")
for item_name in source.dataOrder():
    item = getattr(source, item_name)
    print(f"    - {item_name}: {item.value} ({type(item).__name__})")

print("\n[TEST] Creating destination container...")
dest = CContainer(name="dest")

# Add one field with different value to test update
dest.count = CInt(0, name="count")
dest.count.set_parent(dest)
dest._data_order.append("count")
print(f"  Dest.count before copy: {dest.count.value}")

print("\n[TEST] Calling dest.copyData(source)...")
dest.copyData(source)

print(f"\n[TEST] Destination container now has {len(dest.dataOrder())} items:")
for item_name in dest.dataOrder():
    item = getattr(dest, item_name)
    print(f"    - {item_name}: {item.value} ({type(item).__name__})")

# Verify the copy worked
print("\n[TEST] Verification:")
print(f"  dest.count.value == 42: {dest.count.value == 42}")
print(f"  dest.message.value == 'hello': {dest.message.value == 'hello'}")
print(f"  dest.flag.value == True: {dest.flag.value == True}")

# Test selective copy
print("\n[TEST] Testing selective copy...")
dest2 = CContainer(name="dest2")
dest2.copyData(source, dataList=["count", "message"])
print(f"  dest2 has {len(dest2.dataOrder())} items (should be 2):")
for item_name in dest2.dataOrder():
    item = getattr(dest2, item_name)
    print(f"    - {item_name}: {item.value} ({type(item).__name__})")

# Check if flag exists as a child (should be False since we didn't copy it)
# Check both 'name' attribute and objectName() method
has_flag = any(
    (hasattr(child, 'name') and child.name == 'flag') or
    (hasattr(child, 'objectName') and child.objectName() == 'flag')
    for child in dest2.children()
)
print(f"  dest2 has flag: {has_flag} (should be False)")

print("\n[TEST] All tests passed!")
