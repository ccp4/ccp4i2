#!/usr/bin/env python3
# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""Test that CContainer.copyData() works correctly."""
from ccp4i2.core.base_object.ccontainer import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt, CString, CBoolean


def _make_source_container():
    """Create a source container with test data."""
    source = CContainer(name="source")

    source.count = CInt(42, name="count")
    source.count.set_parent(source)
    source._data_order.append("count")

    source.message = CString("hello", name="message")
    source.message.set_parent(source)
    source._data_order.append("message")

    source.flag = CBoolean(True, name="flag")
    source.flag.set_parent(source)
    source._data_order.append("flag")

    return source


def test_copydata_full():
    """Test that copyData copies all fields from source to destination."""
    source = _make_source_container()

    dest = CContainer(name="dest")
    dest.count = CInt(0, name="count")
    dest.count.set_parent(dest)
    dest._data_order.append("count")

    dest.copyData(source)

    assert dest.count.value == 42, f"Expected count=42, got {dest.count.value}"
    assert dest.message.value == 'hello', f"Expected message='hello', got {dest.message.value}"
    assert dest.flag.value is True, f"Expected flag=True, got {dest.flag.value}"


def test_copydata_selective():
    """Test that copyData with dataList only copies specified fields."""
    source = _make_source_container()

    dest = CContainer(name="dest2")
    dest.copyData(source, dataList=["count", "message"])

    assert len(dest.dataOrder()) == 2, \
        f"Expected 2 items in dest, got {len(dest.dataOrder())}"

    has_flag = any(
        (hasattr(child, 'name') and child.name == 'flag') or
        (hasattr(child, 'objectName') and child.objectName() == 'flag')
        for child in dest.children()
    )
    assert not has_flag, "dest should NOT have 'flag' after selective copy"
