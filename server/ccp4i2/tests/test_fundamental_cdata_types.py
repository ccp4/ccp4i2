"""
Pytest tests for fundamental CData types
"""
import pytest

from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean, CList
from ccp4i2.core.base_object.base_classes import CContainer

def test_cint_basic():
    x = CInt(5, qualifiers={'min': 0, 'max': 10})
    assert int(x) == 5
    assert str(x) == '5'
    assert x == 5
    assert x + 2 == 7
    x.value = 8
    assert int(x) == 8
    # Qualifiers
    assert x.qualifiers['min'] == 0
    assert x.qualifiers['max'] == 10

def test_cfloat_basic():
    y = CFloat(3.14, qualifiers={'min': 0.0, 'max': 10.0})
    assert float(y) == pytest.approx(3.14)
    assert str(y) == '3.14'
    assert y == 3.14
    assert y + 2.0 == pytest.approx(5.14)
    # Qualifiers
    assert y.qualifiers['min'] == 0.0
    assert y.qualifiers['max'] == 10.0

def test_cstring_basic():
    s = CString('hello', qualifiers={'allowedChars': 'abcdehlor'})
    assert str(s) == 'hello'
    assert f"{s}" == 'hello'
    assert s == 'hello'
    assert s + ' world' == 'hello world'
    # Qualifiers
    assert s.qualifiers['allowedChars'] == 'abcdehlor'

def test_cboolean_basic():
    b = CBoolean(True, qualifiers={'default': False})
    assert bool(b) is True
    assert str(b) == 'True'
    assert b == True
    # Qualifiers
    assert b.qualifiers['default'] is False

def test_clist_basic():
    l = CList([CInt(1), CInt(2)], qualifiers={'minLength': 1})
    assert len(l) == 2
    assert int(l[0]) == 1
    l.append(CInt(3))
    assert len(l) == 3
    # Qualifiers
    assert l.qualifiers['minLength'] == 1

def test_ccontainer_basic():
    # Create container and add items via attribute assignment
    c = CContainer(qualifiers={'containerType': 'mixed'})
    c.item1 = CInt(1, name='item1')
    c.item2 = CString('a', name='item2')
    assert len(c) == 2
    assert isinstance(c[0], CInt)
    assert isinstance(c[1], CString)

    # Add another item
    c.item3 = CFloat(2.5, name='item3')
    assert len(c) == 3

    # Qualifiers
    assert c.qualifiers['containerType'] == 'mixed'
