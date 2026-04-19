"""
Ensure CList.set() wraps plain primitive values into the list's declared
subItem class, so the server's canonical form round-trips through the
JSON encoder with proper _class/_objectPath metadata.

Regression test for a bug where deleting an entry from a CList of CString
(e.g. SGALT_TEST in phaser tasks) caused the client to receive an updated
list whose children were bare Python strings — no _objectPath, so the
client's lookup table couldn't render them until a full refetch.
"""

from ccp4i2.core.base_object.fundamental_types import (
    CList,
    CString,
    CInt,
    CFloat,
    CBoolean,
)


def test_set_wraps_plain_strings_into_cstring():
    cl = CList(name="SGALT_TEST")
    cl.set_qualifier("subItem", {"class": CString, "qualifiers": {}})

    cl.set(["P 1", "P 2", "P 21 21 21"])

    assert len(cl) == 3
    for item in cl:
        assert isinstance(item, CString), (
            f"Expected CString, got {type(item).__name__}"
        )
    assert cl[0].value == "P 1"
    assert cl[1].value == "P 2"
    assert cl[2].value == "P 21 21 21"


def test_set_wraps_plain_ints_into_cint():
    cl = CList(name="NUMBERS")
    cl.set_qualifier("subItem", {"class": CInt, "qualifiers": {}})

    cl.set([1, 2, 3])

    assert len(cl) == 3
    for item in cl:
        assert isinstance(item, CInt)
    assert [it.value for it in cl] == [1, 2, 3]


def test_set_wraps_plain_floats_into_cfloat():
    cl = CList(name="VALUES")
    cl.set_qualifier("subItem", {"class": CFloat, "qualifiers": {}})

    cl.set([0.5, 1.5, 2.5])

    assert len(cl) == 3
    for item in cl:
        assert isinstance(item, CFloat)
    assert [it.value for it in cl] == [0.5, 1.5, 2.5]


def test_set_wraps_plain_booleans_into_cboolean():
    cl = CList(name="FLAGS")
    cl.set_qualifier("subItem", {"class": CBoolean, "qualifiers": {}})

    cl.set([True, False, True])

    assert len(cl) == 3
    for item in cl:
        assert isinstance(item, CBoolean)
    assert [it.value for it in cl] == [True, False, True]


def test_set_preserves_cdata_items_unchanged_shape():
    """A list of CString set with existing CString instances still yields CStrings."""
    cl = CList(name="SGALT_TEST")
    cl.set_qualifier("subItem", {"class": CString, "qualifiers": {}})

    a = CString("P 1")
    b = CString("P 2")
    cl.set([a, b])

    assert len(cl) == 2
    for item in cl:
        assert isinstance(item, CString)
    assert cl[0].value == "P 1"
    assert cl[1].value == "P 2"


def test_set_still_wraps_dicts_into_subitem_class():
    """Pre-existing dict-wrapping behaviour must continue to work."""
    cl = CList(name="SGALT_TEST")
    cl.set_qualifier("subItem", {"class": CString, "qualifiers": {}})

    cl.set([{"_value": "P 1"}, {"_value": "P 2"}])

    assert len(cl) == 2
    for item in cl:
        assert isinstance(item, CString)


def test_set_replaces_previous_contents():
    """set() replaces — the classic 'delete one and rewrite the array' flow."""
    cl = CList(name="SGALT_TEST")
    cl.set_qualifier("subItem", {"class": CString, "qualifiers": {}})

    cl.set(["P 1", "P 2", "P 21 21 21"])
    assert len(cl) == 3

    # Simulate a client-side delete: rewrite with one fewer entry
    cl.set(["P 1", "P 21 21 21"])

    assert len(cl) == 2
    for item in cl:
        assert isinstance(item, CString)
    assert cl[0].value == "P 1"
    assert cl[1].value == "P 21 21 21"
