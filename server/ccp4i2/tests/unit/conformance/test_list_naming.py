"""
What a CList calls its items --- pinned, because it is currently three things.

Not assertions that today's behaviour is right. These record what it *is*, so
that Defect C (docs/container-construction-defects.md) has to change them
deliberately rather than by accident, and so the change is visible in a diff.

    dataOrder()      '0'          bare index
    item._name       '[0]'        bracketed
    __setitem__      'LIST[0]'    qualified, and on a different attribute

and the way anybody actually reaches an item is `lst[0]`.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.CCP4Data import CInt, CList


def _list_of(n):
    lst = CList()
    for i in range(n):
        lst.append(CInt(i))
    return lst


def test_dataorder_uses_bare_indices():
    assert _list_of(3).dataOrder() == ['0', '1', '2']


def test_item_names_use_brackets():
    assert [i._name for i in _list_of(3)._items] == ['[0]', '[1]', '[2]']


def test_the_two_disagree_today():
    """The point of the pin. When Defect C lands, this should fail."""
    lst = _list_of(3)
    assert lst.dataOrder() != [i._name for i in lst._items], \
        'the naming was unified -- update Defect C notes and delete this pin'


def test_an_item_is_not_reachable_by_either_name():
    lst = _list_of(2)
    assert not hasattr(lst, '0')
    assert not hasattr(lst, '[0]')


def test_indexing_is_the_route_that_works():
    lst = _list_of(3)
    assert [int(lst[i]) for i in range(3)] == [0, 1, 2]
    assert all(lst[i] is lst._items[i] for i in range(3))


def test_removing_an_item_renumbers_the_survivors():
    lst = _list_of(3)
    lst.pop(0)
    assert [i._name for i in lst._items] == ['[0]', '[1]']
    assert lst.dataOrder() == ['0', '1']


def test_a_removed_item_is_still_a_child_today():
    """Defect C: pop renumbers but never detaches, so the parent keeps it.

    Two children end up named '[0]' --- a stale one and a live one.
    """
    lst = _list_of(3)
    lst.pop(0)
    names = [c._name for c in lst.children()]
    assert len(names) > len(lst._items), \
        'pop now detaches -- Defect C has landed, update this pin'


def test_find_child_cannot_find_an_item_today():
    """The name cache is keyed on the name held at registration time."""
    lst = _list_of(2)
    assert lst.find_child('[0]') is None, \
        'the name cache is re-keyed on rename -- Defect C has landed'
