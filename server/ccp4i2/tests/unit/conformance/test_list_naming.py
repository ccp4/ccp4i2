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


def test_the_two_are_deliberately_different():
    """Not unified, and on purpose.

    dataOrder() is an ordering of positions, used for serialisation; _name is
    a handle, used to build object paths like DICT_LIST[0]. Making dataOrder()
    emit '[0]' would change the shape of every params.xml for no gain, so
    Defect C left the difference and made the *handle* work instead.
    """
    lst = _list_of(3)
    assert lst.dataOrder() == ['0', '1', '2']
    assert [i._name for i in lst._items] == ['[0]', '[1]', '[2]']


def test_an_item_is_reachable_by_its_own_name():
    """It was not: the cache held the name the item had at registration."""
    lst = _list_of(2)
    assert hasattr(lst, '[0]')
    assert not hasattr(lst, '0'), 'dataOrder positions are not handles'


def test_indexing_is_the_route_that_works():
    lst = _list_of(3)
    assert [int(lst[i]) for i in range(3)] == [0, 1, 2]
    assert all(lst[i] is lst._items[i] for i in range(3))


def test_removing_an_item_renumbers_the_survivors():
    lst = _list_of(3)
    lst.pop(0)
    assert [i._name for i in lst._items] == ['[0]', '[1]']
    assert lst.dataOrder() == ['0', '1']


def test_a_removed_item_stops_being_a_child():
    """pop used to renumber the survivors and keep the departed one, so a
    three-item list held two items and three children, two named '[0]'."""
    lst = _list_of(3)
    lst.pop(0)
    assert [c._name for c in lst.children()] == [i._name for i in lst._items]
    assert len(lst.children()) == 2


def test_find_child_finds_an_item_by_its_name():
    lst = _list_of(2)
    assert lst.find_child('[0]') is lst._items[0]
    assert lst.find_child('[1]') is lst._items[1]


def test_assigning_by_index_works_at_all():
    """It did not: `self.name` does not exist on CList, so every CData
    assigned to a list index raised AttributeError. Nothing in the tree
    exercised the path, so nothing said so until the mutation tier did."""
    lst = _list_of(3)
    lst[1] = CInt(99)
    assert int(lst[1]) == 99


def test_assignment_names_the_item_the_way_append_does():
    lst = _list_of(3)
    lst[1] = CInt(99)
    assert lst[1]._name == '[1]', 'it used to write "LIST[1]", a third convention'


def test_assignment_detaches_what_it_displaces():
    lst = _list_of(3)
    displaced = lst[1]
    lst[1] = CInt(99)
    assert not any(c is displaced for c in lst.children())
    assert len(lst.children()) == 3


def test_the_keys_are_the_names_the_items_actually_have():
    """The root of all three: keys were the name held at registration.

    There is now one structure rather than three --- children live in
    ``__dict__``, keyed by name --- so the question is whether renumbering
    after a pop keeps the keys and the items' own names in step.
    """
    lst = _list_of(3)
    keyed = [k for k, v in vars(lst).items() if not k.startswith('_')]
    assert keyed == ['[0]', '[1]', '[2]']
    assert [c._name for c in lst.children()] == ['[0]', '[1]', '[2]']

    lst.pop(0)
    keyed = [k for k, v in vars(lst).items() if not k.startswith('_')]
    assert keyed == ['[0]', '[1]']
    assert [c._name for c in lst.children()] == ['[0]', '[1]']
