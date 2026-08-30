"""
Who is a child, for how long, and what happens when one dies.

Tier 5, written before `_children` changes structure. Today it is a
`Set[weakref.ReferenceType]`; the reevaluation may make it ordered, which is
the point at which these properties stop being free.

They are free today for reasons worth naming, because a different structure
loses them silently:

  - a `set` deduplicates, so registering the same child twice is harmless
  - `weakref.ref` hashes and compares by *referent identity*, so nothing ever
    names an integer, and a dead child cannot be confused with a live one
  - membership and iteration are the same structure, so they cannot disagree

An id-keyed dict would keep the first and third and lose the second: CPython
reuses `id()` values after garbage collection, so between a child dying and
`_cleanup_dead_children` running, a new object can take a dead one's id. These
tests pin the behaviour that must survive whatever is chosen.

Pure Python -- no CCP4 binaries needed.
"""
import gc
import threading

import pytest

from ccp4i2.core.base_object.hierarchy_system import HierarchicalObject


def _child(parent, name):
    child = HierarchicalObject(name=name)
    child.set_parent(parent)
    return child


def _names(parent):
    return [c._name for c in parent.children()]


# --- membership ------------------------------------------------------------

def test_registering_the_same_child_twice_does_not_double_it():
    """Free from a set; not free from a list."""
    parent = HierarchicalObject(name='root')
    child = _child(parent, 'only')
    child.set_parent(parent)

    assert _names(parent) == ['only']
    assert len([c for c in parent.children() if c is child]) == 1


def test_a_name_identifies_at_most_one_child():
    """What `test_hash_collision_fix.py` was really protecting.

    That test asserts `hash(obj1) != hash(obj2)`, pinning the *mechanism* that
    let two same-named children coexist. The guarantee underneath was that
    neither silently vanished. It is now a stronger one: the situation cannot
    arise, because a name identifies at most one child and a second assignment
    replaces the first, as a dict does.

    That is a deliberate narrowing. Coexistence was never usable --- the
    displaced child could not be addressed, serialised, or named in a def.xml
    --- and it cost a real lifetime bug in CList, where every item was '[0]'
    until Defect C numbered them.
    """
    parent = HierarchicalObject(name='root')
    kept = [_child(parent, '[0]') for _ in range(3)]

    returned = parent.children()
    assert len(returned) == 1
    assert returned[0] is kept[-1]
    assert all(c._parent_ref is None for c in kept[:-1])


def test_a_parent_keeps_its_children_alive():
    """Lifetime is governed by the strong reference, not by the weak ones.

    `_add_child` puts the child in `__dict__` keyed by name, so a child with
    its own name cannot be collected while its parent lives. The weakref set
    does not manage lifetime; it manages membership. Worth stating plainly,
    because the weakness is what the structure looks like it is for.

    That strong reference used to live in a separate `_child_storage` dict
    keyed identically --- 19,529 entries, every one already in `__dict__`
    under the same key holding the same object. `__dict__` is now the only
    copy, which is also what ordinary attribute lookup reads.
    """
    parent = HierarchicalObject(name='root')
    _child(parent, 'stays')
    doomed = _child(parent, 'goes')
    del doomed
    gc.collect()

    assert _names(parent) == ['stays', 'goes']
    strong = [k for k, v in vars(parent).items() if v in (parent.children())]
    assert strong == ['stays', 'goes']


def test_replacing_a_child_by_name_releases_the_old_one():
    """The hazard this used to record is now impossible to reach.

    It read: the strong reference is keyed by name, so of several children
    sharing one only the last is held, and dropping the outside references
    collects the rest. That was every CList item before Defect C, when they
    were all named '[0]' --- an item could be collected while still nominally
    a child.

    Replacement removes the hazard rather than documenting it. A displaced
    child is detached deliberately, so being collected afterwards is correct
    and not a leak of something still in use.
    """
    parent = HierarchicalObject(name='root')
    siblings = [_child(parent, 'same') for _ in range(3)]
    assert len(parent.children()) == 1

    survivor = parent.children()[0]
    assert survivor is siblings[-1]

    del siblings
    gc.collect()

    assert parent.children() == [survivor]


def test_a_detached_child_is_no_longer_kept_alive():
    """The counterpart: letting go must actually let go, or removal leaks."""
    parent = HierarchicalObject(name='root')
    child = _child(parent, 'temporary')
    child.set_parent(None)
    del child
    gc.collect()

    assert parent.children() == []
    assert 'temporary' not in vars(parent)


def test_a_new_object_does_not_inherit_a_dead_ones_place():
    """The id-reuse hazard, stated as a test.

    CPython reuses id() values. A structure keyed on id() can let a fresh
    object land on a dead one's key --- shadowing it, or appearing to be a
    child of a parent it was never given to. Creating and dropping many
    objects makes a collision likely if the keying is by id.
    """
    parent = HierarchicalObject(name='root')
    kept = _child(parent, 'stays')

    for i in range(200):
        transient = HierarchicalObject(name=f'transient_{i}')
        transient.set_parent(parent)
        transient.set_parent(None)   # detach, so nothing keeps it alive
        del transient
    gc.collect()

    strangers = [c for c in parent.children() if c is not kept]
    assert not strangers, \
        f'{len(strangers)} objects are children of a parent that never kept them'
    assert _names(parent) == ['stays']


def test_a_child_outliving_its_parent_does_not_keep_it_alive():
    """The reference from child to parent is weak, and must stay weak."""
    parent = HierarchicalObject(name='root')
    child = HierarchicalObject(name='child')
    child.set_parent(parent)
    del parent
    gc.collect()

    assert child.parent() is None


# --- re-parenting ----------------------------------------------------------

def test_a_child_moved_to_another_parent_leaves_the_first():
    first = HierarchicalObject(name='first')
    second = HierarchicalObject(name='second')
    child = _child(first, 'movable')

    child.set_parent(second)

    assert _names(first) == []
    assert _names(second) == ['movable']
    assert child.parent() is second


def test_detaching_leaves_no_trace_behind():
    parent = HierarchicalObject(name='root')
    child = _child(parent, 'temporary')

    child.set_parent(None)

    assert parent.children() == []
    assert parent.find_child('temporary') is None
    assert child.parent() is None


def test_detaching_and_reattaching_restores_one_child_not_two():
    parent = HierarchicalObject(name='root')
    child = _child(parent, 'yo-yo')

    child.set_parent(None)
    child.set_parent(parent)

    assert _names(parent) == ['yo-yo']
    assert parent.find_child('yo-yo') is child


def test_removing_the_first_of_several_keeps_the_others_in_order():
    parent = HierarchicalObject(name='root')
    kept = [_child(parent, n) for n in ('a', 'b', 'c')]
    kept[0].set_parent(None)

    assert _names(parent) == ['b', 'c']


# --- membership and iteration cannot disagree ------------------------------

def test_find_child_agrees_with_children():
    parent = HierarchicalObject(name='root')
    kept = [_child(parent, n) for n in ('a', 'b', 'c')]

    for child in parent.children():
        assert parent.find_child(child._name) is child
    assert len(kept) == 3


def test_children_is_stable_while_nothing_changes():
    parent = HierarchicalObject(name='root')
    kept = [_child(parent, f'c{i}') for i in range(10)]
    first = _names(parent)
    for _ in range(5):
        assert _names(parent) == first
    assert len(kept) == 10


# --- concurrency -----------------------------------------------------------

def test_reading_children_while_they_are_added_does_not_raise():
    """`children()` and `_add_child` are guarded by the same lock. A structure
    that iterates without holding it raises "changed size during iteration"
    under exactly this load."""
    parent = HierarchicalObject(name='root')
    kept = []
    errors = []

    def add():
        try:
            for i in range(150):
                child = HierarchicalObject(name=f'added_{i}')
                child.set_parent(parent)
                kept.append(child)
        except Exception as err:  # pragma: no cover - the failure being tested
            errors.append(err)

    def read():
        try:
            for _ in range(150):
                parent.children()
        except Exception as err:  # pragma: no cover
            errors.append(err)

    threads = [threading.Thread(target=add), threading.Thread(target=read)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert not errors, f'{type(errors[0]).__name__}: {errors[0]}'
