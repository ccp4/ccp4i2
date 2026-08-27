"""
Sibling order is the order things were declared.

`_children` is a set of weak references, and a weakref hashes by its referent's
hash --- the default identity hash for these objects. So `children()` used to
iterate in allocation-address order: arbitrary per process, and arbitrary per
object within a process.

That was observable. i2run breaks a name collision by preferring the shortest
path, and falls through to `children()` order when the candidates are equally
deep. `i2run molrep_pipe --F_SIGF native.mtz` populated an *output* slot on
roughly half its invocations --- same command, same tree, different meaning run
to run --- and nothing downstream corrected it: `checkOutputData` leaves an
explicitly set output alone, and `validity()` filters outputData errors on the
grounds that outputs are set during execution.

Pure Python -- no CCP4 binaries needed.
"""
import gc

import pytest

from ccp4i2.core.base_object.hierarchy_system import HierarchicalObject


def _named(parent, name):
    child = HierarchicalObject(name=name)
    child.set_parent(parent)
    return child


def test_children_come_back_in_the_order_they_were_added():
    parent = HierarchicalObject(name='root')
    kept = [_named(parent, name) for name in ('inputData', 'outputData',
                                              'controlParameters', 'keywords')]
    assert [c._name for c in parent.children()] == \
        ['inputData', 'outputData', 'controlParameters', 'keywords']
    assert len(kept) == 4  # keep them alive


def test_the_order_is_the_same_every_time_it_is_asked():
    parent = HierarchicalObject(name='root')
    kept = [_named(parent, f'child_{i}') for i in range(20)]
    first = [c._name for c in parent.children()]
    for _ in range(5):
        assert [c._name for c in parent.children()] == first
    assert len(kept) == 20


def test_two_trees_built_the_same_way_agree():
    """The defect varied per object, so one tree could not vouch for another."""
    def build():
        parent = HierarchicalObject(name='root')
        kept = [_named(parent, name) for name in
                ('inputData', 'outputData', 'controlParameters')]
        return parent, kept

    a, _ka = build()
    b, _kb = build()
    assert [c._name for c in a.children()] == \
           [c._name for c in b.children()]


def test_a_child_whose_name_collided_is_still_returned():
    """CList items all take the name '[0]', so only the last is in the cache."""
    parent = HierarchicalObject(name='root')
    first = HierarchicalObject(name='[0]')
    first.set_parent(parent)
    second = HierarchicalObject(name='[0]')
    second.set_parent(parent)

    returned = parent.children()

    assert first in returned and second in returned
    assert len(returned) == 2


def test_a_collided_child_is_not_returned_twice():
    parent = HierarchicalObject(name='root')
    kept = [HierarchicalObject(name='[0]') for _ in range(3)]
    for child in kept:
        child.set_parent(parent)
    returned = parent.children()
    assert len({id(c) for c in returned}) == len(returned) == 3


def test_dead_children_are_left_out():
    parent = HierarchicalObject(name='root')
    kept = _named(parent, 'stays')
    doomed = _named(parent, 'goes')
    del doomed
    gc.collect()
    names = [c._name for c in parent.children()]
    assert 'stays' in names
    assert kept is not None


# --- the symptom, on a real task --------------------------------------------

def test_a_tasks_sections_come_out_in_declaration_order():
    """phaser_pipeline's top-level order differed between two instantiations."""
    import tempfile

    from ccp4i2.core.tasks import get_plugin_class

    plugin_class = get_plugin_class('phaser_pipeline')
    if plugin_class is None:
        pytest.skip('phaser_pipeline is not available in this environment')

    orders = []
    for _ in range(3):
        plugin = plugin_class(workDirectory=tempfile.mkdtemp(), name='t')
        orders.append(list(plugin.container.dataOrder()))

    assert orders[0] == orders[1] == orders[2]
    # def.xml declares inputData first, and 158 of 163 task trees do the same;
    # none declares outputData first. So "inputData wins" a name collision
    # without anyone having had to write that rule down.
    assert orders[0][0] == 'inputData'
    assert orders[0].index('inputData') < orders[0].index('outputData')


@pytest.mark.parametrize('task,flag', [
    ('molrep_pipe', 'F_SIGF'),
    ('xia2_ssx_reduce', 'DIALS_INTEGRATED'),
    ('adding_stats_to_mmcif_i2', 'FPHIOUT'),
    ('shelxeMR', 'PERFORMANCE'),
])
def test_a_bare_flag_means_the_same_thing_every_run(task, flag):
    """These four were measured splitting 7/5, 9/3, 8/4 and 11/1 across runs."""
    from ccp4i2.cli.i2run.i2run_components import KeywordExtractor

    landed = set()
    for _ in range(4):
        for keyword in KeywordExtractor.extract_from_task_name(task):
            if keyword.get('simpleName') == flag and \
                    keyword.get('isShortestForSimpleName'):
                landed.add(keyword['path'].split('.')[1])

    assert landed == {'inputData'}, f'--{flag} resolved to {landed}'


def test_a_removed_child_does_not_come_back_through_the_cache():
    """Membership is the set's business; the cache only decides order.

    `_remove_child` always drops the set entry, but clears the cache entry
    only when `_children_by_name[child._name]` still points at that child ---
    so a child renamed since it was registered leaves its cache entry
    orphaned. Reading the cache alone resurrected it: a CFloat phaser keyword
    came back holding a child named after itself, `dataOrder()` reported that
    child, and `phaser_MR.setKeywords` took the leaf for a sub-container and
    asked the float for a BOXS of its own.
    """
    parent = HierarchicalObject(name='root')
    kept = _named(parent, 'stays')
    renamed = _named(parent, 'registered_as_this')
    renamed._name = 'renamed_to_that'      # cache key is now stale
    renamed.set_parent(None)               # removed from the set, not the cache

    returned = parent.children()

    assert renamed not in returned, 'a removed child came back through the cache'
    assert kept in returned
