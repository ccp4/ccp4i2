"""
Signals belong to the thing that has subscribers.

Every HierarchicalObject used to build a SignalManager and four signals ---
destroyed, parent_changed, child_added, child_removed. Nothing in either
process ever connected to any of them: `connectSignal` is used for three names
in the whole tree (`finished` 16 times, `jobId` once, `progressUpdated` once),
all plugin lifecycle, none data change. The frontend learns about changes by
polling the REST API --- useSWR in 18 files, no EventSource, no WebSocket of
ours, no event-stream endpoint.

A full registry build made about 115,000 signal objects for no subscriber.

Pure Python -- no CCP4 binaries needed.
"""
import tempfile

import pytest

from ccp4i2.core.base_object.hierarchy_system import HierarchicalObject


def test_an_object_does_not_build_a_signal_manager_it_will_not_use():
    obj = HierarchicalObject(name='plain')
    assert obj._sigmgr is None


def test_the_four_signals_that_nobody_connected_are_gone():
    obj = HierarchicalObject(name='plain')
    for name in ('destroyed', 'parent_changed', 'child_added', 'child_removed'):
        assert not hasattr(obj, name), f'{name} is still created for every object'


def test_a_manager_appears_when_something_asks_for_a_signal():
    obj = HierarchicalObject(name='wants one')
    signal = obj.create_signal('something', dict)
    assert signal is not None
    assert obj._sigmgr is not None


def test_asking_twice_gives_the_same_manager():
    obj = HierarchicalObject(name='asks twice')
    obj.create_signal('first', dict)
    manager = obj._sigmgr
    obj.create_signal('second', dict)
    assert obj._sigmgr is manager


def test_building_a_tree_creates_no_signals():
    """The point of the change, stated over a hierarchy rather than one object."""
    root = HierarchicalObject(name='root')
    kept = []
    for i in range(20):
        child = HierarchicalObject(name=f'child_{i}')
        child.set_parent(root)
        kept.append(child)

    assert root._sigmgr is None
    assert all(c._sigmgr is None for c in kept)


def test_destroying_an_object_that_never_had_a_manager_is_fine():
    obj = HierarchicalObject(name='never asked')
    obj.destroy()          # must not raise on a manager that was never made


def test_destroying_an_object_that_did_have_one_is_also_fine():
    obj = HierarchicalObject(name='asked once')
    obj.create_signal('something', dict)
    obj.destroy()


# --- the one real consumer --------------------------------------------------

def test_a_plugin_still_has_its_finished_signal():
    """CPluginScript is the only thing that ever wanted a signal.

    Stubbing the four core signals broke every task build with exactly one
    error --- CPluginScript.__init__ creating 'finished' --- which is the
    whole argument for moving them here.
    """
    from ccp4i2.core.CCP4PluginScript import CPluginScript

    plugin = CPluginScript(workDirectory=tempfile.mkdtemp(), name='p')
    assert plugin.finished is not None
    assert plugin._sigmgr is not None, 'the plugin should have made one'


def test_a_pipeline_can_still_hear_a_subjob_finish():
    """16 connectSignal calls depend on this, and it is how pipelines chain."""
    from ccp4i2.core.CCP4PluginScript import CPluginScript

    parent = CPluginScript(workDirectory=tempfile.mkdtemp(), name='parent')
    child = CPluginScript(workDirectory=tempfile.mkdtemp(), name='child')
    heard = []

    parent.connectSignal(child, 'finished', lambda status: heard.append(status))
    child.finished.emit({'finishStatus': CPluginScript.SUCCEEDED})

    assert heard, 'the parent did not hear the child finish'
