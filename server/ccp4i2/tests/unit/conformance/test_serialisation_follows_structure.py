"""
What crosses the API must be what the tree contains.

`children()` answers a structural question --- which objects are registered
beneath this one --- and returns live objects in registration order, from one
mechanism. `value_dict_for_object()` answers a serialisation question and
returns plain data; it is what `set_parameter` sends back, what the digest
reads, and what the job-66 file fix consults.

The second is currently implemented by *guessing at* the first, five ways in
turn --- declared metadata, `children()`, a `__dict__` scan, class `CONTENTS`,
then known property names --- each tried only if the previous produced nothing.
Only about 13% of objects are answered by `children()` at all, so agreement
between the two is presently a coincidence of five fallbacks rather than a
guarantee.

Measured across the registry it does hold: 287 nodes with children agree, and
the only two exceptions are CLists, which serialise to a list rather than a
dict --- correct, not a divergence. This pins that, so a change to how
children are keyed or stored cannot quietly make the API report something the
tree does not contain.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.lib.utils.parameters.value_dict import value_dict_for_object
from ccp4i2.tests.unit.conformance import harness


def _nodes(container, depth=0, seen=None, path='root'):
    """Every node with children, and its path."""
    seen = seen if seen is not None else set()
    if depth > 5 or id(container) in seen:
        return
    seen.add(id(container))
    if callable(getattr(container, 'children', None)) and container.children():
        yield path, container
    order = getattr(container, 'dataOrder', None)
    for name in (order() if callable(order) else []):
        child = getattr(container, name, None)
        if child is not None:
            yield from _nodes(child, depth + 1, seen, f'{path}.{name}')


@pytest.fixture(scope='module')
def nodes():
    found = []
    for task in harness.declared_tasks():
        try:
            plugin = harness.build(task)
        except Exception:
            continue
        if plugin is None:
            continue
        found.append((task, plugin, list(_nodes(plugin.container, path=task))))
    return found


def test_there_are_nodes_to_check(nodes):
    total = sum(len(n) for _t, _p, n in nodes)
    assert total > 200, f'only {total} nodes with children; not measuring much'


def test_serialisation_reports_exactly_the_children_the_tree_has(nodes):
    wrong = []
    for _task, _plugin, found in nodes:
        for path, obj in found:
            names = {c._name for c in obj.children()}
            serialised = value_dict_for_object(obj)
            if isinstance(serialised, list):
                # A CList serialises to a list; its items are positional, and
                # test_list_naming.py holds their naming.
                continue
            if not isinstance(serialised, dict):
                wrong.append(f'{path}: {type(obj).__name__} serialised to '
                             f'{type(serialised).__name__}')
                continue
            if set(serialised) != names:
                extra = sorted(set(serialised) - names)[:4]
                missing = sorted(names - set(serialised))[:4]
                wrong.append(f'{path}: value_dict-only={extra} children-only={missing}')
    assert not wrong, (f'{len(wrong)} nodes serialise something other than their '
                       f'children:\n  ' + '\n  '.join(wrong[:12]))


def test_no_signal_reaches_the_serialised_form(nodes):
    """Machinery is not data.

    The `__dict__` strategy once surfaced every object's four signals as if
    they were parameters, patched by a list of their names --- which went
    stale the moment the signals moved. Nothing named like a signal should
    appear at all.
    """
    leaked = []
    for _task, _plugin, found in nodes:
        for path, obj in found:
            serialised = value_dict_for_object(obj)
            if not isinstance(serialised, dict):
                continue
            for name in ('destroyed', 'parent_changed', 'child_added',
                         'child_removed', 'finished', 'progressUpdated'):
                if name in serialised:
                    leaked.append(f'{path}: {name}')
    assert not leaked, 'signals in the serialised form: ' + ', '.join(leaked[:8])


def test_serialisation_does_not_change_as_values_are_set():
    """The route through the five strategies depends on state; the answer
    should not. Setting a file path and an annotation must not change which
    keys are reported."""
    plugin = harness.build('freerflag')
    if plugin is None:
        pytest.skip('freerflag is not available in this environment')
    obj = plugin.container.inputData.F_SIGF

    fresh = sorted(value_dict_for_object(obj))
    obj.setFullPath('/tmp/some_reflections.mtz')
    after_path = sorted(value_dict_for_object(obj))
    obj.annotation.set('an annotation')
    after_annotation = sorted(value_dict_for_object(obj))

    assert fresh == after_path == after_annotation
