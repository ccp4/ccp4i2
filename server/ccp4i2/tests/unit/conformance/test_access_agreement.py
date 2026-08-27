"""
Every route to a parameter must reach the same object.

Tier 1 of the conformance harness (docs/container-construction-defects.md).
Constructs every registered task and offers each parameter through the
mechanisms plugins use, asserting the mechanisms agree with one another.

It is a whole-registry check by design: the defects it exists to catch ---
ghost containers, a leaf claiming children, a name cache disagreeing with
membership --- appeared in a handful of tasks out of 171 and in none of the
ones anybody thought to test by hand.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.tests.unit.conformance import harness


@pytest.fixture(scope='module')
def walked():
    """Walk every loadable task once; the tests read the result."""
    results, unloadable = {}, []
    for task in harness.declared_tasks():
        try:
            plugin = harness.build(task)
        except Exception as err:
            unloadable.append((task, f'{type(err).__name__}: {err}'))
            continue
        if plugin is None:
            unloadable.append((task, 'not in the registry'))
            continue
        # Hold the plugin: a container outliving its plugin comes back empty.
        results[task] = (plugin, harness.walk(plugin, task))
    return results, unloadable


def test_the_registry_is_walkable(walked):
    results, _unloadable = walked
    assert len(results) > 100, \
        f'only {len(results)} tasks could be built; the harness is not measuring much'


def test_every_route_to_a_parameter_agrees(walked):
    """dataOrder(), getattr and hasattr must not contradict each other."""
    results, _ = walked
    divergences = [d for _plugin, w in results.values() for d in w.divergences]
    assert not divergences, (
        f'{len(divergences)} disagreements between access mechanisms:\n  '
        + '\n  '.join(str(d) for d in divergences[:20]))


def test_the_walk_reaches_a_realistic_number_of_parameters(walked):
    """A guard on the guard: a harness that reaches nothing passes everything."""
    results, _ = walked
    total = sum(len(w.paths) for _plugin, w in results.values())
    assert total > 15000, f'only {total} paths reached'


def test_no_task_became_unloadable(walked):
    """Not a CData property, but this is where it would first show."""
    _results, unloadable = walked
    unexpected = [t for t, _why in unloadable
                  if t not in {'TestObsConversions', 'pisa', 'AMPLE'}]
    assert not unexpected, \
        'tasks that no longer load: ' + ', '.join(f'{t} ({w})' for t, w in unloadable
                                                  if t in unexpected)
