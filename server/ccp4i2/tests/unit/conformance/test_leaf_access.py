"""
What plugins do to a parameter once they have reached it.

Tier 2. Tier 1 asked whether the routes to a parameter agree; this asks
whether the operations plugins perform on it hold, across every leaf in every
registered task rather than the handful anyone tests by hand.

The operations are the measured ones: `str()` (526 explicit conversions across
87 files, together with int() and float()), `isSet()` (719 in 110), and `==`
(102 in 28). These are load-bearing --- `isSet()` decides whether a parameter
reaches a program at all.

Known protocol gaps are *pinned* in test_protocol_pins.py rather than asserted
here, because fixing them is part of the reevaluation, not a precondition of it.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.tests.unit.conformance import harness


@pytest.fixture(scope='module')
def leaves():
    found = []
    for task in harness.declared_tasks():
        try:
            plugin = harness.build(task)
        except Exception:
            continue
        if plugin is None:
            continue
        found.append((task, plugin, harness.walk(plugin, task).leaves))
    return found


def test_there_are_leaves_to_check(leaves):
    total = sum(len(l) for _t, _p, l in leaves)
    assert total > 8000, f'only {total} leaves reached; the tier is not measuring much'


def test_str_never_raises(leaves):
    """Plugins interpolate parameters into command lines constantly."""
    broken = []
    for task, _plugin, found in leaves:
        for path, obj in found:
            try:
                str(obj)
            except Exception as err:
                broken.append(f'{task}: {path}: {type(err).__name__}: {err}')
    assert not broken, f'{len(broken)} parameters cannot be stringified:\n  ' \
        + '\n  '.join(broken[:15])


def test_is_set_never_raises(leaves):
    """719 uses decide whether a parameter reaches a program."""
    broken = []
    for task, _plugin, found in leaves:
        for path, obj in found:
            checker = getattr(obj, 'isSet', None)
            if not callable(checker):
                continue
            try:
                checker()
            except Exception as err:
                broken.append(f'{task}: {path}: {type(err).__name__}: {err}')
    assert not broken, f'{len(broken)} parameters cannot report whether they are set:\n  ' \
        + '\n  '.join(broken[:15])


def test_comparison_with_a_plain_value_never_raises(leaves):
    broken = []
    for task, _plugin, found in leaves:
        for path, obj in found:
            try:
                obj == 'a string that is not the value'
            except Exception as err:
                broken.append(f'{task}: {path}: {type(err).__name__}: {err}')
    assert not broken, f'{len(broken)} parameters cannot be compared:\n  ' \
        + '\n  '.join(broken[:15])


def test_declared_conversions_work(leaves):
    """Where a type says it converts, converting must not raise.

    Not "every numeric type converts" --- that is a gap to be fixed, and it is
    pinned elsewhere. This is the weaker, non-negotiable property: a type that
    advertises __int__ or __float__ must honour it.
    """
    broken = []
    for task, _plugin, found in leaves:
        for path, obj in found:
            for dunder, fn in (('__int__', int), ('__float__', float)):
                if not hasattr(type(obj), dunder):
                    continue
                try:
                    fn(obj)
                except (TypeError, ValueError):
                    # An unset parameter has nothing to convert; that is not a
                    # broken protocol.
                    if getattr(obj, 'isSet', lambda: True)():
                        broken.append(f'{task}: {path}: {type(obj).__name__}.{dunder}')
                except Exception as err:
                    broken.append(f'{task}: {path}: {type(err).__name__}: {err}')
    assert not broken, f'{len(broken)} declared conversions raise:\n  ' \
        + '\n  '.join(broken[:15])
