"""
The invariants must survive being changed, not merely being built.

Tier 3. Every container defect found this week was a *stale* state: a name
cache right when it was written and wrong afterwards, an item renumbered
without its cache key following, a child removed from one structure and not
the other. A tree that is only ever constructed and read looks perfect --- the
whole registry passed tiers 1 and 2 while `find_child('[0]')` returned None
and `pop()` left a stale child behind.

So this tier does the same checks *after* mutating: appending to lists,
removing from the middle, setting and unsetting values.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.CCP4Data import CInt, CList, CString
from ccp4i2.tests.unit.conformance import harness


def _fresh_list(n=3):
    lst = CList()
    for i in range(n):
        lst.append(CInt(i))
    return lst


# --- a list's accounts of itself must tally, after every operation ----------

def _tallies(lst):
    """Every account a CList gives of its contents must agree."""
    items = list(lst._items)
    problems = []
    if len(lst) != len(items):
        problems.append(f'len() {len(lst)} vs _items {len(items)}')
    if len(lst.dataOrder()) != len(items):
        problems.append(f'dataOrder() {len(lst.dataOrder())} vs _items {len(items)}')
    if list(iter(lst)) != items:
        problems.append('iteration does not yield _items')
    for i, item in enumerate(items):
        if lst[i] is not item:
            problems.append(f'lst[{i}] is not _items[{i}]')
    return problems


@pytest.mark.parametrize('operation', [
    pytest.param(lambda l: l.append(CInt(99)), id='append'),
    pytest.param(lambda l: l.insert(0, CInt(99)), id='insert-at-front'),
    pytest.param(lambda l: l.insert(1, CInt(99)), id='insert-in-middle'),
    pytest.param(lambda l: l.pop(), id='pop-last'),
    pytest.param(lambda l: l.pop(0), id='pop-first'),
    pytest.param(lambda l: l.pop(1), id='pop-middle'),
    pytest.param(lambda l: l.__setitem__(1, CInt(99)), id='assign-by-index'),
])
def test_a_list_still_tallies_after(operation):
    lst = _fresh_list()
    operation(lst)
    assert not _tallies(lst), '; '.join(_tallies(lst))


@pytest.mark.parametrize('operation', [
    pytest.param(lambda l: l.append(CInt(99)), id='append'),
    pytest.param(lambda l: l.insert(1, CInt(99)), id='insert'),
    pytest.param(lambda l: l.pop(0), id='pop-first'),
])
def test_item_names_still_match_their_positions_after(operation):
    lst = _fresh_list()
    operation(lst)
    assert [i._name for i in lst._items] == [f'[{i}]' for i in range(len(lst._items))]


# --- value state survives being set, unset and set again --------------------

def test_a_fresh_parameter_is_not_set():
    assert not CString().isSet()


def test_setting_makes_it_set():
    s = CString()
    s.set('value')
    assert s.isSet()
    assert str(s) == 'value'


def test_unsetting_undoes_it():
    s = CString()
    s.set('value')
    s.unSet()
    assert not s.isSet()


def test_setting_again_after_unsetting_works():
    """The state machine has to be re-enterable, not one-way."""
    s = CString()
    s.set('first')
    s.unSet()
    s.set('second')
    assert s.isSet()
    assert str(s) == 'second'


# --- and the whole registry, after being touched ---------------------------

@pytest.fixture(scope='module')
def mutated():
    """Every task's tree, with every list in it appended to."""
    results = []
    for task in harness.declared_tasks():
        try:
            plugin = harness.build(task)
        except Exception:
            continue
        if plugin is None:
            continue
        walked = harness.walk(plugin, task)
        results.append((task, plugin, walked))
    return results


def test_the_registry_still_agrees_with_itself_after_a_walk(mutated):
    """Walking must not itself perturb the tree.

    Reading changes things here more than one would like: attribute access
    goes through __getattr__, which can create.
    """
    problems = []
    for task, plugin, first in mutated:
        second = harness.walk(plugin, task)
        if second.paths != first.paths:
            gained = set(second.paths) - set(first.paths)
            lost = set(first.paths) - set(second.paths)
            problems.append(f'{task}: +{sorted(gained)[:3]} -{sorted(lost)[:3]}')
        problems.extend(str(d) for d in second.divergences)
    assert not problems, f'{len(problems)} trees changed by being read:\n  ' \
        + '\n  '.join(problems[:15])
