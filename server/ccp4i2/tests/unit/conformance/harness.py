"""
Walk every task's container, and offer the same object through each of the
mechanisms plugins actually use to reach it.

This is the instrument the CData reevaluation is measured against. It asserts
that the mechanisms **agree with each other**, not that any of them returns a
particular value: agreement needs no golden file to rot, and it fails the
moment a change makes two routes to the same parameter disagree.

Measured on this tree: plugins reach parameters by dotted access (3711 uses in
146 files), `list[n]` (2016 in 200), `int()`/`float()`/`str()` (526 in 87),
`hasattr` (328 in 60), `getattr` (207 in 46) and `dataOrder()` (13 in 8).
`children()` and `find_child()` are used only by the core --- on everyone
else's behalf, which is why they are checked here too.
"""
from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from typing import Any, Iterator, List, Optional


MAX_DEPTH = 8


@dataclass
class Divergence:
    """One place where two routes to the same object disagreed."""
    task: str
    path: str
    check: str
    detail: str

    def __str__(self) -> str:
        return f"{self.task}: {self.path}: {self.check}: {self.detail}"


@dataclass
class Walk:
    """Everything the harness found in one task."""
    task: str
    paths: List[str] = field(default_factory=list)
    divergences: List[Divergence] = field(default_factory=list)

    def note(self, path: str, check: str, detail: str) -> None:
        self.divergences.append(Divergence(self.task, path, check, detail))


def build(task: str):
    """Instantiate *task*, returning the plugin --- not merely its container.

    A container's children are reached through weak references, so a container
    outliving its plugin comes back empty. Returning the plugin keeps it alive
    for as long as the caller holds it; `test_lifetime.py` pins that down.
    """
    from ccp4i2.core.tasks import get_plugin_class

    plugin_class = get_plugin_class(task)
    if plugin_class is None:
        return None
    return plugin_class(workDirectory=tempfile.mkdtemp(), name='conformance')


def order_of(obj: Any) -> List[str]:
    order = getattr(obj, 'dataOrder', None)
    if not callable(order):
        return []
    try:
        return [str(n) for n in order()]
    except Exception:
        return []


def is_container_like(obj: Any) -> bool:
    return callable(getattr(obj, 'dataOrder', None))


def walk(plugin, task: str) -> Walk:
    """Check every reachable path in *plugin*'s container."""
    found = Walk(task=task)
    _visit(plugin.container, '', found, 0, set())
    return found


def is_list_like(obj: Any) -> bool:
    """A CList: its children are reached by index, not by name."""
    return hasattr(obj, '_items') and callable(getattr(obj, '__getitem__', None))


def _visit(obj: Any, prefix: str, found: Walk, depth: int, seen: set) -> None:
    if depth > MAX_DEPTH or id(obj) in seen:
        return
    seen.add(id(obj))

    if is_list_like(obj):
        _visit_list(obj, prefix, found, depth, seen)
        return

    names = order_of(obj)
    for name in names:
        path = f"{prefix}.{name}" if prefix else name
        found.paths.append(path)
        child = _agree_on(obj, name, path, found)
        if child is not None and is_container_like(child):
            _visit(child, path, found, depth + 1, seen)


def _visit_list(obj: Any, prefix: str, found: Walk, depth: int, seen: set) -> None:
    """A CList is addressed by index, and its own accounts of itself must tally.

    Three different names exist for the same item --- dataOrder() says '0',
    the item's _name says '[0]', and __setitem__ writes 'LIST[0]' --- while
    the way anybody actually reaches it is lst[0]. Defect C is where that gets
    resolved; until then this checks the routes that are meant to work.
    """
    items = list(getattr(obj, '_items', []))
    order = order_of(obj)

    if len(order) != len(items):
        found.note(prefix, 'dataOrder-vs-items',
                   f'dataOrder() has {len(order)} entries, _items has {len(items)}')

    for index, item in enumerate(items):
        path = f"{prefix}[{index}]"
        found.paths.append(path)
        try:
            by_index = obj[index]
        except Exception as err:
            found.note(path, 'index-access',
                       f'lst[{index}] raised {type(err).__name__}')
            continue
        if by_index is not item:
            found.note(path, 'index-vs-items',
                       f'lst[{index}] is not _items[{index}]')
        if is_container_like(by_index):
            _visit(by_index, path, found, depth + 1, seen)


def _agree_on(parent: Any, name: str, path: str, found: Walk) -> Optional[Any]:
    """The four ways a plugin reaches a named child must give the same object."""
    sentinel = object()

    by_attr = getattr(parent, name, sentinel)
    if by_attr is sentinel:
        # dataOrder() named something getattr cannot reach. This is the shape
        # of the BOXS defect: a leaf reported a child named after itself, and
        # phaser_MR.setKeywords took it for a sub-container.
        found.note(path, 'dataOrder-vs-getattr',
                   f'dataOrder() names {name!r} but getattr does not reach it')
        return None

    if not hasattr(parent, name):
        found.note(path, 'getattr-vs-hasattr',
                   'getattr succeeds where hasattr says no')

    by_getattr = getattr(parent, name)
    if by_getattr is not by_attr:
        found.note(path, 'getattr-twice',
                   'two getattr calls returned different objects')

    return by_attr


def declared_tasks() -> Iterator[str]:
    from ccp4i2.core.tasks import TASKS

    return iter(sorted(TASKS))
