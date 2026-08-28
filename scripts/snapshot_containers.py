#!/usr/bin/env ccp4-python
"""Snapshot every task's container tree, for before/after comparison.

The instrument the container-construction work is measured with. It produced
the numbers in docs/container-construction-defects.md --- "161 trees
byte-identical, 0 paths gained, 1,358 removed" --- and it lives here rather
than in someone's scratch directory so those numbers can be reproduced by
anyone, and so the next structural change to CData can be held to the same
standard.

Records two things, because the two questions differ:

- ``paths``: the sorted set of full paths in each task's container. A change
  that removes ghosts *removes* paths; the acceptance test is that nothing is
  ever gained.
- ``order``: the top-level child order from several fresh instantiations of
  the same task in one process. Non-determinism shows up as these differing
  from each other.

Usage
-----
    cd server
    ccp4-python ../scripts/snapshot_containers.py before.json
    # ... make the change ...
    ccp4-python ../scripts/snapshot_containers.py after.json
    ccp4-python ../scripts/snapshot_containers.py --diff before.json after.json

The diff is the thing to read: it reports trees changed, paths gained (which
should be zero for a removal), paths lost, and any task that stopped loading.
"""
import argparse
import json
import os
import sys
import tempfile

INSTANTIATIONS = 3
MAX_DEPTH = 12


def _setup_django():
    sys.path.insert(0, os.getcwd())
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4i2.config.test_settings')
    import django
    django.setup()


def walk(obj, prefix, out, depth=0):
    if depth > MAX_DEPTH:
        return
    order = getattr(obj, 'dataOrder', None)
    names = order() if callable(order) else []
    for name in names:
        child = getattr(obj, name, None)
        path = f"{prefix}.{name}" if prefix else name
        out.append(path)
        if child is not None and callable(getattr(child, 'dataOrder', None)):
            walk(child, path, out, depth + 1)


def snapshot_task(task):
    from ccp4i2.core.tasks import get_plugin_class

    plugin_class = get_plugin_class(task)
    if plugin_class is None:
        return {'error': 'did not import'}
    orders, paths = [], None
    for _ in range(INSTANTIATIONS):
        # Hold the plugin: a container outliving its plugin comes back empty,
        # because children are reached through weak references.
        plugin = plugin_class(workDirectory=tempfile.mkdtemp(), name='snapshot')
        orders.append(list(plugin.container.dataOrder()))
        if paths is None:
            collected = []
            walk(plugin.container, '', collected)
            paths = sorted(set(collected))
    return {'paths': paths, 'order': orders}


def take(out_path):
    from ccp4i2.core.tasks import TASKS

    out = {}
    for task in sorted(TASKS):
        try:
            out[task] = snapshot_task(task)
        except Exception as err:
            out[task] = {'error': f'{type(err).__name__}: {err}'}
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True)

    good = [t for t, v in out.items() if 'paths' in v]
    unstable = [t for t in good
                if len({tuple(o) for o in out[t]['order']}) > 1]
    print(f'tasks snapshotted: {len(good)} of {len(out)}')
    print(f'total paths:       {sum(len(out[t]["paths"]) for t in good)}')
    print(f'unstable order:    {len(unstable)} tasks')
    for task in unstable[:10]:
        print(f'  {task}')
    return out


def diff(before_path, after_path):
    with open(before_path) as f:
        before = json.load(f)
    with open(after_path) as f:
        after = json.load(f)

    same = gained = lost = 0
    changed = []
    for task in sorted(before):
        if 'paths' not in before[task] or 'paths' not in after.get(task, {}):
            continue
        b, a = set(before[task]['paths']), set(after[task]['paths'])
        g, l = a - b, b - a
        if not g and not l:
            same += 1
        else:
            changed.append((task, len(g), len(l)))
            gained += len(g)
            lost += len(l)

    was_ok = {t for t, v in before.items() if 'paths' in v}
    now_ok = {t for t, v in after.items() if 'paths' in v}

    print(f'trees byte-identical : {same}')
    print(f'trees that change    : {len(changed)}')
    print(f'paths GAINED anywhere: {gained}')
    print(f'paths removed        : {lost}')
    for task, g, l in sorted(changed, key=lambda x: -x[2]):
        print(f'   {task:26} +{g} -{l}')
    if was_ok - now_ok:
        print(f'STOPPED LOADING      : {sorted(was_ok - now_ok)}')
    if now_ok - was_ok:
        print(f'started loading      : {sorted(now_ok - was_ok)}')

    for label, snap in (('before', before), ('after', after)):
        bad = [t for t, v in snap.items()
               if 'order' in v and len({tuple(o) for o in v['order']}) > 1]
        print(f'unstable order {label:6}: {len(bad)} tasks')

    return 1 if gained else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+',
                        help='output file, or two snapshots with --diff')
    parser.add_argument('--diff', action='store_true',
                        help='compare two snapshots instead of taking one')
    args = parser.parse_args()

    if args.diff:
        if len(args.paths) != 2:
            parser.error('--diff needs exactly two snapshot files')
        return diff(*args.paths)

    _setup_django()
    take(args.paths[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())
