#!/usr/bin/env ccp4-python
"""Snapshot every task's container tree, for before/after comparison.

The instrument the container-construction work is measured with. It produced
the numbers in docs/container-construction-defects.md --- "161 trees
byte-identical, 0 paths gained, 1,358 removed" --- and it lives here rather
than in someone's scratch directory so those numbers can be reproduced by
anyone, and so the next structural change to CData can be held to the same
standard.

Records what the container *is*, and what the rest of the system *derives*
from it. The second half matters more for a structural change: tiers 1 to 5 of
the conformance suite assert that the access mechanisms agree with each other,
which is a different claim from "nothing any consumer sees has changed".

  paths     the sorted set of full paths in each task's container. A change
            that removes ghosts *removes* paths; the acceptance test for a
            removal is that nothing is ever gained.
  order     the top-level child order from several fresh instantiations in one
            process. Non-determinism shows up as these differing.
  keywords  every task's i2run addressing table --- path, minimumPath, and
            which candidate a bare --FLAG resolves to. Derived entirely from
            container structure, so this is the cheap proxy for what the i2run
            suite checks about addressing, without running a job. It would have
            caught the sibling-ordering defect on its own.
  json      the shape the GUI renders, from CCP4i2JsonEncoder, whose ordering
            comes from dataOrder(). If membership or order shifts, every task
            interface shifts.
  validity  the error report validity() produces, which is what the Run dialog
            shows. It recurses children.

Usage
-----
    cd server
    ccp4-python ../scripts/snapshot_containers.py before.json
    # ... make the change ...
    ccp4-python ../scripts/snapshot_containers.py after.json
    ccp4-python ../scripts/snapshot_containers.py --diff before.json after.json

The diff is the thing to read: it reports trees changed, paths gained (which
should be zero for a removal), paths lost, any task that stopped loading, and
whether the i2run addressing, GUI-rendered shape or validity report differ.

**Every difference must be predicted before it is observed.** An unpredicted
difference is a defect; a predicted one is the deliverable. Defect A's
acceptance test was not "nothing changed" --- it was "0 paths gained, 1,358
removed, in these ten tasks", a shape stated in advance. So a gain fails the
run outright, since nothing should invent a parameter; losses and derived
changes are reported for a human to check against their prediction; and
``--strict`` fails on any difference at all, which is what a pure refactor
should assert.

Not every change is visible here. Content cached on *read* --- a file's cell
and spacegroup --- never appears, because these containers are freshly built
and have read nothing. Match the instrument to the change: that one belongs to
the persistence tier and to a before/after over real ``params.xml`` files.
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


def keyword_table(task):
    """Every parameter's i2run addressing, derived from container structure."""
    from ccp4i2.cli.i2run.i2run_components import KeywordExtractor

    rows = []
    for keyword in KeywordExtractor.extract_from_task_name(task):
        rows.append({
            'path': str(keyword.get('path', '')),
            'simpleName': str(keyword.get('simpleName', '')),
            'minimumPath': str(keyword.get('minimumPath', '')),
            'shortest': bool(keyword.get('isShortestForSimpleName', False)),
            'ambiguous': bool(keyword.get('isAmbiguousSimpleName', False)),
        })
    return sorted(rows, key=lambda r: (r['path'], r['simpleName']))


def json_shape(container):
    """What the GUI renders: class and child order at every node.

    Values are deliberately excluded --- they carry temp directories and other
    per-run noise. The shape is what a structural change would disturb.
    """
    import json as _json

    from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder

    encoded = _json.loads(_json.dumps(container, cls=CCP4i2JsonEncoder))

    def shape(node, depth=0):
        if depth > MAX_DEPTH or not isinstance(node, dict):
            return None
        out = {'_class': node.get('_class')}
        value = node.get('_value')
        if isinstance(value, dict):
            out['_children'] = [[k, shape(v, depth + 1)] for k, v in value.items()]
        elif isinstance(value, list):
            out['_items'] = [shape(v, depth + 1) for v in value]
        return out

    return shape(encoded)


def validity_report(plugin):
    """The errors the Run dialog would show, as codes and names."""
    try:
        error = plugin.validity()
    except Exception as err:
        return [{'error': f'{type(err).__name__}: {err}'}]
    rows = []
    for entry in error.entries() if hasattr(error, 'entries') else []:
        rows.append({
            'class': str(entry.get('class', '')),
            'code': entry.get('code'),
            'name': str(entry.get('name', '')),
            'severity': entry.get('severity'),
        })
    return sorted(rows, key=lambda r: (r['name'], str(r['code'])))


def snapshot_task(task):
    from ccp4i2.core.tasks import get_plugin_class

    plugin_class = get_plugin_class(task)
    if plugin_class is None:
        return {'error': 'did not import'}
    orders, paths, shape, validity = [], None, None, None
    for _ in range(INSTANTIATIONS):
        # Hold the plugin: a container outliving its plugin comes back empty,
        # because children are reached through weak references.
        plugin = plugin_class(workDirectory=tempfile.mkdtemp(), name='snapshot')
        orders.append(list(plugin.container.dataOrder()))
        if paths is None:
            collected = []
            walk(plugin.container, '', collected)
            paths = sorted(set(collected))
            shape = json_shape(plugin.container)
            validity = validity_report(plugin)

    out = {'paths': paths, 'order': orders, 'json': shape, 'validity': validity}
    try:
        out['keywords'] = keyword_table(task)
    except Exception as err:
        out['keywords'] = [{'error': f'{type(err).__name__}: {err}'}]
    return out


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
    print(f'keywords:          {sum(len(out[t].get("keywords", [])) for t in good)}')
    print(f'validity entries:  {sum(len(out[t].get("validity", [])) for t in good)}')
    print(f'unstable order:    {len(unstable)} tasks')
    for task in unstable[:10]:
        print(f'  {task}')
    return out


def diff(before_path, after_path, strict=False):
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

    # What the rest of the system derives from the container. Any difference
    # here is a difference a user or a script can see.
    derived_changed = 0
    for field, label in (('keywords', 'i2run addressing'),
                         ('json', 'GUI-rendered shape'),
                         ('validity', 'validity report')):
        differing = [t for t in sorted(before)
                     if field in before[t] and field in after.get(t, {})
                     and before[t][field] != after[t][field]]
        derived_changed += len(differing)
        print(f'{label:20}: {len(differing)} tasks differ'
              + (f'  {differing[:6]}' if differing else ''))

    # Gains and losses are not the same kind of event. Inventing a parameter
    # that was not there is almost always wrong; removing one is frequently the
    # entire point --- defect A removed 1,358 ghost paths on purpose. So a gain
    # is fatal by default and everything else is reported for a human to check
    # against what they predicted. --strict fails on any difference at all,
    # which is what a pure refactor should assert.
    if gained:
        print('\nFAIL: paths were gained. Nothing should invent a parameter.')
        return 1
    if strict and (lost or derived_changed):
        print('\nFAIL: --strict, and something differs. For a refactor that '
              'should be nothing; for a fix, drop --strict and check the diff '
              'against what you predicted.')
        return 1
    if lost or derived_changed:
        print('\nDifferences above are not failures by themselves. They are '
              'the deliverable if you predicted them, and a defect if you did '
              'not.')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+',
                        help='output file, or two snapshots with --diff')
    parser.add_argument('--diff', action='store_true',
                        help='compare two snapshots instead of taking one')
    parser.add_argument('--strict', action='store_true',
                        help='fail on any difference, not only on gained paths. '
                             'What a pure refactor should assert.')
    args = parser.parse_args()

    if args.diff:
        if len(args.paths) != 2:
            parser.error('--diff needs exactly two snapshot files')
        return diff(*args.paths, strict=args.strict)

    _setup_django()
    take(args.paths[0])
    return 0


if __name__ == '__main__':
    sys.exit(main())
