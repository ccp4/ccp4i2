"""Which tasks' help is current: the status page, regenerated on every build.

For each task the task chooser offers (client task-chooser.tsx, in its
categories) the page is one of:

  current   its figures are captured from this app (a shots.json beside it)
  draft     the same, but written new and not yet reviewed by someone who
            knows the program ("draft": true in its shots.json)
  Qt        it still shows the Qt interface
  none      there is no page

A current or draft page is also marked *stale* when the files it describes
have changed since it was checked against them (tools/stamp.py).

conf.py calls write_status(); run directly to print the counts, or with
`todo` to list the tasks still Qt or without a page, by chooser category
(each task once, under its first category): the menu routes are chosen from.
"""
import json
import re
from functools import cache
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent
REPO = DOCS.parent.parent
CHOOSER = REPO / "client/renderer/components/task/task-chooser.tsx"
TASKS_PY = REPO / "server/ccp4i2/core/tasks.py"
TASK_PAGES = DOCS / "source/tasks"

# Pages whose directory is not the task's name.
ALIASES = {
    "AMPLE": "ample",
    "i2Dimple": "dimple",
    "import_serial_pipe": "import_serial",
    "add_fractional_coords": "fractional_coordinates",
    "ProvideTLS": "tls",
    "sheetbend": "shift_field",
    "qtpisa": "pisapipe",
    "phaser_ensembler": "ensemble_phaser",
    "molrep_selfrot": "srf",
    "coot_rsr_morph": "coot_refinement",
    "dials_image": "dials",
    "phaser_pipeline_phil": "phaser_pipeline",
    "adding_stats_to_mmcif_i2": "PrepareDeposit",
    # One page for Phaser MR over PHIL, whole or in steps.
    "phaser_mr_auto_phil": "phaser_mr_phil",
    "phaser_mr_frf_phil": "phaser_mr_phil",
    "phaser_mr_ftf_phil": "phaser_mr_phil",
    "phaser_mr_pak_phil": "phaser_mr_phil",
    "phaser_mr_rnp_phil": "phaser_mr_phil",
    # One page for the single-file import tasks.
    "ImportCoordinate": "import_files",
    "ImportSequence": "import_files",
    "ImportAsuContent": "import_files",
    "ImportDictionary": "import_files",
    "ImportMap": "import_files",
    "ImportObs": "import_files",
    "ImportUnmerged": "import_files",
    "ImportMapCoeffs": "import_files",
    "ImportFreeR": "import_files",
    "ImportPhases": "import_files",
}


@cache
def superseded():
    """Tasks with a successor: the chooser hides them, so they need no page."""
    text = TASKS_PY.read_text(encoding="utf-8")
    return set(re.findall(
        r'"(\w+)":\s*Task\((?:(?!\n    "\w+":\s*Task\()[\s\S])*?successor="\w+"', text))


def chooser_categories():
    """[(category title, [task names])], as the task chooser shows them."""
    text = CHOOSER.read_text(encoding="utf-8")
    block = text[text.index("const TASK_CATEGORIES"):]
    block = block[:block.index("\n];")]
    hidden = superseded()
    return [(title, [t for t in re.findall(r'"([^"]+)"', tasks) if t not in hidden])
            for title, tasks in re.findall(
                r'title:\s*"([^"]+)"[\s\S]*?tasks:\s*\[([^\]]*)\]', block)]


@cache
def task_documents():
    """{page directory: its document}, from the tasks index's toctrees
    (a page is not always index.rst: crank2/crank2, dials/image_viewer)."""
    docs = {}
    index = (TASK_PAGES / "index.rst").read_text(encoding="utf-8")
    for entry in re.findall(r"^   (\S+/\S+)\s*$", index, flags=re.M):
        docs.setdefault(entry.split("/")[0], entry)
    return docs


def state(task):
    directory = ALIASES.get(task, task)
    doc = task_documents().get(directory)
    if doc is None or not (TASK_PAGES / f"{doc}.rst").exists():
        return "none", None
    shots = TASK_PAGES / directory / "shots.json"
    if not shots.exists():
        return "Qt", doc
    draft = json.loads(shots.read_text(encoding="utf-8")).get("draft", False)
    return ("draft" if draft else "current"), doc


def is_stale(task):
    from stamp import staleness  # stamp imports this module
    return bool(staleness(ALIASES.get(task, task)))


def write_status(out: Path):
    counts = {"current": 0, "draft": 0, "Qt": 0, "none": 0, "stale": 0}
    seen = set()
    rows = []
    for title, tasks in chooser_categories():
        rows.append(f"\n{title}\n{'-' * len(title)}\n\n"
                    ".. list-table::\n   :widths: 40 20\n\n")
        for task in tasks:
            s, page = state(task)
            stale = s in ("current", "draft") and is_stale(task)
            if task not in seen:
                counts[s] += 1
                counts["stale"] += stale
                seen.add(task)
            name = f":doc:`{task} <tasks/{page}>`" if page else task
            shown = f"{s}, **stale**" if stale else s
            rows.append(f"   * - {name}\n     - {shown}\n")
    total = sum(v for k, v in counts.items() if k != "stale")
    head = (
        "####################\n"
        "Documentation status\n"
        "####################\n\n"
        "Generated by ``tools/status.py`` on every build, from the tasks the\n"
        "task chooser offers. *current*: pictures captured from this app;\n"
        "*draft*: the same, newly written and awaiting review by someone who\n"
        "knows the program; *Qt*: the page still shows the Qt interface;\n"
        "*none*: no page yet. *stale*: the task's interface, def.xml, script,\n"
        "report or scenario has changed since the page was checked against\n"
        "them, so it may no longer be right (``tools/stamp.py``).\n\n"
        f"Of {total} tasks: **{counts['current']} current**, "
        f"{counts['draft']} draft, {counts['Qt']} Qt, "
        f"{counts['none']} with no page"
        + (f"; {counts['stale']} stale" if counts["stale"] else "") + ".\n")
    out.write_text(head + "".join(rows), encoding="utf-8")
    return counts


def todo():
    """[(category, [(task, state)])] for tasks still Qt or without a page."""
    seen, out = set(), []
    for title, tasks in chooser_categories():
        left = [(t, state(t)[0]) for t in tasks
                if t not in seen and state(t)[0] in ("Qt", "none")]
        seen.update(tasks)
        if left:
            out.append((title, left))
    return out


if __name__ == "__main__":
    import sys
    import tempfile
    if sys.argv[1:] == ["todo"]:
        for title, left in todo():
            print(f"{title}:\n    " + ", ".join(f"{t} [{s}]" for t, s in left))
    else:
        with tempfile.TemporaryDirectory() as d:
            print(write_status(Path(d) / "status.rst"))
