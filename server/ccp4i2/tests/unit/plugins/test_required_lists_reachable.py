"""A list a task requires (listMinLength >= 1) must be one its interface shows.

Validity is the server's: a required list the interface never offers fails
every job, and in the app disables Confirm with no field to fix. SliceNDice
(ENSEMBLES), xia2_multiplex and xia2_ssx_reduce (XIA2_RUN, two required)
could not be run at all, each list a leftover the wrapper never reads;
mrbump_basic had the same list and a validity() override hiding its error.

A task without an interface of its own gets the generic one, which renders
every parameter, so only tasks with their own interface are checked.
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from ccp4i2.core.tasks import TASKS

PACKAGE = Path(__file__).resolve().parents[3]          # server/ccp4i2
INTERFACES = PACKAGE.parents[1] / "client/renderer/components/task/task-interfaces"
CONTAINER = INTERFACES / "task-container.tsx"

pytestmark = pytest.mark.skipif(not CONTAINER.is_file(), reason="needs the client tree")


def interface_text(task):
    """The task's own interface and the local files it imports, or None."""
    text = CONTAINER.read_text(encoding="utf-8")
    m = re.search(rf'^\s*"?{re.escape(task)}"?:\s*(\w+),', text, flags=re.M)
    if not m:
        return None
    imp = re.search(rf'import\s+{m.group(1)}\s+from\s+"\./([^"]+)"', text)
    if not imp:
        return None
    found, todo = [], [INTERFACES / imp.group(1)]
    while todo:
        base = todo.pop()
        f = next((p for p in (base, base.with_suffix(".tsx"), base.with_suffix(".ts"),
                              base / "index.tsx") if p.is_file()), None)
        # Not the container: it imports every task's interface.
        if f is None or f in found or INTERFACES not in f.parents or f == CONTAINER:
            continue
        found.append(f)
        for spec in re.findall(r'from\s+"(\.\.?/[^"]+)"', f.read_text(encoding="utf-8")):
            todo.append((f.parent / spec).resolve())
    return "".join(f.read_text(encoding="utf-8") for f in found)


def required_lists(def_xml):
    for content in ET.parse(def_xml).getroot().iter("content"):
        minimum = content.findtext("qualifiers/listMinLength", "").strip()
        if minimum.isdigit() and int(minimum) >= 1:
            yield content.get("id")


def cases():
    for task, entry in TASKS.items():
        if not entry.defXmlPath:
            continue
        def_xml = PACKAGE / entry.defXmlPath
        if def_xml.is_file():
            for name in required_lists(def_xml):
                yield task, name, def_xml


@pytest.mark.parametrize("task,name,def_xml", list(cases()))
def test_required_list_is_in_the_interface(task, name, def_xml):
    text = interface_text(task)
    if text is None:
        pytest.skip("generic interface: every parameter is shown")
    assert re.search(rf'\b{re.escape(name)}\b', text), (
        f"{task} requires {name} (listMinLength in {def_xml.name}) but its "
        "interface never shows it: every job fails validity")
