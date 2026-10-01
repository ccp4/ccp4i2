"""What a help page needs to know about a task, in one place.

    python3 docs/user/tools/taskcard.py <task> [<task>...]

Prints, from the source alone (no server, no CCP4):

- the registry entry: title, plugin, def.xml, report, flags;
- where the chooser offers it, and the state of its page;
- the def.xml parameters by container: class, label, default, menu,
  whether it may be left unset;
- the interface: its file(s), tabs, sections, and every element with the
  label the interface gives it (its own guiLabel, or the def.xml's);
- the i2run tests that run it.

It is a reading aid, built by pattern rather than by running the interface:
an element the interface builds in a loop, or labels computed at run time
(a PHIL task's parameters), will not appear. The capture's outline mode
(capture.mjs, "outline": true) shows the page as it renders.
"""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import status
from stamp import SERVER, interface_files, task_entries

REPO = Path(__file__).resolve().parents[3]
TESTS = SERVER / "tests" / "i2run"


def registry(task):
    text = status.TASKS_PY.read_text(encoding="utf-8")
    m = re.search(rf'^    "{re.escape(task)}":\s*Task\(([\s\S]*?)^    \),', text, flags=re.M)
    if not m:
        return {}
    return dict(re.findall(r'(\w+)=("[^"]*"|True|False)', m.group(1)))


def parameters(def_xml: Path):
    """[(container, id, class, label, default, menu, allowUndefined)]"""
    rows = []
    root = ET.parse(def_xml).getroot()

    def walk(node, where):
        for child in node:
            if child.tag == "container":
                walk(child, where + [child.get("id") or "?"])
            elif child.tag == "content":
                q = child.find("qualifiers")
                get = (lambda k: (q.findtext(k) or "").strip()) if q is not None else (lambda k: "")
                menu = get("menuText") or get("enumerators")
                rows.append(("/".join(where[1:]) or where[0], child.get("id"),
                             child.findtext("className") or "", get("guiLabel"),
                             get("default"), menu, get("allowUndefined")))
    walk(root.find(".//ccp4i2_body") if root.find(".//ccp4i2_body") is not None else root, ["body"])
    return rows


# An element: itemName="X" ... optional qualifiers={{ ... guiLabel: "Y" ... }}
ELEMENT = re.compile(
    r'<CCP4i2(TaskElement|ContainerElement)\b(?P<attrs>(?:[^<>]|=>)*?)/?>', re.S)


def interface(tsx: Path):
    """[(kind, text)] in source order: tabs, sections, elements."""
    text = tsx.read_text(encoding="utf-8")
    found = []
    for m in re.finditer(r'<CCP4i2Tab\b[^>]*?label="([^"]+)"', text):
        found.append((m.start(), "tab", m.group(1)))
    for m in ELEMENT.finditer(text):
        attrs = m.group("attrs")
        item = re.search(r'itemName="([^"]*)"', attrs)
        label = re.search(r'guiLabel:\s*"([^"]*)"', attrs)
        mode = re.search(r'guiMode:\s*"([^"]*)"', attrs)
        name = item.group(1) if item else "?"
        if m.group(1) == "ContainerElement" and not name:
            found.append((m.start(), "section", label.group(1) if label else "(container)"))
        else:
            shown = label.group(1) if label else ""
            found.append((m.start(), "element",
                          f"{name}" + (f'  "{shown}"' if shown else "")
                          + (f"  [{mode.group(1)}]" if mode else "")))
    for m in re.finditer(r'<InlineField\b[^>]*?label="([^"]+)"', text):
        found.append((m.start(), "inline", m.group(1)))
    return [(kind, value) for _, kind, value in sorted(found)]


def tests(task):
    hits = []
    for f in sorted(TESTS.glob("test_*.py")):
        text = f.read_text(encoding="utf-8", errors="replace")
        if re.search(rf'["\']{re.escape(task)}["\']', text):
            n = len(re.findall(r"^def test_", text, flags=re.M))
            hits.append(f"{f.name} ({n} test{'s' if n != 1 else ''})")
    return hits


def card(task):
    out = [f"=== {task}"]
    reg = registry(task)
    if not reg:
        return out + ["  not in core/tasks.py"]
    for key in ("title", "shortTitle", "successor", "interactive", "ccp4_free"):
        if key in reg:
            out.append(f"  {key}: {reg[key].strip(chr(34))}")
    entry = task_entries().get(task, {})
    for key in ("pluginPath", "reportPath", "defXmlPath"):
        if key in entry:
            out.append(f"  {key}: {entry[key]}")

    cats = [title for title, tasks in status.chooser_categories() if task in tasks]
    state, page = status.state(task)
    directory = status.ALIASES.get(task, task)
    out.append(f"  chooser: {', '.join(cats) or 'not offered'}")
    out.append(f"  page: {state}" + (f" ({page})" if page else f" (directory {directory})")
               + (", STALE" if state in ("current", "draft") and status.is_stale(task) else ""))

    def_xml = SERVER / entry["defXmlPath"] if "defXmlPath" in entry else None
    if def_xml and def_xml.exists():
        out.append("  parameters (def.xml):")
        for where, pid, cls, label, default, menu, undef in parameters(def_xml):
            if where.startswith("outputData"):
                continue
            bits = [f"{where}.{pid}", cls]
            if label:
                bits.append(f'"{label}"')
            if default:
                bits.append(f"default={default[:30]}")
            if menu:
                bits.append(f"menu={menu[:70]}")
            if undef.lower() in ("true", "1"):
                bits.append("optional")
            out.append("    " + "  ".join(bits))
        outputs = [pid for where, pid, *_ in parameters(def_xml) if where.startswith("outputData")]
        if outputs:
            out.append("  outputs: " + ", ".join(outputs))

    files = interface_files(task)
    if files:
        for f in files[:1]:
            out.append(f"  interface: {f.relative_to(REPO)}")
            for kind, value in interface(f):
                out.append(f"    {kind:8} {value}")
        if len(files) > 1:
            out.append("  also: " + ", ".join(str(f.relative_to(REPO)) for f in files[1:]))
    else:
        out.append("  interface: none registered (generic interface from def.xml)")

    found = tests(task)
    out.append("  i2run tests: " + (", ".join(found) if found else "none"))
    return out


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for name in sys.argv[1:]:
        print("\n".join(card(name)))
