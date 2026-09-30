"""Is a help page still true to its task? Fingerprints of what it describes.

A page captured from the app describes the task as it was then: its interface
(the React task interface, the def.xml behind it), what it does (the wrapper
script) and what it reports (the report class). When any of those changes the
page may no longer be right, and nothing about the page says so. So each page's
shots.json records, under "sources", the git blob id of every file it was
checked against; a page whose files no longer match is *stale*.

The files are found, not listed by hand: the page's task (its directory, or
status.ALIASES, or "task" in shots.json), its def.xml, plugin and report modules
from server/ccp4i2/core/tasks.py, its interface from task-container.tsx plus the
files that interface imports from its own directory tree, and the scenario that
built the project the pictures were taken from.

It also checks what can be checked without a browser: that every label a shot
finds its way by (a section, a callout, from/until/through) still appears in
those sources. A renamed label breaks the next capture, and usually the prose.

Blob ids, not commit ids: a squash merge discards the branch's commits, but the
file contents it merged are in the repository, so `git cat-file -p <blob>`
recovers what the page was checked against, and `diff` shows what changed.

    stamp.py check [page ...]    list stale pages and why (--github: annotations)
    stamp.py diff page           what changed in its sources since it was stamped
    stamp.py stamp page ...      record the sources as they are now: after
                                 re-capturing, or after reading the diff and
                                 finding the page still right

Standard library only, so it runs in the docs CI job and from conf.py.
"""
import argparse
import difflib
import hashlib
import json
import re
import subprocess
import sys
from functools import cache
from pathlib import Path

from status import ALIASES, REPO, TASK_PAGES, TASKS_PY

INTERFACES = REPO / "client/renderer/components/task/task-interfaces"
CONTAINER = INTERFACES / "task-container.tsx"
ELEMENTS = REPO / "client/renderer/components/task/task-elements"
COMPONENTS = REPO / "client/renderer/components"
SERVER = REPO / "server/ccp4i2"


def blob_id(path: Path) -> str:
    """The id git gives this file's content (as `git hash-object` would,
    line endings normalised so a Windows checkout agrees)."""
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()


def rel(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


@cache
def task_entries():
    """{task: {pluginPath, defXmlPath, reportPath}} from tasks.py, by pattern
    (importing it would need the server's environment)."""
    text = TASKS_PY.read_text(encoding="utf-8")
    entries = {}
    for m in re.finditer(r'^    "(\w+)":\s*Task\(([\s\S]*?)^    \),', text, flags=re.M):
        entries[m.group(1)] = dict(re.findall(r'(\w+Path)="([^"]+)"', m.group(2)))
    return entries


def module_file(dotted: str) -> Path | None:
    """ccp4i2.pipelines.crank2.script.crank2_script:crank2 -> its .py file."""
    parts = dotted.split(":")[0].split(".")
    path = REPO / "server" / Path(*parts).with_suffix(".py")
    return path if path.exists() else None


def resolve_import(base: Path, spec: str) -> Path | None:
    target = (base.parent / spec).resolve()
    for candidate in (target.with_suffix(".tsx"), target.with_suffix(".ts"),
                      target / "index.tsx", target / "index.ts"):
        if candidate.is_file():
            return candidate
    return target if target.is_file() else None


def interface_files(task: str) -> list[Path]:
    """The task's registered interface and the components it imports (MakeLink's
    monomer editor), but not the shared task elements, which are every page's."""
    text = CONTAINER.read_text(encoding="utf-8")
    m = re.search(rf'^\s*"?{re.escape(task)}"?:\s*(\w+),', text, flags=re.M)
    if not m:
        return []
    imp = re.search(rf'import\s+{m.group(1)}\s+from\s+"(\./[^"]+)"', text)
    if not imp:
        return []
    start = resolve_import(CONTAINER, imp.group(1))
    found, todo = [], [start] if start else []
    root = COMPONENTS.resolve()
    shared = ELEMENTS.resolve()
    while todo:
        f = todo.pop()
        if (f in found or root not in f.parents or shared in f.parents
                or f == CONTAINER.resolve()):
            continue
        found.append(f)
        for spec in re.findall(r'from\s+"(\.{1,2}/[^"]+)"',
                               f.read_text(encoding="utf-8")):
            dep = resolve_import(f, spec)
            if dep:
                todo.append(dep)
    return sorted(found)


def page_task(directory: str, shots: dict) -> str:
    if "task" in shots:
        return shots["task"]
    reverse = {page: task for task, page in ALIASES.items()}
    return reverse.get(directory, directory)


def shots_path(page: str) -> Path:
    return TASK_PAGES / page / "shots.json"


def load(page: str) -> dict:
    return json.loads(shots_path(page).read_text(encoding="utf-8"))


def sources(page: str, shots: dict | None = None) -> list[Path]:
    shots = shots if shots is not None else load(page)
    task = page_task(page, shots)
    entry = task_entries().get(task, {})
    files = []
    if "defXmlPath" in entry:
        files.append(SERVER / entry["defXmlPath"])
    for key in ("pluginPath", "reportPath"):
        if key in entry and (f := module_file(entry[key])):
            files.append(f)
    files += interface_files(task)
    if "scenario" in shots:
        files.append((shots_path(page).parent / shots["scenario"]).resolve())
    return [f for f in files if f.exists()]


def fingerprint(page: str, shots: dict | None = None) -> dict:
    return {rel(f): blob_id(f) for f in sources(page, shots)}


def normalise(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace("&nbsp;", " ")).strip()


def labels(shots: dict):
    """(shot, label) for every label a shot finds its way by."""
    for shot in shots["shots"]:
        out = shot.get("out", "?")
        section = shot.get("section")
        if section:
            yield out, section["text"] if isinstance(section, dict) else section
        for key in ("from", "until", "through"):
            if isinstance(shot.get(key), str):
                yield out, shot[key]
        for callout in shot.get("callouts", []):
            if isinstance(callout, str):
                yield out, callout
            elif not callout.get("dynamic"):
                # "dynamic": text the page makes from data ("delete O4A"),
                # which no source contains.
                yield out, callout.get("field") or callout.get("text")


@cache
def element_text() -> str:
    """The shared widgets' and data classes' own words. Searched for
    labels but not fingerprinted: every page uses them, so a change to one
    would mark every page stale and the warning would stop meaning anything."""
    shared = sorted(ELEMENTS.rglob("*.tsx"))
    # And the data classes' own labels ("Observed data", "Free R set"), which
    # a field shows when neither its def.xml nor its interface names it.
    shared += sorted((SERVER / "core").glob("CCP4*Data*.py"))
    return " ".join(f.read_text(encoding="utf-8", errors="replace")
                    for f in shared)


def missing_labels(page: str, shots: dict) -> list[str]:
    """Labels no source mentions. Only interface shots are checked: a report
    shows what the program wrote, which is in no source here."""
    text = normalise(" ".join(f.read_text(encoding="utf-8", errors="replace")
                              for f in sources(page, shots)) + " " + element_text())
    interface_outs = {s.get("out") for s in shots["shots"]
                      if "Task interface" in s.get("tabs", [])}
    missing = []
    for out, label in labels(shots):
        if out in interface_outs and label and normalise(label).rstrip(":") not in text:
            missing.append(f"{out}: {label!r}")
    return missing


def pages() -> list[str]:
    return sorted(p.parent.name for p in TASK_PAGES.glob("*/shots.json"))


def staleness(page: str) -> list[str]:
    """Why the page may be out of date; empty if it is not."""
    shots = load(page)
    recorded = shots.get("sources")
    if recorded is None:
        return ["never stamped"]
    now = fingerprint(page, shots)
    reasons = [f"{f} changed" for f in now if f in recorded and recorded[f] != now[f]]
    reasons += [f"{f} is new" for f in now if f not in recorded]
    reasons += [f"{f} is gone" for f in recorded if f not in now]
    reasons += [f"label gone from the interface: {m}" for m in missing_labels(page, shots)]
    return reasons


def stamp(page: str):
    path = shots_path(page)
    shots = load(page)
    shots["sources"] = fingerprint(page, shots)
    # Only the "sources" block is (re)written: the rest of the file is laid
    # out by hand, and a reformat would bury the stamp in noise.
    block = '"sources": ' + json.dumps(shots["sources"], indent=2).replace("\n", "\n  ")
    text = path.read_text(encoding="utf-8")
    text, n = re.subn(r'"sources":\s*\{[^{}]*\}', lambda m: block, text)
    if not n:
        text = re.sub(r"\s*\}\s*$", lambda m: ",\n  " + block + "\n}\n", text)
    assert json.loads(text)["sources"] == shots["sources"]
    path.write_text(text, encoding="utf-8")
    missing = missing_labels(page, shots)
    print(f"{page}: stamped {len(shots['sources'])} sources")
    for m in missing:
        print(f"  but a label is not in them: {m}")


def diff(page: str):
    shots = load(page)
    recorded = shots.get("sources", {})
    for f, now in fingerprint(page, shots).items():
        then = recorded.get(f)
        if then == now:
            continue
        new = (REPO / f).read_text(encoding="utf-8").splitlines(keepends=True)
        old = []
        if then:
            got = subprocess.run(["git", "cat-file", "-p", then], cwd=REPO,
                                 capture_output=True, text=True)
            if got.returncode:
                print(f"{f}: stamped as {then[:12]}, not in this clone")
                continue
            old = got.stdout.splitlines(keepends=True)
        sys.stdout.writelines(difflib.unified_diff(
            old, new, f"{f} (stamped)", f"{f} (now)"))


def check(names: list[str], github: bool) -> int:
    stale = 0
    for page in names or pages():
        reasons = staleness(page)
        if not reasons:
            continue
        stale += 1
        if github:
            where = rel(shots_path(page))
            print(f"::warning file={where},title=Help page {page} may be out of date::"
                  + "%0A".join(reasons)
                  + f"%0ARead: docs/user/tools/stamp.py diff {page}; then re-capture,"
                  f" or stamp it if the page still holds.")
        else:
            print(f"{page}:")
            for r in reasons:
                print(f"  {r}")
    if not github:
        print(f"{stale} of {len(names or pages())} pages stale")
    return stale


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("pages", nargs="*")
    c.add_argument("--github", action="store_true")
    d = sub.add_parser("diff")
    d.add_argument("page")
    s = sub.add_parser("stamp")
    s.add_argument("pages", nargs="+")
    args = parser.parse_args()
    if args.cmd == "check":
        check(args.pages, args.github)  # a warning, never a failure
    elif args.cmd == "diff":
        diff(args.page)
    else:
        for page in args.pages:
            stamp(page)


if __name__ == "__main__":
    main()
