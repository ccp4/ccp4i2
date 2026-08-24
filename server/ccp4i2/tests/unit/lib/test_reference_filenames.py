"""Bibliography files must resolve on a case-sensitive filesystem.

Two callers name the same file by different conventions: a report class asks by
its ``TASKNAME`` (``Acedrg``) and the bibliography builder asks by citation key
(``acedrg``). macOS and Windows resolve both, because the filesystem ignores
case. Linux resolves only the one that matches on disk, and the other gets *no*
references with no error -- the report simply renders without its citation.
That is how the AceDRG citation went missing from every Linux report until the
report fixtures first ran on a CI runner.

``CCP4Utils.findReferenceFile`` matches against the real directory listing so
both conventions work everywhere. These tests guard that, and they must never
use ``Path.exists()`` / ``is_file()``: those consult the filesystem, so on a
developer's Mac they would answer yes to a name that fails on Linux -- passing
on exactly the machines where the bug is invisible.
"""

import re
from pathlib import Path

import pytest

from ccp4i2 import I2_TOP
from ccp4i2.core.CCP4Utils import findReferenceFile

I2_ROOT = Path(I2_TOP)
REFERENCES = I2_ROOT / "references"

TASKNAME_RE = re.compile(r"^\s*TASKNAME\s*=\s*['\"]([^'\"]+)['\"]", re.MULTILINE)


def _declared_tasknames():
    """Every TASKNAME declared by a *_report.py module, with its source file.

    Scanning the source rather than importing keeps this fast and covers tasks
    whose report modules need CCP4 to import.
    """
    found = {}
    for report in sorted(I2_ROOT.rglob("*_report.py")):
        try:
            text = report.read_text(encoding="utf-8", errors="replace")
        except OSError:  # pragma: no cover - unreadable file
            continue
        for name in TASKNAME_RE.findall(text):
            found.setdefault(name, report)
    return found


def _real_names():
    """Filenames as the directory actually spells them."""
    return {p.name for p in REFERENCES.glob("*.medline.txt")}


def test_some_report_tasknames_were_found():
    # Guards the scan itself: a regex that silently matched nothing would make
    # the real test below vacuously pass.
    assert len(_declared_tasknames()) > 20


def test_reference_files_resolve_case_sensitively():
    """A TASKNAME whose file exists only in a different case must still resolve.

    Without the resolver this is the Linux-only failure: the file is there, the
    lookup misses it, and the bibliography comes out empty.
    """
    real = _real_names()
    unresolved = []
    for taskname, report in sorted(_declared_tasknames().items()):
        wanted = f"{taskname}.medline.txt"
        if wanted in real:
            continue  # exact match, fine on every platform
        if wanted.lower() not in {n.lower() for n in real}:
            continue  # no bibliography for this task at all, which is allowed
        if findReferenceFile(taskname) is None:
            unresolved.append(
                f"  {report.relative_to(I2_ROOT)}: TASKNAME={taskname!r} cannot "
                f"resolve {wanted!r}"
            )
    assert not unresolved, (
        "reference files present on disk but unreachable on a case-sensitive "
        "filesystem:\n" + "\n".join(unresolved)
    )


def test_resolver_ignores_case_and_reports_absence():
    # 'Acedrg' is the report TASKNAME, 'acedrg' the citation key, one file.
    by_taskname = findReferenceFile("Acedrg")
    by_key = findReferenceFile("acedrg")
    assert by_taskname is not None and by_key is not None
    assert by_taskname.name == by_key.name
    assert findReferenceFile("no_such_task_xyz") is None


def test_reference_directory_has_no_case_collisions():
    """Two files differing only in case cannot both survive a checkout on a
    case-insensitive filesystem, and would make the resolver's answer arbitrary."""
    seen = {}
    collisions = []
    for name in sorted(_real_names()):
        previous = seen.setdefault(name.lower(), name)
        if previous != name:
            collisions.append(f"  {previous} vs {name}")
    assert not collisions, (
        "reference filenames collide case-insensitively:\n" + "\n".join(collisions)
    )
