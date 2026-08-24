"""Bibliography files must match their TASKNAME exactly, not just case-insensitively.

``ReferenceGroup.loadFromMedLine(taskName)`` opens
``references/{taskName}.medline.txt``. macOS and Windows resolve that
case-insensitively, so a file whose case does not match its TASKNAME works on a
developer's machine and silently produces an *empty* bibliography on Linux --
no error, just a missing citation in every report for that task. That is how
``acedrg``/``acedrgNew`` (TASKNAME ``Acedrg``, file ``acedrg.medline.txt``) went
unnoticed: it only showed up when the report fixtures first ran on a Linux CI
runner.

Scanning the source rather than importing the report classes keeps this fast and
lets it cover tasks whose modules need CCP4 to import.
"""

import re
from pathlib import Path

import pytest

from ccp4i2 import I2_TOP

I2_ROOT = Path(I2_TOP)
REFERENCES = I2_ROOT / "references"

TASKNAME_RE = re.compile(r"^\s*TASKNAME\s*=\s*['\"]([^'\"]+)['\"]", re.MULTILINE)


def _declared_tasknames():
    """Every TASKNAME declared by a *_report.py module, with its source file."""
    found = {}
    for report in sorted(I2_ROOT.rglob("*_report.py")):
        try:
            text = report.read_text(encoding="utf-8", errors="replace")
        except OSError:  # pragma: no cover - unreadable file
            continue
        for name in TASKNAME_RE.findall(text):
            found.setdefault(name, report)
    return found


def test_some_report_tasknames_were_found():
    # Guards the scan itself: a regex that silently matches nothing would make
    # the real test below vacuously pass.
    assert len(_declared_tasknames()) > 20


def test_reference_files_match_taskname_case():
    # Compare against the real directory listing, never Path.exists()/is_file():
    # those go through the filesystem, which on macOS and Windows answers
    # case-insensitively -- so the check would pass on exactly the machines
    # where the bug is invisible, which is the whole problem being guarded.
    actual_names = {p.name for p in REFERENCES.glob("*.medline.txt")}
    by_lower = {name.lower(): name for name in actual_names}

    mismatches = []
    for taskname, report in sorted(_declared_tasknames().items()):
        wanted = f"{taskname}.medline.txt"
        if wanted in actual_names:
            continue
        actual = by_lower.get(wanted.lower())
        if actual is not None:
            mismatches.append(
                f"  {report.relative_to(I2_ROOT)}: TASKNAME={taskname!r} looks for "
                f"{wanted!r} but the file is named {actual!r}"
            )

    assert not mismatches, (
        "Bibliography filenames differ from their TASKNAME only by case. These "
        "resolve on macOS/Windows and yield an empty bibliography on Linux:\n"
        + "\n".join(mismatches)
    )


def test_reference_directory_has_no_case_collisions():
    """Two files differing only in case cannot both survive a checkout on a
    case-insensitive filesystem, so one would arrive missing or overwritten."""
    seen = {}
    collisions = []
    for path in sorted(REFERENCES.glob("*.medline.txt")):
        previous = seen.setdefault(path.name.lower(), path.name)
        if previous != path.name:
            collisions.append(f"  {previous} vs {path.name}")
    assert not collisions, "reference filenames collide case-insensitively:\n" + "\n".join(
        collisions
    )
