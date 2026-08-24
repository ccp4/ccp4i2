"""A missing reference file must not be silently empty.

`ReferenceGroup.loadFromMedLine` used to record the miss on `self.errReport`
and return, producing an empty <CCP4i2ReportReferenceGroup/>.  Nothing reads
that errReport, so the AceDRG-on-Linux bug rendered a normal-looking report
with no citation and no complaint.  These tests pin the warning.

No CCP4 binaries or DB needed.
"""

import logging

import pytest

from ccp4i2.report.metadata import ReferenceGroup
from ccp4i2.core.citations import NON_CITABLE


def test_missing_reference_warns_naming_the_task(caplog):
    group = ReferenceGroup()
    with caplog.at_level(logging.WARNING, logger="ccp4i2.report.metadata"):
        group.loadFromMedLine("NoSuchTaskWhatsoever")

    warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert warnings, "a missing reference file produced no warning"
    message = warnings[0].getMessage()
    assert "NoSuchTaskWhatsoever" in message, (
        "the warning must name the task, or it is not actionable"
    )
    assert len(group) == 0


def test_non_citable_task_stays_quiet(caplog):
    """Tasks with no citable upstream program are declared, not warned about."""
    quiet = next(iter(sorted(NON_CITABLE)))
    group = ReferenceGroup()
    with caplog.at_level(logging.WARNING, logger="ccp4i2.report.metadata"):
        group.loadFromMedLine(quiet)

    assert not [r for r in caplog.records if r.levelno >= logging.WARNING], (
        f"{quiet} is in NON_CITABLE and must not warn"
    )


def test_resolvable_task_neither_warns_nor_is_empty(caplog):
    group = ReferenceGroup()
    with caplog.at_level(logging.WARNING, logger="ccp4i2.report.metadata"):
        group.loadFromMedLine("refmac")

    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert len(group) > 0, "refmac has a reference file; the group must not be empty"


def test_errreport_still_records_the_miss():
    """The warning is additive - the existing errReport entry stays, so any
    future consumer of errReport keeps working."""
    group = ReferenceGroup()
    group.loadFromMedLine("NoSuchTaskWhatsoever")
    assert len(group.errReport) > 0


def test_report_metadata_stays_django_free():
    """The allowlist lives in core.citations, not bibliography, precisely so
    that importing it here does not drag Django into report rendering.
    bibliography imports ccp4i2.db.models at module level; report/metadata.py
    must stay usable without a configured Django."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; import ccp4i2.report.metadata; "
         "sys.exit(1 if any(m == 'django' or m.startswith('django.') "
         "for m in sys.modules) else 0)"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        "importing ccp4i2.report.metadata pulled in Django:\n" + result.stderr
    )


# ---------------------------------------------------------------------------
# loadFromTask: one key->file resolver, shared with the Bibliography button
# ---------------------------------------------------------------------------

def _report_tasknames():
    """Every TASKNAME declared by a *_report.py, scanned rather than imported."""
    import re
    from pathlib import Path
    from ccp4i2 import I2_TOP

    pattern = re.compile(r"^\s*TASKNAME\s*=\s*['\"]([^'\"]+)['\"]", re.MULTILINE)
    names = set()
    for path in Path(I2_TOP).rglob("*_report.py"):
        names.update(pattern.findall(path.read_text(encoding="utf-8", errors="replace")))
    return names


def test_no_report_task_has_an_empty_bibliography(caplog):
    """The report side must resolve citations the same way the Bibliography
    button does.

    It used not to: the button expanded aliases via TASK_CITES while the report
    went straight to ``{TASKNAME}.medline.txt``, so 52 of 143 report classes
    rendered an empty references section for programs cited perfectly well
    elsewhere in the same application. Nothing may regress that to >0.
    """
    from ccp4i2.core.citations import NON_CITABLE

    caplog.set_level(logging.CRITICAL)  # silence the burn-down warning here
    empty = []
    for taskname in sorted(_report_tasknames()):
        if taskname in NON_CITABLE:
            continue
        group = ReferenceGroup()
        group.loadFromTask(taskname)
        if len(group) == 0:
            empty.append(taskname)
    assert not empty, (
        "report classes whose bibliography renders empty; either map them in "
        f"core.citations.TASK_CITES or declare them NON_CITABLE: {empty}"
    )


def test_load_from_task_expands_aliases():
    """A task whose citation is filed under other keys still gets them."""
    group = ReferenceGroup()
    group.loadFromTask("ShelxCD")  # -> shelxc + shelxd
    assert len(group) >= 2, "alias expansion did not pull in both programs"


def test_software_citation_without_a_journal_line_is_kept():
    """BUSTER's entry is a program, not a paper: AU/TI/URL and no SO. Requiring
    SO silently dropped it, which is the same two-parser divergence as the alias
    map, one layer down."""
    group = ReferenceGroup()
    group.loadFromTask("buster")
    assert len(group) > 0, "a software citation with no journal line was dropped"


def test_rendering_a_bibliography_does_not_import_django():
    """Stronger than the import-time check above: TASK_CITES and NON_CITABLE
    both live in core.citations precisely so that *calling* loadFromTask does
    not reach into bibliography, which imports ccp4i2.db.models. A lazy import
    inside the method would pass the import-time test and still drag Django into
    report rendering."""
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; from ccp4i2.report.metadata import ReferenceGroup; "
         "g = ReferenceGroup(); g.loadFromTask('phaser_MR_AUTO'); "
         "assert len(g) > 0, 'expected references'; "
         "sys.exit(1 if any(m == 'django' or m.startswith('django.') "
         "for m in sys.modules) else 0)"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        "rendering a bibliography pulled in Django:\n" + result.stdout + result.stderr
    )
