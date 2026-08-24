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
