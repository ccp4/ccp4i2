"""A failed regeneration must not destroy a report that once rendered.

report_xml.xml is written once, when a job finishes, and served from then on.
It is less a cache than a record: the rendering a report class produced from
the program.xml its wrapper wrote at the time.

Report classes and the program.xml they read change together, so a newer
report class can fail on an older job's output through no fault of the job.
When that happens the old rendering is the only one that job will ever have,
and must survive.
"""

import uuid
from types import SimpleNamespace

import pytest

from ccp4i2.db import models
from ccp4i2.lib.utils.jobs import reports
from ccp4i2.lib.utils.reporting.i2_report import failed_report

GOOD_REPORT = b"<report><CCP4i2ReportFold label='xia2 report'/></report>"


@pytest.fixture(name="finished_job")
def finished_job_fixture(tmp_path):
    """Enough of a job for get_job_report_xml, without touching the database."""
    return SimpleNamespace(
        status=models.Job.Status.FINISHED,
        directory=tmp_path,
        uuid=uuid.uuid4(),
    )


def test_cached_report_survives_a_regeneration_that_fails(
    finished_job, monkeypatch
):
    cached = finished_job.directory / "report_xml.xml"
    cached.write_bytes(GOOD_REPORT)

    def failing_report(_job):
        # A real failure panel, so this test tracks how a failure is actually
        # signalled rather than restating it. It used to build one whose error
        # text said "Report generation failed", which is what the predicate
        # searched the markup for.
        return failed_report("Report generation failed", "xia2")

    monkeypatch.setattr(reports, "generate_job_report", failing_report)

    result = reports.get_job_report_xml(finished_job, regenerate=True)

    assert result.success
    assert cached.exists(), "the only rendering this job has was deleted"
    assert result.data == GOOD_REPORT, "an error report was served over a good one"


def test_regeneration_that_succeeds_replaces_the_cache(finished_job, monkeypatch):
    cached = finished_job.directory / "report_xml.xml"
    cached.write_bytes(GOOD_REPORT)

    def better_report(_job):
        import xml.etree.ElementTree as ET

        root = ET.Element("report")
        ET.SubElement(root, "CCP4i2ReportFileLink", relativePath="xia2_i2.html")
        return root

    monkeypatch.setattr(reports, "generate_job_report", better_report)

    result = reports.get_job_report_xml(finished_job, regenerate=True)

    assert result.success
    assert b"xia2_i2.html" in result.data
    assert b"xia2_i2.html" in cached.read_bytes()


def test_failed_regeneration_with_no_cache_returns_the_error(
    finished_job, monkeypatch
):
    """With nothing to fall back on, the error must still reach the caller."""

    def failing_report(_job):
        return failed_report(
            "No program XML found", "xia2", details="of missing data"
        )

    monkeypatch.setattr(reports, "generate_job_report", failing_report)

    result = reports.get_job_report_xml(finished_job, regenerate=True)

    assert result.success
    assert b"No program XML found" in result.data
    assert b"of missing data" in result.data, "the panel still says why"
    assert not (finished_job.directory / "report_xml.xml").exists()


def test_a_report_whose_own_text_mentions_failure_is_still_cached(
    finished_job, monkeypatch
):
    """The predicate is an attribute now, not a substring of the prose.

    A working report is entitled to contain the words the old check looked
    for -- "no report because", say, in a fold explaining a skipped step --
    and was previously refused a cache entry for saying so.
    """

    def wordy_report(_job):
        import xml.etree.ElementTree as ET

        root = ET.Element("report")
        ET.SubElement(root, "CCP4i2ReportPre").text = (
            "No report because the anomalous step was skipped"
        )
        return root

    monkeypatch.setattr(reports, "generate_job_report", wordy_report)

    result = reports.get_job_report_xml(finished_job, regenerate=True)

    assert result.success
    assert (finished_job.directory / "report_xml.xml").exists(), (
        "a good report was refused a cache entry for its choice of words"
    )
