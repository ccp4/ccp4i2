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
        import xml.etree.ElementTree as ET

        root = ET.Element("report")
        ET.SubElement(root, "error").text = "Report generation failed"
        return root

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
        import xml.etree.ElementTree as ET

        root = ET.Element("report")
        ET.SubElement(root, "error").text = "No report because of missing data"
        return root

    monkeypatch.setattr(reports, "generate_job_report", failing_report)

    result = reports.get_job_report_xml(finished_job, regenerate=True)

    assert result.success
    assert b"No report because" in result.data
    assert not (finished_job.directory / "report_xml.xml").exists()
