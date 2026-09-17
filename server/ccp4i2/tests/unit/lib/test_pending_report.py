"""A job that is under way and has not written program.xml gets a calm
"running" report, not the failure panel."""

import xml.etree.ElementTree as ET
from types import SimpleNamespace

from ccp4i2.db.models import Job
from ccp4i2.lib.utils.reporting.i2_report import (
    failed_report,
    pending_report,
    report_is_failure,
    report_is_pending,
)


def job(status, number="5"):
    # No database: _file_sections tolerates a job it cannot look up.
    return SimpleNamespace(status=status, number=number, uuid="not-a-real-uuid")


def test_it_is_marked_pending_and_not_a_failure():
    report = pending_report("molrep_map", job(Job.Status.RUNNING))
    assert report.get("reportPending") == "true"
    assert report_is_pending(report) is True
    assert report_is_failure(report) is False
    assert report.find("CCP4i2ReportErrorReports") is None


def test_it_says_the_job_is_running_in_its_own_voice():
    report = pending_report("molrep_map", job(Job.Status.RUNNING))
    assert report.find("CCP4i2ReportTitle").get("title1") == "Job 5 is running"
    text = report.find("CCP4i2ReportText").text
    assert text.startswith("Job 5 is running.")
    assert "not written a report yet" in text


def test_queued_and_remote_jobs_say_so():
    assert "is queued" in pending_report("x", job(Job.Status.QUEUED)).find("CCP4i2ReportTitle").get("title1")
    assert "running remotely" in pending_report("x", job(Job.Status.RUNNING_REMOTELY)).find("CCP4i2ReportTitle").get("title1")


def test_the_failure_panel_is_untouched_for_finished_jobs():
    report = failed_report("No program XML found", "molrep_map", code="PROGRAM_XML_NOT_FOUND")
    assert report_is_failure(report) is True
    assert report_is_pending(report) is False
