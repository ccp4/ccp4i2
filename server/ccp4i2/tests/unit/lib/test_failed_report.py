"""What the user sees when a report class raises.

Previously: a panel saying "Report generation failed ... See server logs for
full traceback", and nothing else — no line number, and no sign of the output
files, which were gleaned long before the report class ran and are the reason
anyone opened the job. `phaser_MR_PAK_report` rendered exactly that for every
job it ever ran on.

Now the panel names the line that raised, carries the traceback, and keeps the
ordinary input/output file sections.
"""

import xml.etree.ElementTree as ET

import pytest

from ccp4i2.lib.utils.reporting import i2_report
from ccp4i2.lib.utils.reporting.i2_report import (
    failed_report,
    report_is_failure,
    simple_failed_report,
)


def raising(kind=AttributeError, message="no attribute 'xpath'"):
    """An exception with a real traceback, raised at a known line."""
    try:
        raise kind(message)
    except kind as err:
        return err


def diagnostic(root):
    return root.find("CCP4i2ReportDiagnostics/CCP4i2ReportDiagnostic")


def traceback_text(root):
    node = root.find("CCP4i2ReportFold/CCP4i2ReportPre")
    return node.text if node is not None else None


class TestTheFailureIsMarked:
    def test_the_root_says_so(self):
        assert failed_report("boom", "acorn").get("reportFailed") == "true"

    def test_and_the_predicate_reads_it(self):
        assert report_is_failure(failed_report("boom", "acorn")) is True

    def test_a_report_merely_containing_the_old_phrase_is_not_a_failure(self):
        """The predicate used to be a substring search over the markup."""
        innocent = ET.fromstring(
            "<CCP4i2Reportacorn_report>"
            "<CCP4i2ReportPre>No report because the map was empty</CCP4i2ReportPre>"
            "</CCP4i2Reportacorn_report>"
        )
        assert report_is_failure(innocent) is False


class TestTheDiagnosticCard:
    def test_it_names_the_error(self):
        card = diagnostic(failed_report("Report generation failed", "acorn",
                                        exception=raising()))
        assert card is not None
        assert card.get("level") == "error"
        assert "no attribute 'xpath'" in card.get("message")

    def test_it_names_the_line_that_raised(self):
        card = diagnostic(failed_report("boom", "acorn", exception=raising()))
        location = card.get("location")
        assert "test_failed_report.py" in location
        assert " in raising()" in location, "the deepest frame, not the shallowest"

    def test_the_code_is_machine_readable(self):
        card = diagnostic(failed_report("boom", "acorn", code="PROGRAM_XML_NOT_FOUND"))
        assert card.get("code") == "PROGRAM_XML_NOT_FOUND"

    def test_there_is_still_a_card_with_no_exception(self):
        card = diagnostic(failed_report("No program XML found", "acorn"))
        assert card.get("message") == "No program XML found"
        assert card.get("location") is None


class TestTheTraceback:
    def test_it_is_in_the_panel_not_the_server_log(self):
        text = traceback_text(failed_report("boom", "acorn", exception=raising()))
        assert "Traceback (most recent call last)" in text
        assert "AttributeError" in text

    def test_details_come_through_too(self):
        text = traceback_text(
            failed_report("boom", "acorn", details="Task: acorn\nJob: 12")
        )
        assert "Job: 12" in text

    def test_markup_in_the_traceback_survives_serialisation(self):
        """The Pre element renders innerHTML, so `<module>` must be escaped."""
        report = failed_report("boom", "acorn", exception=raising(message="<module>"))
        serialised = ET.tostring(report, encoding="unicode")
        assert "&lt;module&gt;" in serialised

    def test_no_fold_when_there_is_nothing_to_put_in_it(self):
        assert traceback_text(failed_report("boom", "acorn")) is None


class TestTheFileSections:
    """The outputs were gleaned before the report class ran."""

    @pytest.fixture
    def job(self):
        class FakeJob:
            uuid = "1234"

        return FakeJob()

    def test_outputs_and_inputs_are_offered(self, monkeypatch, job):
        monkeypatch.setattr(
            i2_report,
            "get_report_job_info",
            lambda _uuid: {
                "outputfiles": [{"fileId": "out-1"}],
                "inputfiles": [{"fileId": "in-1"}],
                "importedfiles": [],
            },
        )
        report = failed_report("boom", "acorn", exception=raising(), job=job)
        serialised = ET.tostring(report, encoding="unicode")

        assert report.find("CCP4i2ReportOutputData") is not None
        assert report.find("CCP4i2ReportInputData") is not None
        assert "input_file_out-1" in serialised
        assert "input_file_in-1" in serialised

    def test_an_empty_section_is_left_out(self, monkeypatch, job):
        """A heading over nothing is worse than no heading."""
        monkeypatch.setattr(
            i2_report,
            "get_report_job_info",
            lambda _uuid: {"outputfiles": [], "inputfiles": [], "importedfiles": []},
        )
        report = failed_report("boom", "acorn", job=job)
        assert report.find("CCP4i2ReportOutputData") is None
        assert report.find("CCP4i2ReportInputData") is None

    def test_a_failure_in_here_does_not_become_the_failure(self, monkeypatch, job):
        """This runs on the error path. It is not allowed to raise on it."""

        def unavailable(_uuid):
            raise RuntimeError("database has gone away")

        monkeypatch.setattr(i2_report, "get_report_job_info", unavailable)
        report = failed_report("boom", "acorn", exception=raising(), job=job)

        assert report_is_failure(report)
        assert diagnostic(report) is not None, "the original failure still reported"

    def test_no_job_means_no_sections(self):
        report = failed_report("boom", "acorn")
        assert report.find("CCP4i2ReportOutputData") is None


def test_the_deprecated_entry_point_still_works():
    report = simple_failed_report("boom", "acorn", details="context")
    assert report_is_failure(report)
    assert "context" in traceback_text(report)
