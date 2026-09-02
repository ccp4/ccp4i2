"""What the user sees when a report class raises.

Previously: a panel saying "Report generation failed ... See server logs for
full traceback", and nothing else — no line number, and no sign of the output
files, which were gleaned long before the report class ran and are the reason
anyone opened the job. `phaser_MR_PAK_report` rendered exactly that for every
job it ever ran on.

The failure is now carried as an ``errorReportList``, the same shape a job's
``diagnostic.xml`` has, so the Diagnostics panel's own renderer draws it and
there is one such renderer rather than two. The tags themselves are pinned
against that panel in ``tests/unit/validation/test_diagnostic_xml_contract.py``;
what is asserted here is that the failure path produces them at all.
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


def error_reports(root):
    return root.findall("CCP4i2ReportErrorReports/errorReportList/errorReport")


def field(report, tag):
    node = report.find(tag)
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


class TestTheErrorReport:
    def test_there_is_exactly_one_when_the_job_left_none(self):
        assert len(error_reports(failed_report("boom", "acorn"))) == 1

    def test_it_names_the_error(self):
        report = error_reports(failed_report("boom", "acorn", exception=raising()))[0]
        assert field(report, "severityName") == "ERROR"
        assert "no attribute 'xpath'" in field(report, "stack")

    def test_it_names_the_line_that_raised(self):
        report = error_reports(failed_report("boom", "acorn", exception=raising()))[0]
        location = field(report, "name")
        assert "test_failed_report.py" in location
        assert " in raising()" in location, "the deepest frame, not the shallowest"

    def test_it_prefers_the_deepest_ccp4i2_frame_to_the_deepest_frame(self):
        """An exception raised inside a library must still name our line.

        A zero-length program.xml reports `ElementTree.py:573 in parse()`,
        which is true and of no use to anyone; the caller in i2_report.py is
        the line to go and read.
        """
        import xml.etree.ElementTree as element_tree

        try:
            element_tree.fromstring("")
        except element_tree.ParseError as err:
            reports = error_reports(failed_report("b", "acorn", exception=err))
            location = field(reports[0], "name")

        assert "ElementTree.py" not in location
        assert "test_failed_report.py" in location
        # The checkout path is not part of the answer, and a checkout that is
        # itself called ccp4i2 contains the marker twice.
        assert location.startswith("ccp4i2/tests/"), location

    def test_the_code_is_machine_readable_and_explained(self):
        report = error_reports(
            failed_report("boom", "acorn", code="PROGRAM_XML_NOT_FOUND")
        )[0]
        assert field(report, "code") == "PROGRAM_XML_NOT_FOUND"
        assert "wrote no program.xml" in field(report, "description")

    def test_an_unknown_code_falls_back_to_the_reason(self):
        report = error_reports(failed_report("Something odd", "acorn", code="XYZ"))[0]
        assert field(report, "description") == "Something odd"

    def test_details_come_through(self):
        report = error_reports(
            failed_report("boom", "acorn", details="Task: acorn\nJob: 12")
        )[0]
        assert "Job: 12" in field(report, "details")

    def test_the_traceback_is_in_the_panel_not_the_server_log(self):
        report = error_reports(failed_report("boom", "acorn", exception=raising()))[0]
        assert "Traceback (most recent call last)" in field(report, "stack")

    def test_no_stack_when_there_was_no_exception(self):
        """CErrorReport omits an empty stack, and the panel reads that as none.

        "No program XML found" has no traceback to show: there was no
        exception, only an absence.
        """
        assert not field(error_reports(failed_report("boom", "acorn"))[0], "stack")

    def test_markup_in_the_traceback_survives_serialisation(self):
        report = failed_report("boom", "acorn", exception=raising(message="<module>"))
        assert "&lt;module&gt;" in ET.tostring(report, encoding="unicode")


class TestTheJobsOwnDiagnostics:
    """A failing report and a failing job are different events."""

    @pytest.fixture
    def job(self, tmp_path, monkeypatch):
        monkeypatch.setattr(i2_report, "get_report_job_info", lambda _uuid: {})

        class FakeJob:
            uuid = "1234"
            directory = tmp_path

        return FakeJob()

    def test_they_are_folded_in_beside_the_rendering_failure(self, job):
        (job.directory / "diagnostic.xml").write_text(
            "<errorReportList><errorReport>"
            "<class>cmapcoeff</class><code>45</code>"
            "<details>Error in processing output files</details>"
            "</errorReport></errorReportList>"
        )

        reports = error_reports(
            failed_report("boom", "acorn", exception=raising(), job=job)
        )

        assert len(reports) == 2
        assert field(reports[0], "code") != "45", "the rendering failure comes first"
        assert field(reports[1], "code") == "45"

    def test_a_job_that_left_none_adds_none(self, job):
        (job.directory / "diagnostic.xml").write_text("<errorReportList />")
        assert len(error_reports(failed_report("boom", "acorn", job=job))) == 1

    def test_no_diagnostic_file_at_all_is_fine(self, job):
        assert len(error_reports(failed_report("boom", "acorn", job=job))) == 1

    def test_an_unreadable_diagnostic_xml_is_not_the_end_of_it(self, job):
        (job.directory / "diagnostic.xml").write_text("<not-xml")

        reports = error_reports(failed_report("boom", "acorn", job=job))

        assert len(reports) == 1, "the rendering failure still reported"


class TestTheFileSections:
    """The outputs were gleaned before the report class ran."""

    @pytest.fixture
    def job(self, tmp_path):
        class FakeJob:
            uuid = "1234"
            directory = tmp_path

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
        assert error_reports(report), "the original failure still reported"

    def test_no_job_means_no_sections(self):
        report = failed_report("boom", "acorn")
        assert report.find("CCP4i2ReportOutputData") is None


def test_the_deprecated_entry_point_still_works():
    report = simple_failed_report("boom", "acorn", details="context")
    assert report_is_failure(report)
    assert "context" in field(error_reports(report)[0], "details")
