"""GenericReport — the fallback report for tasks with no report class.

Tasks such as the Import* family ship no ``*_report.py``.  They used to render
as "No report because: No report class found for task"; they now fall back to
GenericReport, which shows the same standard sections Qt-i2 showed: a title,
the input files, the output files and the job details.
"""

import xml.etree.ElementTree as ET

from ccp4i2.report.CCP4ReportParser import GenericReport


JOB_INFO = {
    "status": "Finished",
    "taskname": "ImportCoordinate",
    "tasktitle": "Import coordinate file",
    "jobnumber": 3,
    "jobtitle": "",
    "creationtime": 1_700_000_000.0,
    "finishtime": 1_700_000_060.0,
    "projectname": "test",
    "fileroot": "/tmp/job_3/",
    "inputfiles": [{"fileId": "aaaa-1111"}],
    "outputfiles": [{"fileId": "bbbb-2222"}],
    "importedfiles": [],
}


def _report_xml(standardise=True, **kw):
    report = GenericReport(jobInfo=dict(JOB_INFO), standardise=standardise, **kw)
    return ET.tostring(report.as_data_etree(), encoding="unicode")


def test_generic_report_needs_no_program_xml():
    # Import tasks write no program.xml; the report must not demand one.
    assert GenericReport.USEPROGRAMXML is False


def test_generic_report_has_the_standard_sections():
    xml = _report_xml()
    for element in (
        "CCP4i2ReportTitle",
        "CCP4i2ReportInputData",
        "CCP4i2ReportOutputData",
        "CCP4i2ReportJobDetails",
    ):
        assert element in xml, f"{element} missing from generic report"


def test_generic_report_lists_input_and_output_file_ids():
    xml = _report_xml()
    # The frontend picks file UUIDs out of input_file_<uuid> div ids.
    assert "input_file_aaaa-1111" in xml
    assert "input_file_bbbb-2222" in xml


def test_generic_report_names_the_task():
    # Standardised, the Title bar carries the name; unstandardised, a text line.
    assert "Import coordinate file" in _report_xml()
    assert "Import coordinate file" in _report_xml(standardise=False)


def test_generic_report_reports_a_non_finished_status():
    # A running job is not standardised, so the status line is all it has.
    xml = _report_xml(standardise=False, jobStatus="Running")
    assert "Running" in xml
