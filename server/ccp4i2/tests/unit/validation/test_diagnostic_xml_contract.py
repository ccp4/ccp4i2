"""
The contract between diagnostic.xml and the Diagnostics panel.

C6 of docs/error-handling-remediation.md. The backend has always written
``<class>``; the panel read ``<className>`` and so drew an empty heading over
every error CCP4i2 has ever reported --- the traceback was there, in
``<details>``, under a title reading " - 993".

Both ends are now asserted, here and in
``client/renderer/__tests__/diagnostic-parse.test.ts``. The two tests name the
same tags on purpose: if either side is changed alone, one of them fails.

Pure Python -- no CCP4 binaries needed.
"""
import xml.etree.ElementTree as ET

from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR

#: Every tag the panel reads. Keep in step with the TypeScript test.
TAGS_THE_PANEL_READS = {
    "class",        # the task; shown in the heading
    "code",         # shown beside it
    "name",         # which parameter the report is about
    "description",  # what the code means in general
    "details",      # what went wrong this time
    "stack",        # the traceback, folded away
    "severityName", # ERROR / WARNING, drives the icon and colour
}


def _report(**kwargs):
    report = CErrorReport()
    report.append(**{
        "klass": "refmac",
        "code": 201,
        "details": "Exit code: 5",
        "name": "XYZIN",
        "severity": SEVERITY_ERROR,
        **kwargs,
    })
    return ET.fromstring(ET.tostring(report.getEtree())).find("errorReport")


def _text(element, tag):
    found = element.find(tag)
    return found.text if found is not None else None


def test_a_full_report_carries_everything_the_panel_reads():
    element = _report(description="No geometry dictionary for a ligand",
                      stack="Traceback (most recent call last):\n  ...")
    tags = {child.tag for child in element}
    assert TAGS_THE_PANEL_READS <= tags, (
        f"missing from diagnostic.xml: {TAGS_THE_PANEL_READS - tags}"
    )


def test_the_task_and_the_parameter_are_both_recorded():
    """The heading is built from these two; without them it read ' - 201'."""
    element = _report()
    assert _text(element, "class") == "refmac"
    assert _text(element, "name") == "XYZIN"


def test_details_and_description_are_separate_things():
    """One says what this code means; the other what happened this time."""
    element = _report(description="Ligand has no geometry dictionary")
    assert _text(element, "description") == "Ligand has no geometry dictionary"
    assert _text(element, "details") == "Exit code: 5"


def test_the_traceback_is_kept_apart_from_the_details():
    element = _report(stack="Traceback (most recent call last):\n  boom")
    assert "Traceback" in _text(element, "stack")
    assert "Traceback" not in _text(element, "details")


def test_optional_fields_are_omitted_rather_than_left_empty():
    """An empty <description> would draw a heading over nothing."""
    element = _report()
    assert element.find("description") is None
    assert element.find("stack") is None


def test_severity_is_readable_as_a_word():
    element = _report()
    assert _text(element, "severityName") == "ERROR"
    assert _text(element, "severity") == str(SEVERITY_ERROR)


# --- a code the task does not declare ---------------------------------------
#
# ERROR_CODES is how a task states what its codes mean, and 38 calls in the tree
# pass a code and nothing else -- so an undeclared code there leaves the reader
# a card with no words on it. It must not stop the job either: the report being
# made matters more than its being well-formed.

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript


def _plugin(error_codes):
    plugin = CPluginScript.__new__(CPluginScript)
    plugin.errorReport = CErrorReport()
    plugin.TASKNAME = "sometask"
    plugin.ERROR_CODES = error_codes
    CPluginScript._undeclared_codes_seen.clear()
    return plugin


def _only(plugin):
    return plugin.errorReport.getErrors()[0]


def test_a_declared_code_supplies_its_description():
    plugin = _plugin({201: {"description": "No geometry dictionary"}})
    plugin.appendErrorReport(201, "Exit code: 5")
    assert _only(plugin)["description"] == "No geometry dictionary"


def test_an_undeclared_code_still_records_the_report():
    """Graceful: a fault in the task must not cost the user the report."""
    plugin = _plugin({})
    plugin.appendErrorReport(999, "The program wrote no output")
    error = _only(plugin)
    assert error["code"] == 999
    assert error["details"] == "The program wrote no output"
    assert error["description"] == ""


def test_an_undeclared_code_with_no_details_still_says_something():
    """Helpful: the alternative is a card with no words on it."""
    plugin = _plugin({})
    plugin.appendErrorReport(999)
    details = _only(plugin)["details"]
    assert "999" in details
    assert "CPluginScript" in details or "sometask" in details
    assert "report it" in details


def test_a_declared_code_with_no_details_uses_its_description():
    """The 38 call sites that pass a code alone are the reason for this."""
    plugin = _plugin({200: {"description": "Input file could not be read"}})
    plugin.appendErrorReport(200)
    assert _only(plugin)["description"] == "Input file could not be read"


@pytest.mark.parametrize("codes", [
    {},                                   # not declared at all
    {303: "just a string"},               # declared as the wrong shape
    {303: {"description": ""}},           # declared, described with nothing
])
def test_every_shape_of_missing_description_is_survivable(codes):
    plugin = _plugin(codes)
    plugin.appendErrorReport(303, "something went wrong")
    assert _only(plugin)["details"] == "something went wrong"
    assert _only(plugin)["description"] == ""


def test_the_same_undeclared_code_is_complained_about_once(caplog):
    """A task that reports it on every run should not fill the log."""
    plugin = _plugin({})
    with caplog.at_level("WARNING"):
        for _ in range(5):
            plugin.appendErrorReport(999, "again")
    complaints = [r for r in caplog.records if "ERROR_CODES" in r.getMessage()]
    assert len(complaints) == 1


# --- codes the base class reports on a task's behalf -------------------------

def test_a_shared_code_is_described_even_though_the_task_never_declared_it():
    """The panel showed "servalcat_pipe - job_4/... - 350" with no words above it.

    checkMonomeCoverage is a base-class check reporting code 350; the task did
    not raise it and has no business declaring it, so the reader got a card
    with only the specifics on it.
    """
    plugin = _plugin({})
    plugin.appendErrorReport(350, "No dictionary at all for: DRG")
    assert _only(plugin)["description"] == \
        "A ligand or modified residue has no geometry dictionary"


def test_the_task_still_wins_where_it_has_something_to_say():
    plugin = _plugin({350: {"description": "This task's own words"}})
    plugin.appendErrorReport(350, "No dictionary at all for: DRG")
    assert _only(plugin)["description"] == "This task's own words"


def test_a_code_nobody_declares_is_still_only_a_warning_in_the_log():
    plugin = _plugin({})
    plugin.appendErrorReport(4242, "something happened")
    assert _only(plugin)["description"] == ""
    assert _only(plugin)["details"] == "something happened"
