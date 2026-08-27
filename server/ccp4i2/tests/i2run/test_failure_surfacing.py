"""
Deliberate failures against real wrappers: does the user learn *why*?

The rest of this directory runs tasks that succeed. That is the right default
--- but it means the error path, which is the one a user meets on a bad day,
was never exercised on purpose. These tests hand real wrappers input they
cannot use and assert on the artefact the Diagnostics panel reads.

They are deliberately the cheap end of failure: a job given a file it cannot
parse fails in its first seconds, so nothing here refines, phases or scales.
The inputs are written on the spot --- no downloads.

What is asserted is the contract, not any program's exact wording:

- the job records a severity-4 report (it does not finish quietly);
- the report says something specific, not only that the task failed;
- when the failure happened in a subjob, the report names that subjob.
"""
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from .utils import demoData, i2run


def _reports(job: Path):
    """The severity-4 entries of a job's diagnostic.xml."""
    path = job / "diagnostic.xml"
    assert path.exists(), f"a job that failed wrote no diagnostic.xml in {job}"
    out = []
    for report in ET.parse(path).findall(".//errorReport"):
        def text(tag):
            element = report.find(tag)
            return (element.text or "") if element is not None else ""
        if int(text("severity") or 0) >= 4:
            out.append({tag: text(tag) for tag in
                        ("class", "code", "details", "name", "description", "stack")})
    return out


def _said(reports):
    return " ".join(r["details"] + r["description"] + r["stack"] for r in reports)


@pytest.fixture(name="unparsable_pdb")
def unparsable_pdb_fixture(tmp_path):
    path = tmp_path / "not_really_a_model.pdb"
    path.write_text("this file is not a model\nand never was\n", encoding="utf-8")
    return str(path)


@pytest.fixture(name="unparsable_mtz")
def unparsable_mtz_fixture(tmp_path):
    path = tmp_path / "not_really_reflections.mtz"
    path.write_bytes(b"MTZ but not really\x00\x00nonsense\n")
    return str(path)


def test_a_wrapper_given_an_unreadable_model_says_so(unparsable_pdb):
    """A single wrapper: the failure is its own, and it must explain itself."""
    args = ["editbfac", "--XYZIN", unparsable_pdb]
    with i2run(args, allow_errors=True) as job:
        reports = _reports(job)
        assert reports, "the job failed and reported nothing"
        assert len(_said(reports)) > 30, (
            f"nothing said beyond the fact of failure: {reports}")


def test_a_failure_in_a_subjob_names_the_subjob(unparsable_mtz):
    """C3: aimless_pipe's first step cannot read the file; the pipeline must
    report that step's complaint, not merely that a step failed."""
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={unparsable_mtz}"]
    with i2run(args, allow_errors=True) as job:
        reports = _reports(job)
        assert reports, "the pipeline failed and reported nothing"
        attributed = [r for r in reports if r["name"].startswith("job_")]
        assert attributed, (
            "no report named the subjob it came from; the panel would say "
            f"only that the pipeline failed: {reports}")


def test_a_failing_subjob_leaves_its_own_diagnostic(unparsable_mtz):
    """The child's diagnostics live in the child's directory, where the job
    tree in the GUI can find them."""
    args = ["aimless_pipe", "--UNMERGEDFILES", f"file={unparsable_mtz}"]
    with i2run(args, allow_errors=True) as job:
        children = sorted(job.glob("job_*/diagnostic.xml"))
        assert children, f"no subjob wrote a diagnostic.xml under {job}"
