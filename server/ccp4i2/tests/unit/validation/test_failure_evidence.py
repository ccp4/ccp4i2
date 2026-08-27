"""
What a failed program's output contributes to the error report.

C7 of docs/error-handling-remediation.md. The framework read stderr when a
program exited non-zero and print()ed it to the server console, where no user
can see it, while the Diagnostics panel showed "Process refmac5 exited with
code 1" and nothing more.

The rule these tests hold to: **the report carries evidence, not artefacts**. A
bounded, labelled tail and the name of the file it came from --- never the file.
Seven wrappers embed whole logs into program.xml and produce 0.9 MB files for
it; this must not become the eighth.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.CCP4PluginScript import (
    EVIDENCE_MAX_CHARS,
    _candidate_lines,
    _collapse_progress,
    failure_evidence,
)


def _write(tmp_path, name, text, newline="\n"):
    """Write *text* with the given line ending, byte for byte.

    newline="" on the write side too: these tests are about what survives the
    round trip through a file, and Windows is an equal peer target, so a log
    with CRLF endings has to be exercised as a file rather than as a string.
    """
    path = tmp_path / name
    with open(path, "w", newline="") as handle:
        handle.write(text.replace("\n", newline))
    return str(path)


# --- what a reader actually needs -------------------------------------------

def test_a_traceback_is_surfaced_as_a_candidate(tmp_path):
    """modelcraft's real failure: the cause is in stderr."""
    stderr = _write(tmp_path, "stderr.txt",
                    "Traceback (most recent call last):\n"
                    "  File \"/x/bin/modelcraft\", line 7, in <module>\n"
                    "ModuleNotFoundError: No module named 'modelcraft.scripts'\n")
    evidence = failure_evidence(stderr, None)
    assert "Possibly relevant" in evidence
    assert "Traceback" in evidence
    assert "stderr.txt" in evidence, "the reader must be told where it came from"


def test_a_missing_file_is_surfaced(tmp_path):
    """findmyseq's real failure."""
    stderr = _write(tmp_path, "stderr.txt",
                    "FileNotFoundError: [Errno 2] No such file or directory: "
                    "'/jobs/job_1/hklin.mtz'\n")
    assert "hklin.mtz" in failure_evidence(stderr, None)


def test_candidates_never_claim_to_be_the_cause(tmp_path):
    """crank2's log ends with a warning about SHELXC parameters on a run that
    died for want of the shelxc binary --- adjacent to the truth, not it."""
    log = _write(tmp_path, "log.txt",
                 "*** Starting FA estimation ***\n"
                 "Warning: Virtual parameter binary has not been used by SHELXC.\n")
    evidence = failure_evidence(None, log)
    assert "Possibly relevant" in evidence or "Last lines" in evidence
    for forbidden in ("caused by", "the error was", "because"):
        assert forbidden not in evidence.lower()


# --- noise that would crowd it out ------------------------------------------

def test_progress_frames_are_dropped(tmp_path):
    """nucleofind writes thirty of these to stderr."""
    frames = "".join(
        f"Predicting:  {i*3}%|###   | {i}/30 [00:{i:02d}<00:47,  2.39s/it]\n"
        for i in range(1, 31)
    )
    stderr = _write(tmp_path, "stderr.txt", frames + "Saving grid to .\n")
    evidence = failure_evidence(stderr, None)
    assert "Predicting" not in evidence
    assert "Saving grid" in evidence


def test_carriage_returns_are_collapsed():
    text = "Working:  1%\rWorking: 50%\rWorking: 100% done\nnext line\n"
    assert _collapse_progress(text).splitlines() == ["Working: 100% done",
                                                     "next line"]


# --- Windows is an equal peer target ----------------------------------------
#
# open() defaults to universal newlines, which translates \r, \n and \r\n all
# to \n before any of this code sees them --- so the collapsing above could
# only ever work on a string, never on a file. Reading with newline="" fixes
# that and introduces the opposite hazard: a CRLF line is "text\r", and taking
# the last \r-separated segment of it yields "".

def test_a_crlf_log_is_not_blanked(tmp_path):
    """The whole evidence would be empty lines on Windows."""
    log = _write(tmp_path, "log.txt",
                 "Refinement cycle 1\nRefinement cycle 2\nerror: bad geometry\n",
                 newline="\r\n")
    evidence = failure_evidence(None, log)
    assert "Refinement cycle 2" in evidence
    assert "error: bad geometry" in evidence


def test_a_crlf_log_leaves_no_stray_carriage_returns(tmp_path):
    """They would reach the panel, and diagnostic.xml, as literal characters."""
    log = _write(tmp_path, "log.txt", "one\ntwo\nthree\n", newline="\r\n")
    assert "\r" not in failure_evidence(None, log)


def test_a_progress_bar_written_with_crlf_endings(tmp_path):
    """Both mechanisms at once: \r overwriting within CRLF-terminated lines."""
    text = ("Predicting: 10%|# | 3/30 [00:07<00:64, 2.39s/it]\r"
            "Predicting: 20%|## | 6/30 [00:14<00:57, 2.39s/it]\n"
            "done writing maps\n")
    stderr = _write(tmp_path, "stderr.txt", text, newline="\r\n")
    evidence = failure_evidence(stderr, None)
    assert "Predicting" not in evidence
    assert "done writing maps" in evidence


def test_a_file_read_as_a_file_collapses_overwriting(tmp_path):
    """The case the original test could not reach, because it passed a string."""
    stderr = _write(tmp_path, "stderr.txt", "Working: 1%\rWorking: 99% done\n")
    evidence = failure_evidence(stderr, None)
    assert "Working: 99% done" in evidence
    assert "Working: 1%" not in evidence


def test_the_same_line_is_not_repeated(tmp_path):
    stderr = _write(tmp_path, "stderr.txt", "error: bad input\n" * 20)
    assert _candidate_lines(stderr and open(stderr).read()) == ["error: bad input"]


# --- evidence, not artefacts ------------------------------------------------

def test_a_huge_log_contributes_a_bounded_amount(tmp_path):
    """The whole point. A megabyte of log must not enter the report."""
    log = _write(tmp_path, "log.txt", "refinement cycle line\n" * 50_000)
    evidence = failure_evidence(None, log)
    assert 0 < len(evidence) < 4 * EVIDENCE_MAX_CHARS, len(evidence)


def test_files_are_named_not_copied(tmp_path):
    log = _write(tmp_path, "log.txt", "\n".join(f"line {i}" for i in range(500)))
    evidence = failure_evidence(None, log)
    assert "log.txt" in evidence
    assert "line 0" not in evidence, "only the tail, not the file"
    assert "line 499" in evidence


# --- nothing to say ---------------------------------------------------------

def test_absent_files_produce_nothing(tmp_path):
    assert failure_evidence(str(tmp_path / "nope.txt"), str(tmp_path / "also.txt")) == ""


def test_an_empty_stderr_is_skipped(tmp_path):
    """stderr is empty in 86% of jobs -- CCP4 programs write to the log."""
    stderr = _write(tmp_path, "stderr.txt", "")
    log = _write(tmp_path, "log.txt", "something happened\n")
    evidence = failure_evidence(stderr, log)
    assert "stderr.txt" not in evidence
    assert "log.txt" in evidence


def test_no_arguments_at_all_is_survivable():
    assert failure_evidence(None, None) == ""
