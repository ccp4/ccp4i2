"""Tests for distinguishing a failed MrParse search from an empty one.

MrParse writes a report either way, so a search program that could not be
executed is otherwise indistinguishable from one that ran and matched
nothing. These pin the three outcomes the wrappers act on.
"""

from ccp4i2.lib.utils.logs.mrparse_log import search_program_failure

# The failure that prompted this: an x86_64-only phmmer in an arm64 CCP4
# build, as MrParse's own log records it.
BAD_ARCHITECTURE_LOG = """\
2026-08-24 14:22:19,001 - root - INFO - Running phmmer search
Traceback (most recent call last):
  File "/opt/ccp4/lib/python3.11/site-packages/mrparse/mr_util.py", line 84, in run_cmd
    out = subprocess.check_output(cmd, **optd)
  File "/opt/ccp4/lib/python3.11/subprocess.py", line 1955, in _execute_child
    raise child_exception_type(errno_num, err_msg, err_filename)
OSError: [Errno 86] Bad CPU type in executable: '/opt/ccp4/libexec/phmmer'
2026-08-24 14:22:20,774 - root - DEBUG - No homologues found
2026-08-24 14:22:20,777 - root - INFO - Wrote MrParse output file: mrparse.html
"""

EMPTY_BUT_HEALTHY_LOG = """\
2026-08-24 14:22:19,001 - root - INFO - Running phmmer search
2026-08-24 14:22:20,774 - root - DEBUG - No homologues found
2026-08-24 14:22:20,777 - root - INFO - Wrote MrParse output file: mrparse.html
"""


def test_reports_the_failure_detail(tmp_path):
    log = tmp_path / "mrparse.log"
    log.write_text(BAD_ARCHITECTURE_LOG)
    assert search_program_failure(log) == (
        "[Errno 86] Bad CPU type in executable: '/opt/ccp4/libexec/phmmer'"
    )


def test_a_search_that_found_nothing_is_not_a_failure(tmp_path):
    log = tmp_path / "mrparse.log"
    log.write_text(EMPTY_BUT_HEALTHY_LOG)
    assert search_program_failure(log) is None


def test_a_missing_binary_is_a_failure(tmp_path):
    log = tmp_path / "mrparse.log"
    log.write_text(
        "INFO - Running phmmer search\n"
        "FileNotFoundError: [Errno 2] No such file or directory: 'phmmer'\n"
    )
    assert search_program_failure(log) == (
        "[Errno 2] No such file or directory: 'phmmer'"
    )


def test_a_nonzero_exit_is_a_failure(tmp_path):
    log = tmp_path / "mrparse.log"
    log.write_text(
        "subprocess.CalledProcessError: Command '['phmmer']' returned non-zero exit status 1.\n"
    )
    assert "returned non-zero exit status 1" in search_program_failure(log)


def test_an_absent_log_is_not_evidence_of_failure(tmp_path):
    assert search_program_failure(tmp_path / "never-written.log") is None


def test_mention_of_an_exception_in_prose_is_not_a_failure(tmp_path):
    """Only a traceback's final line counts, not any line naming an exception."""
    log = tmp_path / "mrparse.log"
    log.write_text(
        "INFO - retrying after OSError: transient network problem, recovered\n"
        "DEBUG - No homologues found\n"
    )
    assert search_program_failure(log) is None
