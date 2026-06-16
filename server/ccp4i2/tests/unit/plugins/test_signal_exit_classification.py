"""
Unit tests for describe_signal_exit() -- the classifier that distinguishes a
native crash (signal death / Windows NTSTATUS) from a clean non-zero program
exit. Pure Python, no CCP4 binaries needed.

This is the mechanism that lets every subprocess boundary in CCP4i2 (external
binary execution, and the worker's per-job subprocess) report "the process
crashed in native code" as its own class of failure, instead of laundering a
SIGSEGV into a generic "exited non-zero".
"""
import signal

import pytest

from ccp4i2.core.CCP4PluginScript import describe_signal_exit


def test_clean_exit_is_not_a_crash():
    assert describe_signal_exit(0) is None


@pytest.mark.parametrize("code", [1, 2, 3, 42, 255])
def test_ordinary_nonzero_exit_is_not_a_crash(code):
    """A program that ran and reported a problem is NOT a crash."""
    assert describe_signal_exit(code) is None


def test_none_returncode_is_not_a_crash():
    assert describe_signal_exit(None) is None


@pytest.mark.parametrize(
    "sig,expected_name",
    [
        (signal.SIGSEGV, "SIGSEGV"),
        (signal.SIGABRT, "SIGABRT"),
        (signal.SIGILL, "SIGILL"),
    ],
)
def test_posix_signal_death_is_classified(sig, expected_name):
    """POSIX subprocess reports signal death as a negative return code."""
    desc = describe_signal_exit(-int(sig))
    assert desc is not None
    assert expected_name in desc
    assert f"signal {int(sig)}" in desc


def test_unknown_negative_code_still_flagged():
    """An unrecognised signal number is still reported as a crash."""
    desc = describe_signal_exit(-99)
    assert desc is not None
    assert "99" in desc


@pytest.mark.parametrize(
    "code,fragment",
    [
        (0xC0000005, "access violation"),
        (0xC000001D, "illegal instruction"),
        (0xC00000FD, "stack overflow"),
    ],
)
def test_windows_ntstatus_crash_is_classified(code, fragment):
    """Windows native crashes surface as large unsigned NTSTATUS exit codes."""
    desc = describe_signal_exit(code)
    assert desc is not None
    assert fragment in desc
    assert f"0x{code:08X}" in desc


def test_unknown_ntstatus_code_still_flagged():
    desc = describe_signal_exit(0xC0000000)
    assert desc is not None
    assert "abnormal termination" in desc
