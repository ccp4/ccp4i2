"""The server notices when the desktop app that launched it is gone.

An orphaned uvicorn tree cannot be told to stop, so it watches the launcher's
pid itself (parent_watchdog). These tests use a throwaway subprocess as the
"launcher" and an injected reaction, so nothing here exits.
"""
import os
import subprocess
import sys
import threading
import time

import pytest

from ccp4i2.config import parent_watchdog as pw


def test_our_own_pid_is_alive():
    assert pw.pid_is_alive(os.getpid())


def test_a_finished_process_is_not_alive():
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    # On POSIX a zombie stays "alive" until reaped; wait() reaped it.
    assert not pw.pid_is_alive(proc.pid)


def test_nonsense_pids_are_not_alive():
    assert not pw.pid_is_alive(0)
    assert not pw.pid_is_alive(-1)


@pytest.fixture(autouse=True)
def desktop_marker(monkeypatch):
    """Every test here runs 'as the desktop', except the one that checks the gate."""
    monkeypatch.setenv(pw.DESKTOP_MARKER, "test-token")


def test_nothing_starts_without_the_variable(monkeypatch):
    monkeypatch.delenv(pw.ENV_VAR, raising=False)
    assert pw.start_from_environment() is None


def test_nothing_starts_outside_the_desktop(monkeypatch):
    """A deployed server (Azure, Materia) has no session token: even a live
    pid in the variable must not arm a watchdog there."""
    monkeypatch.delenv(pw.DESKTOP_MARKER, raising=False)
    monkeypatch.setenv(pw.ENV_VAR, str(os.getpid()))
    assert pw.start_from_environment() is None


def test_a_stale_variable_does_not_start_a_watch(monkeypatch):
    """A dead pid in an inherited environment must not kill a server someone
    started deliberately, so it is ignored rather than acted on."""
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    monkeypatch.setenv(pw.ENV_VAR, str(proc.pid))
    assert pw.start_from_environment() is None


def test_garbage_in_the_variable_is_ignored(monkeypatch):
    monkeypatch.setenv(pw.ENV_VAR, "not-a-pid")
    assert pw.start_from_environment() is None


def test_the_reaction_fires_when_the_launcher_dies(monkeypatch):
    launcher = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
    monkeypatch.setenv(pw.ENV_VAR, str(launcher.pid))
    fired = threading.Event()
    seen = []

    def on_gone(pid):
        seen.append(pid)
        fired.set()

    thread = pw.start_from_environment(poll_seconds=0.05, on_gone=on_gone)
    assert thread is not None and thread.is_alive()
    assert not fired.wait(0.3), "must not fire while the launcher lives"

    launcher.kill()
    launcher.wait()
    assert fired.wait(3.0), "the watchdog did not notice the launcher die"
    assert seen == [launcher.pid]


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signal semantics")
def test_the_default_reaction_signals_the_supervisor(monkeypatch):
    """The default asks its parent (the uvicorn supervisor) to stop, then exits."""
    calls = []
    monkeypatch.setattr(pw.os, "getppid", lambda: 4242)
    monkeypatch.setattr(pw.os, "kill", lambda pid, sig: calls.append((pid, sig)))
    monkeypatch.setattr(pw.os, "_exit", lambda code: calls.append(("exit", code)))
    pw._on_parent_gone(999)
    assert calls == [(4242, pw.signal.SIGTERM), ("exit", 0)]
