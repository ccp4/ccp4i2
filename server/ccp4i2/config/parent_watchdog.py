"""Exit with the process that launched us, when it dies without saying so.

The desktop app spawns uvicorn (two workers) and kills it on quit. That covers
a quit. It does not cover Electron crashing, being killed from a debugger, or
SIGKILL: no handler runs, nothing is signalled, and the uvicorn tree lives on,
re-parented to init or ``systemd --user``, bound to a port nobody will use
again. Do that a few times in a day and a machine carries a stack of orphaned
backends, each with two workers, eating memory (Paul Bond's report).

An orphan cannot be told to die, so it has to notice for itself. When the
launcher passes its pid in ``CCP4I2_PARENT_PID``, every uvicorn worker polls
that pid. The moment it is gone the worker asks the uvicorn supervisor -- its
own parent -- to shut down, which ends the whole tree in order; and exits
itself, so a supervisor that has already gone does not matter either.

Only the ASGI entry point starts this, so ``i2run`` and the job subprocesses
the server launches never do: those are meant to outlive a closed window.
"""
import logging
import os
import signal
import sys
import threading
import time

logger = logging.getLogger(f"ccp4i2:{__name__}")

ENV_VAR = "CCP4I2_PARENT_PID"
DESKTOP_MARKER = "CCP4I2_LOCAL_SESSION_TOKEN"
POLL_SECONDS = 2.0


def pid_is_alive(pid: int) -> bool:
    """Whether a process with this pid exists. Works on POSIX and Windows."""
    if pid <= 0:
        return False
    if sys.platform == "win32":
        import ctypes

        SYNCHRONIZE = 0x00100000
        STILL_ACTIVE = 259
        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        handle = kernel32.OpenProcess(SYNCHRONIZE | 0x0400, False, pid)
        if not handle:
            return False
        try:
            code = ctypes.c_ulong()
            if kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return code.value == STILL_ACTIVE
            return True
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists, owned by someone else
    return True


def _on_parent_gone(parent_pid: int) -> None:
    logger.warning(
        "Launcher process %d is gone; shutting this server down with it", parent_pid
    )
    supervisor = os.getppid()
    # Ask the uvicorn supervisor to stop: it ends every worker in order and
    # exits. Guard against the reparented case, where ppid is init.
    if supervisor > 1:
        try:
            os.kill(supervisor, signal.SIGTERM if sys.platform != "win32" else signal.CTRL_BREAK_EVENT)  # type: ignore[attr-defined]
        except Exception:  # noqa: BLE001 -- best effort, we exit regardless
            pass
    os._exit(0)


def _watch(parent_pid: int, poll_seconds: float, on_gone) -> None:
    while True:
        time.sleep(poll_seconds)
        if not pid_is_alive(parent_pid):
            on_gone(parent_pid)
            return


def start_from_environment(
    poll_seconds: float = POLL_SECONDS, on_gone=_on_parent_gone
) -> threading.Thread | None:
    """Start the watchdog if ``CCP4I2_PARENT_PID`` names a pid; else do nothing.

    ``on_gone`` is what to do when the pid disappears; the default shuts the
    server down. Injectable so the watching can be tested without exiting.
    """
    # Two locks, so no deployed server (Azure, Materia, a lab's shared host)
    # can ever arm this: the launcher's pid, which only the desktop spawn
    # sets, AND the desktop session token, which marks the process as one the
    # desktop app owns. A container has neither.
    if not os.environ.get(DESKTOP_MARKER):
        return None
    raw = os.environ.get(ENV_VAR, "").strip()
    if not raw:
        return None
    try:
        parent_pid = int(raw)
    except ValueError:
        logger.warning("%s=%r is not a pid; not watching", ENV_VAR, raw)
        return None
    if parent_pid == os.getpid() or not pid_is_alive(parent_pid):
        # Already gone before we started, or nonsense: nothing to wait for. Do
        # not exit here -- a stale variable in an inherited environment must
        # not kill a server someone started deliberately.
        return None
    thread = threading.Thread(
        target=_watch,
        args=(parent_pid, poll_seconds, on_gone),
        name="ccp4i2-parent-watchdog",
        daemon=True,
    )
    thread.start()
    logger.info("Watching launcher pid %d; this server exits when it does", parent_pid)
    return thread
