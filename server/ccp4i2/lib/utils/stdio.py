"""Keep print() from failing a job over a character the console cannot encode.

On Windows a job's stdout is cp1252, so ``print()`` of any text holding a
Greek letter or an emoji raises UnicodeEncodeError. Inside a wrapper that is
a failed job: aimless_pipe (emoji in a status line), then freerflag, whose
annotation carries ``cell.guiLabel()`` and so the alpha, beta and gamma of
the cell (2026-10-22). Fixing each print does not stop the next one, so job
entry points call :func:`tolerate_unencodable_output` before running any
plugin: output keeps its encoding and an unencodable character is written as
an escape (``\\u03b1``) instead of raising.
"""
import sys


def tolerate_unencodable_output():
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:  # replaced by a non-TextIOWrapper (pytest capture)
            continue
        try:
            reconfigure(errors="backslashreplace")
        except (ValueError, OSError):  # closed or detached stream
            pass
