"""Tell a MrParse search that failed from one that found nothing.

MrParse runs its sequence searches through ``mr_util.run_cmd``. When the
search program cannot be executed at all — a missing binary, or one built
for the wrong architecture — the exception is logged and MrParse carries
on to write a report with empty result sets. To the user the outcome is
indistinguishable from an honest "no homologues found": a finished job
with an empty report.

The distinction matters. One means there is nothing to find; the other
means the search never ran and the answer is unknown. The observed case
was an x86_64-only ``phmmer`` in an otherwise arm64 CCP4 build on Apple
Silicon, failing with ``OSError: [Errno 86] Bad CPU type in executable``.
"""

import logging
import re
from pathlib import Path
from typing import Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")

# Exceptions that mean a program could not be run, as opposed to one that
# ran and reported nothing. Matched at the head of a traceback's final line.
EXECUTION_FAILURES = (
    "OSError",
    "FileNotFoundError",
    "PermissionError",
    "NotADirectoryError",
    "subprocess.CalledProcessError",
    "CalledProcessError",
)

_FAILURE_LINE = re.compile(
    r"^(?:%s):\s*(?P<detail>.+)$" % "|".join(re.escape(e) for e in EXECUTION_FAILURES)
)


def search_program_failure(log_path: Path) -> Optional[str]:
    """Return a description of the first search program that failed to run.

    Returns None when the log shows no such failure — including when the
    log is missing or unreadable, since an absent log is not evidence of
    one.
    """
    log_path = Path(log_path)
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as err:
        logger.debug("Could not read %s: %s", log_path, err)
        return None

    for line in lines:
        match = _FAILURE_LINE.match(line.strip())
        if match:
            return match.group("detail").strip()
    return None
