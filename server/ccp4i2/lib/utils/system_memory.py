"""How much memory this process can really use.

``psutil.virtual_memory().total`` is the host's memory. Inside a container
that is the node, not the cgroup limit the process will be killed at, so a
"will this fit" estimate computed against it is wrong exactly where it
matters (Materia, first Batch run, 2026-09-27). The cgroup limit, when there
is one, is the smaller number and the true one.
"""
from pathlib import Path
from typing import Optional

# cgroup v2, then v1. "max" (v2) or a huge number (v1) means unlimited.
_CGROUP_LIMIT_FILES = (
    Path("/sys/fs/cgroup/memory.max"),
    Path("/sys/fs/cgroup/memory/memory.limit_in_bytes"),
)
_UNLIMITED_ABOVE = 1 << 60


def cgroup_limit_bytes(files=_CGROUP_LIMIT_FILES) -> Optional[int]:
    """The container's memory limit in bytes, or None when unlimited/absent."""
    for path in files:
        try:
            text = Path(path).read_text().strip()
        except OSError:
            continue
        if text == "max":
            return None
        try:
            value = int(text)
        except ValueError:
            continue
        if value <= 0 or value >= _UNLIMITED_ABOVE:
            return None
        return value
    return None


def host_total_bytes() -> Optional[int]:
    try:
        import psutil
        return int(psutil.virtual_memory().total)
    except Exception:  # noqa: BLE001 -- no psutil, no answer
        return None


def usable_gib(files=_CGROUP_LIMIT_FILES) -> Optional[float]:
    """The smaller of the host's memory and the cgroup limit, in GiB; None if unknown."""
    candidates = [b for b in (host_total_bytes(), cgroup_limit_bytes(files)) if b]
    if not candidates:
        return None
    return min(candidates) / 2 ** 30
