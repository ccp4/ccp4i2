"""Content-based detection of reflection-file format.

The legacy path decided format from the filename extension alone
(``CGenericReflDataFile.getFormat()``), which is unsafe: ``.hkl`` is used by
*both* XDS_ASCII (which gemmi reads) and SHELX (which it does not), and a
mis-named file routes wrongly with no content check. This peeks at the actual
bytes instead — the first step of the single Python ``diagnose_reflection_file``
authority the import_merged rework is built around.

Pure gemmi/numpy-free stdlib, so it runs unchanged on the slim server.
"""
from __future__ import annotations

import re
from pathlib import Path


# The formats import_merged accepts, plus "unknown".
FORMAT_MTZ = "mtz"
FORMAT_MMCIF = "mmcif"
FORMAT_XDS = "xds"
FORMAT_SHELX = "shelx"
FORMAT_SCALEPACK = "scalepack"
FORMAT_UNKNOWN = "unknown"

_MTZ_MAGIC = b"MTZ "
# A SHELX HKLF record is 3I4 then 2F8.x: the first 28 chars are all in
# [0-9 .+-] with the integer h,k,l in fixed 4-col fields.
_SHELX_RE = re.compile(rb"^[ \d+-]{12}[ \d.+\-eE]{16}")
# A Scalepack merged header: line 1 a small int, line 2 an int (often negative),
# line 3 six floats then a space-group token.
_SCALEPACK_CELL_RE = re.compile(
    rb"^\s*(?:[-+]?\d+\.\d+\s+){6}\S"  # 6 floats then a non-space (SG symbol)
)


def _head_bytes(path, n: int = 4096) -> bytes:
    with open(path, "rb") as fh:
        return fh.read(n)


def detect_format(path) -> str:
    """Classify a reflection file by its content (not its extension)."""
    path = Path(path)
    head = _head_bytes(path)
    if not head:
        return FORMAT_UNKNOWN

    # 1. MTZ — binary, starts with the "MTZ " stamp.
    if head[:4] == _MTZ_MAGIC:
        return FORMAT_MTZ

    # From here treat as text.
    try:
        text = head.decode("utf-8", errors="replace")
    except Exception:
        return FORMAT_UNKNOWN
    lines = text.splitlines()
    if not lines:
        return FORMAT_UNKNOWN

    # 2. XDS_ASCII — a "!FORMAT=XDS_ASCII" (or "!OUTPUT_FILE") banner up top.
    #    This is the crucial disambiguation of the .hkl extension.
    for ln in lines[:5]:
        if "XDS_ASCII" in ln or ln.startswith("!FORMAT") or ln.startswith("!OUTPUT_FILE"):
            return FORMAT_XDS

    # 3. mmCIF — a data_ block that carries reflection categories.
    if any(ln.lstrip().startswith("data_") for ln in lines[:20]):
        if "_refln" in text or "_diffrn_refln" in text:
            return FORMAT_MMCIF

    # 4. Scalepack — 3-line header whose 3rd line is "6 floats + SG symbol".
    if len(lines) >= 3 and _SCALEPACK_CELL_RE.match(lines[2].encode("utf-8", "replace")):
        # guard: the first two lines are short integer-ish headers
        if lines[0].strip().lstrip("-").isdigit():
            return FORMAT_SCALEPACK

    # 5. SHELX — bare fixed-width h k l F sig records from the first line.
    if _SHELX_RE.match(lines[0].encode("utf-8", "replace")):
        return FORMAT_SHELX

    return FORMAT_UNKNOWN


# Which formats gemmi reads natively vs need a pure-Python reader here.
GEMMI_NATIVE = frozenset({FORMAT_MTZ, FORMAT_MMCIF, FORMAT_XDS})
NEEDS_READER = frozenset({FORMAT_SHELX, FORMAT_SCALEPACK})

# Which formats carry cell / space group inside the file (so the resolver need
# not demand them from the user).
CARRIES_CELL_SG = frozenset({FORMAT_MTZ, FORMAT_MMCIF, FORMAT_XDS, FORMAT_SCALEPACK})
