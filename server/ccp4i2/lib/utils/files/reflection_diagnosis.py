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

import os
import re
from functools import lru_cache
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
# Unmerged scalepack has a different header: line 1 is "<nsym> <SGname>" and is
# followed by rows of symmetry-operator integers -- quite unlike the merged
# 3-line header. We never READ unmerged scalepack (it is rejected as unmerged),
# but we must still recognise it so the merged-check can reject it rather than
# letting it fall through as "unknown".
_SCALEPACK_UNMERGED_HEAD_RE = re.compile(r"^\s*\d+\s+[A-Za-z][A-Za-z0-9 /\-]*$")


def _looks_unmerged_scalepack(line0: str, line1: str) -> bool:
    if not _SCALEPACK_UNMERGED_HEAD_RE.match(line0.rstrip()):
        return False
    toks = line1.split()
    return len(toks) == 9 and all(t.lstrip("-").isdigit() for t in toks)


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

    # 4. Scalepack (merged) — 3-line header whose 3rd line is "6 floats + SG".
    if len(lines) >= 3 and _SCALEPACK_CELL_RE.match(lines[2].encode("utf-8", "replace")):
        # guard: the first two lines are short integer-ish headers
        if lines[0].strip().lstrip("-").isdigit():
            return FORMAT_SCALEPACK

    # 4b. Scalepack (unmerged) — "<nsym> <SGname>" header then symop rows.
    if len(lines) >= 2 and _looks_unmerged_scalepack(lines[0], lines[1]):
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


# ---------------------------------------------------------------------------
# The unified diagnosis.
# ---------------------------------------------------------------------------
# diagnose_reflection_file() is the single Python authority: it detects the
# format, then reports the cross-format facts the resolver needs — cell, space
# group, wavelength, resolution, whether the data are merged and anomalous, and
# crucially a `needs` list naming any metadata the file does NOT carry and the
# user must therefore supply (cell/SG/data-type for SHELX). gemmi reads
# MTZ/mmCIF/XDS; SHELX/Scalepack are peeked/handled here. gemmi is slim-safe.


def _cell_list(cell) -> list:
    return [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma]


@lru_cache(maxsize=64)
def _diagnose_cached(path_str, mtime, size):
    return diagnose_reflection_file(path_str)


def diagnose_reflection_file_cached(path) -> dict:
    """Like ``diagnose_reflection_file`` but memoised on (path, mtime, size).

    ``validity()`` is polled on every parameter edit, so an uncached read of
    the reflection file each time would be wasteful. The cache key includes
    mtime and size so a replaced file is re-diagnosed. Returns a shallow copy so
    a caller cannot mutate the cached dict.
    """
    p = str(path)
    try:
        st = os.stat(p)
    except OSError:
        return diagnose_reflection_file(p)
    return dict(_diagnose_cached(p, st.st_mtime, st.st_size))


def diagnose_reflection_file(path) -> dict:
    """Detect the format and report the cross-format diagnosis + `needs`."""
    fmt = detect_format(path)
    d = {
        "format": fmt,
        "merged": None,
        "anomalous": None,
        "staraniso": False,   # data from the STARANISO anisotropy server?
        "cell": None,
        "spaceGroup": None,
        "spaceGroupNumber": None,
        "wavelength": None,
        "resolutionHigh": None,
        "resolutionLow": None,
        "needs": [],          # metadata absent from the file, user must supply
        "warnings": [],
    }
    try:
        if fmt == FORMAT_MTZ:
            _diagnose_mtz(path, d)
        elif fmt == FORMAT_XDS:
            _diagnose_xds(path, d)
        elif fmt == FORMAT_MMCIF:
            _diagnose_mmcif(path, d)
        elif fmt == FORMAT_SCALEPACK:
            _diagnose_scalepack(path, d)
        elif fmt == FORMAT_SHELX:
            _diagnose_shelx(path, d)
    except Exception as err:  # never let diagnosis raise on the request path
        d["warnings"].append(f"diagnosis error: {err}")
    return d


def _diagnose_mtz(path, d):
    import gemmi

    mtz = gemmi.read_mtz_file(str(path))
    d["cell"] = _cell_list(mtz.cell)
    if mtz.spacegroup is not None:
        d["spaceGroup"] = mtz.spacegroup.hm
        d["spaceGroupNumber"] = mtz.spacegroup.number
    d["resolutionHigh"] = mtz.resolution_high()
    d["resolutionLow"] = mtz.resolution_low()
    # unmerged iff there are batches or an M/ISYM column; merged otherwise.
    has_msym = any(c.type == "Y" for c in mtz.columns)  # M/ISYM
    d["merged"] = not (len(mtz.batches) > 0 or has_msym)
    # anomalous iff any (+)/(-) column is present.
    d["anomalous"] = any("(+)" in c.label or "(-)" in c.label for c in mtz.columns)
    # StarAniso: the server writes an "SA_flag" column, and stamps "STARANISO"
    # into the MTZ history. Either cue is enough (Qt used the SA_flag column).
    labels_up = [c.label.upper() for c in mtz.columns]
    hist_up = " ".join(mtz.history).upper() if mtz.history else ""
    d["staraniso"] = ("SA_FLAG" in labels_up) or ("STARANISO" in hist_up)
    # first non-base dataset wavelength
    for ds in mtz.datasets:
        if ds.id != 0 and ds.wavelength:
            d["wavelength"] = ds.wavelength
            break


def _diagnose_xds(path, d):
    import gemmi

    xds = gemmi.read_xds_ascii(str(path))
    cc = xds.cell_constants
    if cc and any(cc):
        d["cell"] = list(cc)
    if xds.spacegroup_number:
        d["spaceGroupNumber"] = xds.spacegroup_number
        sg = gemmi.find_spacegroup_by_number(xds.spacegroup_number)
        if sg is not None:
            d["spaceGroup"] = sg.hm
    if xds.wavelength:
        d["wavelength"] = xds.wavelength
    # header MERGE= flag; XdsAscii from CORRECT/INTEGRATE is unmerged.
    head = _head_bytes(path).decode("utf-8", "replace")
    m = re.search(r"MERGE\s*=\s*(TRUE|FALSE)", head)
    d["merged"] = (m.group(1) == "TRUE") if m else False
    fl = re.search(r"FRIEDEL'?S_LAW\s*=\s*(TRUE|FALSE)", head)
    if fl:
        d["anomalous"] = (fl.group(1) == "FALSE")  # Friedel's law FALSE -> anomalous kept
    # If cell/SG missing from an unusual XDS file, flag them.
    if not d["cell"]:
        d["needs"].append("cell")
    if not d["spaceGroupNumber"]:
        d["needs"].append("spaceGroup")


def _diagnose_mmcif(path, d):
    import gemmi
    from ccp4i2.pipelines.import_merged.script import mmcifutils

    doc = gemmi.cif.read(str(path))
    # StarAniso first, independent of the (fragile) per-block extraction below,
    # so the flag is still set if CifBlockInfo chokes on an unusual file.
    d["staraniso"] = _mmcif_names_staraniso(doc)
    rblocks = gemmi.as_refln_blocks(doc)
    if not rblocks:
        d["warnings"].append("no reflection blocks in mmCIF")
        return
    info = mmcifutils.CifBlockInfo(rblocks[0])
    cell = getattr(info, "cell", None)
    if cell and any(cell):
        d["cell"] = list(cell)
    sg = getattr(info, "spacegroupname", None) or getattr(info, "spacegroup", None)
    if sg:
        d["spaceGroup"] = str(sg)
        try:
            d["spaceGroupNumber"] = gemmi.SpaceGroup(str(sg)).number
        except Exception:
            pass
    if getattr(info, "wavelength", None):
        d["wavelength"] = info.wavelength
    d["merged"] = not getattr(info, "unmerged", False)


def _mmcif_names_staraniso(doc) -> bool:
    # StarAniso: named in a _software.name / _computing.* record. Scan those
    # items across all blocks rather than the raw bytes, so an unrelated
    # occurrence of the word elsewhere cannot trip the flag.
    items = ("_software.name", "_software.description",
             "_computing.data_reduction", "_computing.data_scaling")
    for block in doc:
        for item in items:
            for val in block.find_loop(item):
                if "STARANISO" in str(val).upper():
                    return True
            v = block.find_value(item)
            if v and "STARANISO" in str(v).upper():
                return True
    return False


# Cap the merged-case scan; an unmerged file repeats an hkl long before this, so
# the cap only bounds pathologically large *merged* files (assume merged then).
_MAX_HKL_SCAN = 300000


def _scan_hkl_records(path, skip=0, stop_on_zero=False):
    """One pass over fixed-width ``3I4`` hkl records (scalepack/SHELX).

    Returns ``(anomalous, unmerged)``:
      - ``anomalous`` -- any record carries more than two F8 value fields after
        the hkl (the I+/sigma+/I-/sigma- width). Only meaningful for scalepack.
      - ``unmerged`` -- some ``(h, k, l)`` repeats (multiple observations).
        Merged data lists each reflection once (anomalous mates share a row), so
        a repeat means unmerged. Early-exits on the first repeat.
    """
    seen = set()
    anomalous = False
    unmerged = False
    with open(path, "r", errors="replace") as fh:
        for i, rec in enumerate(fh):
            if i < skip:
                continue
            if len(rec) < 12:
                continue
            try:
                h = int(rec[0:4]); k = int(rec[4:8]); l = int(rec[8:12])
            except ValueError:
                continue
            if stop_on_zero and h == 0 and k == 0 and l == 0:
                break  # SHELX end-of-data marker
            if not anomalous and len(rec[12:].rstrip()) > 2 * 8:
                anomalous = True
            key = (h, k, l)
            if key in seen:
                unmerged = True
                break
            seen.add(key)
            if len(seen) >= _MAX_HKL_SCAN:
                break  # enormous unique-hkl set, no repeat -> treat as merged
    return anomalous, unmerged


def _diagnose_scalepack(path, d):
    import gemmi

    # 3-line header; line 3 = 6 cell floats + a space-group symbol.
    with open(path, "r", errors="replace") as fh:
        lines = [next(fh, "") for _ in range(4)]

    # Unmerged scalepack ("<nsym> <SGname>" + symop rows) has no 3-line merged
    # header. We do not read unmerged data (it is rejected as unmerged); just
    # flag it so the merged-check blocks it. The space group is on line 1.
    if len(lines) >= 2 and _looks_unmerged_scalepack(lines[0], lines[1]):
        d["merged"] = False
        toks = lines[0].split()
        if len(toks) >= 2:
            sg_sym = toks[1]
            d["spaceGroup"] = sg_sym
            try:
                sg = gemmi.SpaceGroup(sg_sym)
                d["spaceGroup"] = sg.hm
                d["spaceGroupNumber"] = sg.number
            except Exception:
                pass
        return

    line3 = lines[2]
    m = re.match(r"\s*((?:[-+]?\d+\.\d+\s+){6})(.+)", line3)
    if m:
        d["cell"] = [float(x) for x in m.group(1).split()]
        sg_sym = m.group(2).strip()
        d["spaceGroup"] = sg_sym
        try:
            sg = gemmi.SpaceGroup(sg_sym)
            d["spaceGroup"] = sg.hm
            d["spaceGroupNumber"] = sg.number
        except Exception:
            d["warnings"].append(f"unrecognised space group '{sg_sym}'")
    # Scalepack carries no MERGE flag, so decide merged/unmerged by content:
    # unmerged data repeats (h,k,l) (multiple observations), merged data lists
    # each reflection once (anomalous I+/I- share one row, so they do not
    # repeat). One pass also settles anomalous (a 7-field record). Early-exit on
    # the first duplicate makes the unmerged case cheap; a cap bounds the merged
    # case for very large files.
    anomalous, unmerged = _scan_hkl_records(path, skip=3)
    d["anomalous"] = anomalous
    d["merged"] = not unmerged


def _diagnose_shelx(path, d):
    # SHELX carries no metadata at all — everything must be supplied.
    d["needs"] = ["cell", "spaceGroup", "dataType"]
    # ...but merged/unmerged is still decidable from the records: an unmerged
    # .hkl repeats reflections. (anomalous is not inferable from HKLF width, so
    # it is left undecided here.)
    _, unmerged = _scan_hkl_records(path, skip=0, stop_on_zero=True)
    d["merged"] = not unmerged
