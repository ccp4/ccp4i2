"""Pure-Python/gemmi readers for reflection formats gemmi cannot read natively.

The slim server has no CCP4 binaries, so the classic converters ``f2mtz`` and
``scalepack2mtz`` are unavailable on any path the slim server must run. gemmi
0.7.5 reads MTZ, mmCIF and XDS_ASCII natively; SHELX (``.hkl``) and Scalepack
(``.sca``) it does not. These readers fill that gap in pure Python + gemmi,
building a :class:`gemmi.Mtz` from the parsed data, and are pinned against the
binaries they replace by ``tests/parity/`` (the same discipline as the
``freerflag``/``matthews``/``chltofom`` ports).

Both formats are **fixed-width Fortran** — fields can abut with no separating
space when a value fills its column — so parsing is by *column position*, never
``str.split()`` (verified: ``scalepack2mtz`` reads a crafted abutting line
correctly, which whitespace-splitting cannot).
"""
from __future__ import annotations

import gemmi
import numpy as np


def _as_cell(cell) -> gemmi.UnitCell:
    if isinstance(cell, gemmi.UnitCell):
        return cell
    a, b, c, al, be, ga = (float(x) for x in cell)
    return gemmi.UnitCell(a, b, c, al, be, ga)


def _as_spacegroup(sg) -> gemmi.SpaceGroup:
    if isinstance(sg, gemmi.SpaceGroup):
        return sg
    if isinstance(sg, int):
        return gemmi.find_spacegroup_by_number(sg)
    return gemmi.SpaceGroup(str(sg))


def _slice_float(line: str, start: int, width: int):
    """Fixed-width float field; blank/short field -> None (missing)."""
    field = line[start:start + width]
    if not field.strip():
        return None
    return float(field)


def _slice_int(line: str, start: int, width: int) -> int:
    return int(line[start:start + width])


# --------------------------------------------------------------------------
# Scalepack (.sca) — merged output of the DENZO/HKL Scalepack program.
# --------------------------------------------------------------------------
# Format (CCP4 scalepack2mtz.html + empirical verification):
#   3-line header; line 3 carries the cell (and a lattice/point-group string).
#   Reflection records, fixed width 3I4 then F8.1 fields:
#     non-anomalous (5 cols): H K L  I SIGI
#     anomalous     (7 cols): H K L  I(+) SIGI(+) I(-) SIGI(-)
# scalepack2mtz is invoked with CELL/SYMMETRY keywords, so cell + space group
# are supplied here (not trusted from the header), matching the reference.

_SCA_HKL_W = 4          # I4 for each of H, K, L
_SCA_F_W = 8            # F8.1 for each intensity/sigma field


def read_scalepack(path, cell, spacegroup, anomalous: bool = True, scff: float = 1.0) -> gemmi.Mtz:
    """Read a merged Scalepack ``.sca`` file into a gemmi.Mtz.

    A line-for-line port of ``scalepack2mtz`` (Fortran source, Winn/CCP4). Each
    record is read ``(3I4, 4F8.0)`` — H K L then FPLUS SIGFPLUS FNEG SIGFNEG,
    with blank fields reading as 0.0. A **missing** mate is one whose sigma is
    ``<= 0`` (a blank sigma field reads as 0). hkl are written straight through
    (no ASU reduction — scalepack has already reduced them). Output:

      anomalous  (``KANOM != 0``): H K L IMEAN SIGIMEAN I(+) SIGI(+) I(-) SIGI(-)
      merged     (``KANOM == 0``): H K L I SIGI

    with (both mates present)::

        SIGIMEAN = SIGFPLUS*SIGFNEG / (SIGFPLUS + SIGFNEG)
        IMEAN    = (FPLUS/SIGFPLUS + FNEG/SIGFNEG) * SIGIMEAN

    ``scff`` mirrors the SCALE keyword (default 1.0).
    """
    cell = _as_cell(cell)
    sg = _as_spacegroup(spacegroup)

    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    if len(lines) < 4:
        raise ValueError(f"{path}: too short to be a Scalepack file")
    records = lines[3:]  # skip the 3-line header (cell/SG supplied as args)

    NAN = float("nan")

    def _f(rec, i):
        """i-th F8.0 field after HKL; blank reads as 0.0 (Fortran F8.0)."""
        field = rec[3 * _SCA_HKL_W + i * _SCA_F_W: 3 * _SCA_HKL_W + (i + 1) * _SCA_F_W]
        field = field.strip()
        return float(field) if field else 0.0

    rows = []
    for rec in records:
        if not rec.strip():
            continue
        try:
            h = _slice_int(rec, 0, _SCA_HKL_W)
            k = _slice_int(rec, _SCA_HKL_W, _SCA_HKL_W)
            l = _slice_int(rec, 2 * _SCA_HKL_W, _SCA_HKL_W)
        except ValueError:
            continue
        # terminators / bad indices (source: |index| >= 999 skipped)
        if abs(h) >= 999 or abs(k) >= 999 or abs(l) >= 999:
            continue
        fp, sfp, fn, sfn = _f(rec, 0), _f(rec, 1), _f(rec, 2), _f(rec, 3)

        if not anomalous:
            if sfp > 0.0:
                rows.append([h, k, l, scff * fp, scff * sfp])
            else:
                rows.append([h, k, l, NAN, NAN])
            continue

        # anomalous output: IMEAN SIGIMEAN I(+) SIGI(+) I(-) SIGI(-)
        if sfp <= 0.0 and sfn <= 0.0:
            rows.append([h, k, l, NAN, NAN, NAN, NAN, NAN, NAN])
        elif sfn <= 0.0:  # only I(+)
            rows.append([h, k, l, scff * fp, scff * sfp, scff * fp, scff * sfp, NAN, NAN])
        elif sfp <= 0.0:  # only I(-)
            rows.append([h, k, l, scff * fn, scff * sfn, NAN, NAN, scff * fn, scff * sfn])
        else:             # both mates
            sigimean = sfp * sfn / (sfp + sfn)
            imean = (fp / sfp + fn / sfn) * sigimean
            rows.append([h, k, l, scff * imean, scff * sigimean,
                         scff * fp, scff * sfp, scff * fn, scff * sfn])

    if anomalous:
        labels = ["H", "K", "L", "IMEAN", "SIGIMEAN", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
        types = ["H", "H", "H", "J", "Q", "K", "M", "K", "M"]
    else:
        labels = ["H", "K", "L", "IMEAN", "SIGIMEAN"]
        types = ["H", "H", "H", "J", "Q"]

    return _build_mtz(cell, sg, labels, types, rows)


# --------------------------------------------------------------------------
# SHELX (.hkl) — HKLF 4 (intensities) fixed format 3I4,2F8.2 (+ optional I4 batch)
# --------------------------------------------------------------------------
_SHELX_HKL_W = 4
_SHELX_F_W = 8


def read_shelx(path, cell, spacegroup, intensities: bool = True) -> gemmi.Mtz:
    """Read a SHELX ``.hkl`` (HKLF 4 = intensities, or HKLF 3 = amplitudes).

    Fixed format ``3I4,2F8.2`` then an optional ``I4`` batch/scale number.
    SHELX carries no cell or space group, so both are supplied. A trailing
    ``0 0 0`` record terminates the data (SHELX convention). Output columns are
    ``H K L I SIGI`` (intensities) or ``H K L F SIGF`` (amplitudes).
    """
    cell = _as_cell(cell)
    sg = _as_spacegroup(spacegroup)

    rows = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if len(line.rstrip("\n")) < 3 * _SHELX_HKL_W + 2 * _SHELX_F_W:
                # too short for a full h k l val sig record
                if not line.strip():
                    continue
            try:
                h = _slice_int(line, 0, _SHELX_HKL_W)
                k = _slice_int(line, _SHELX_HKL_W, _SHELX_HKL_W)
                l = _slice_int(line, 2 * _SHELX_HKL_W, _SHELX_HKL_W)
            except ValueError:
                continue
            if h == 0 and k == 0 and l == 0:
                break  # SHELX end-of-data marker
            base = 3 * _SHELX_HKL_W
            val = _slice_float(line, base, _SHELX_F_W)
            sig = _slice_float(line, base + _SHELX_F_W, _SHELX_F_W)
            if val is None or sig is None:
                continue
            rows.append([h, k, l, val, sig])

    if intensities:
        labels, types = ["H", "K", "L", "I", "SIGI"], ["H", "H", "H", "J", "Q"]
    else:
        labels, types = ["H", "K", "L", "F", "SIGF"], ["H", "H", "H", "F", "Q"]
    return _build_mtz(cell, sg, labels, types, rows)


# --------------------------------------------------------------------------
def _build_mtz(cell, sg, labels, types, rows) -> gemmi.Mtz:
    mtz = gemmi.Mtz(with_base=False)
    mtz.cell = cell
    mtz.spacegroup = sg
    ds = mtz.add_dataset("HKL_base")
    ds.cell = cell
    base = mtz.add_dataset("crystal")
    base.cell = cell
    for label, ctype in zip(labels, types):
        col = mtz.add_column(label, ctype)
        col.dataset_id = 0 if label in ("H", "K", "L") else 1
    arr = np.array(rows, dtype=np.float32) if rows else np.empty((0, len(labels)), np.float32)
    mtz.set_data(arr)
    return mtz
