"""Parity test: the pure-Python SHELX reader vs the CCP4 f2mtz binary.

reflection_formats.read_shelx reproduces f2mtz's SHELX HKLF-4 conversion so a
.hkl can be imported on a CCP4-free (slim) server, where f2mtz is unavailable.

There is no real-world SHELX .hkl in the tree (every .hkl here is XDS_ASCII), so
this test *generates* a valid HKLF-4 file — deliberately including records where
a filled field abuts the next with no separating space (the classic fixed-width
footgun) — and pins the reader against f2mtz reading the same file. That both
agree on the abut records is the point: whitespace-splitting cannot read them.

Runs only where f2mtz is on PATH (under CCP4); skipped on a slim interpreter.
"""
import shutil
import subprocess

import pytest
import gemmi
import numpy as np

from ccp4i2.lib.utils.files.reflection_formats import read_shelx

F2MTZ = shutil.which("f2mtz")
pytestmark = pytest.mark.skipif(F2MTZ is None, reason="f2mtz binary not on PATH")

CELL = (50.0, 50.0, 50.0, 90.0, 90.0, 90.0)

# (h, k, l, I, SIGI). Values chosen to fill F8.2 without overflowing
# (|F| <= 99999.99, negatives >= -9999.99) and to force abutting fields.
RECORDS = [
    (1, 2, 3, 1000.00, 10.00),
    (4, 5, 6, 12345.67, 88.20),
    (-10, -20, -30, -1234.56, 99.99),   # negatives
    (5, 6, -100, -9999.99, 99.99),      # l=-100 fills I4, I fills F8.2 -> abut
    (-100, 7, 8, 99999.99, 12.34),      # h=-100 fills I4
]


def _hklf4(h, k, l, f, s):
    line = "%4d%4d%4d%8.2f%8.2f" % (h, k, l, f, s)
    assert len(line) == 28, (line, len(line))  # 3I4 + 2F8.2, no overflow
    return line


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    d = tmp_path_factory.mktemp("shelx")
    hkl = d / "gen.hkl"
    with hkl.open("w") as fh:
        for r in RECORDS:
            fh.write(_hklf4(*r) + "\n")
        fh.write(_hklf4(0, 0, 0, 0, 0) + "\n")  # SHELX end-of-data marker
    out = d / "ref.mtz"
    script = (
        "TITLE shelx-parity\n"
        f"CELL {' '.join(str(x) for x in CELL)}\n"
        "SYMMETRY 1\n"
        "LABOUT H K L I SIGI\n"
        "CTYPOUT H H H J Q\n"
        "FORMAT '(3I4,2F8.2)'\n"
        "END\n"
    )
    subprocess.run(
        [F2MTZ, "HKLIN", str(hkl), "HKLOUT", str(out)],
        input=script, text=True, capture_output=True, check=True,
    )
    return hkl, gemmi.read_mtz_file(str(out))


def _by_hkl(mtz):
    arr = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]
    return {(int(r[0]), int(r[1]), int(r[2])): dict(zip(labels, r)) for r in arr}


def test_abut_records_defeat_whitespace_split():
    # Guard: at least one record must actually abut, or the test proves nothing.
    abut = [r for r in RECORDS if len(_hklf4(*r).split()) != 5]
    assert abut, "no abutting record present — parity test would be vacuous"


def test_shelx_reader_matches_f2mtz(generated):
    hkl, ref_mtz = generated
    mine = read_shelx(hkl, CELL, 1, intensities=True)

    got = _by_hkl(mine)
    ref = _by_hkl(ref_mtz)
    assert set(got) == set(ref), (
        f"hkl set differs: only-mine={set(got)-set(ref)} only-ref={set(ref)-set(got)}"
    )
    for hkl_key, rrow in ref.items():
        grow = got[hkl_key]
        for lab in ("I", "SIGI"):
            assert abs(grow[lab] - rrow[lab]) <= 1e-3, (
                f"{hkl_key} {lab}: {grow[lab]} vs {rrow[lab]}"
            )
