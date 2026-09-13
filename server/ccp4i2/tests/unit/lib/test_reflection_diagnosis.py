"""Unit tests for content-based reflection-format detection (CCP4-free).

detect_format must classify by *content*, not extension — most importantly it
must call a ``.hkl`` that is actually XDS_ASCII "xds", not "shelx", and must
recognise Scalepack's 3-line header. Pure stdlib, so this runs on the slim
interpreter and in CI.
"""
import os

import pytest

from ccp4i2 import I2_TOP
from ccp4i2.lib.utils.files.reflection_diagnosis import (
    detect_format,
    FORMAT_MTZ,
    FORMAT_MMCIF,
    FORMAT_XDS,
    FORMAT_SHELX,
    FORMAT_SCALEPACK,
    FORMAT_UNKNOWN,
)

DEMO = I2_TOP / "demo_data"


# --- synthetic content (self-contained, format-defining) ----------------------

SHELX_HKLF4 = "   1   2   3 1000.00   10.00\n   0   0   0    0.00    0.00\n"

SCALEPACK = (
    "    1\n"
    " -987\n"
    "    83.090    96.790    57.950    90.000    90.000    90.000 c 2 2 21\n"
    "   0   0   2   887.6    17.3\n"
)

XDS = (
    "!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=FALSE\n"
    "!OUTPUT_FILE=XDS_ASCII.HKL\n"
    "   -47   -11    -2 -6.014E-01  7.599E+00\n"
)

MMCIF_SF = (
    "data_r1abcsf\n"
    "loop_\n"
    "_refln.index_h\n_refln.index_k\n_refln.index_l\n_refln.F_meas_au\n"
    "0 0 2 123.4\n"
)


@pytest.mark.parametrize(
    "text,expected",
    [
        (SHELX_HKLF4, FORMAT_SHELX),
        (SCALEPACK, FORMAT_SCALEPACK),
        (XDS, FORMAT_XDS),
        (MMCIF_SF, FORMAT_MMCIF),
        ("garbage not a reflection file\nsecond line\n", FORMAT_UNKNOWN),
        ("", FORMAT_UNKNOWN),
    ],
)
def test_detect_format_synthetic(tmp_path, text, expected):
    p = tmp_path / "sample.dat"
    p.write_text(text)
    assert detect_format(p) == expected


def test_mtz_magic(tmp_path):
    p = tmp_path / "x.mtz"
    p.write_bytes(b"MTZ \x00\x00\x00\x00rest of a binary mtz")
    assert detect_format(p) == FORMAT_MTZ


# --- real demo files (prove it on the wild data) ------------------------------

_DEMO_CASES = {
    FORMAT_MTZ: DEMO / "gamma" / "freeR.mtz",
    FORMAT_SCALEPACK: (
        DEMO / "baz2b" / "BAZ2BA_x828.xia2" / "3daii-run" / "DataFiles"
        / "nt5073v16_xBAZ2BAx8281_scaled.sca"
    ),
    # The decisive case: a .hkl that is really XDS_ASCII must NOT read as shelx.
    FORMAT_XDS: DEMO / "ceue" / "apo-ceue-sad-sweep1.hkl",
}


@pytest.mark.parametrize("expected,path", list(_DEMO_CASES.items()), ids=lambda v: str(v))
def test_detect_format_demo(expected, path):
    if not isinstance(path, str) and not os.path.exists(path):
        pytest.skip(f"demo file absent: {path}")
    assert detect_format(path) == expected
