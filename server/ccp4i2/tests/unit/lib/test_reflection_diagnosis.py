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
    diagnose_reflection_file,
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


# --- diagnose_reflection_file (needs gemmi; CCP4-free) ------------------------


def test_diagnose_shelx_reports_needs(tmp_path):
    p = tmp_path / "s.hkl"
    p.write_text(SHELX_HKLF4)
    d = diagnose_reflection_file(p)
    assert d["format"] == FORMAT_SHELX
    # SHELX carries no metadata: the resolver must demand it.
    assert d["needs"] == ["cell", "spaceGroup", "dataType"]
    assert d["merged"] is True
    assert d["cell"] is None and d["spaceGroup"] is None


@pytest.mark.skipif(
    not (DEMO / "gamma" / "freeR.mtz").exists(), reason="demo mtz absent"
)
def test_diagnose_mtz_demo():
    d = diagnose_reflection_file(DEMO / "gamma" / "freeR.mtz")
    assert d["format"] == FORMAT_MTZ
    assert d["merged"] is True         # a merged MTZ, not the getMerged() stub
    assert d["needs"] == []
    assert d["cell"] and len(d["cell"]) == 6
    assert d["spaceGroupNumber"] == 19  # P 21 21 21


_SCA_DEMO = (
    DEMO / "baz2b" / "BAZ2BA_x828.xia2" / "3daii-run" / "DataFiles"
    / "nt5073v16_xBAZ2BAx8281_scaled.sca"
)


@pytest.mark.skipif(not _SCA_DEMO.exists(), reason="demo sca absent")
def test_diagnose_scalepack_demo():
    d = diagnose_reflection_file(_SCA_DEMO)
    assert d["format"] == FORMAT_SCALEPACK
    assert d["merged"] is True
    assert d["anomalous"] is True       # baz2b .sca is anomalous
    assert d["spaceGroupNumber"] == 20  # C 2 2 21
    assert d["needs"] == []


_XDS_DEMO = DEMO / "ceue" / "apo-ceue-sad-sweep1.hkl"


@pytest.mark.skipif(not _XDS_DEMO.exists(), reason="demo xds absent")
def test_diagnose_xds_demo_is_unmerged():
    d = diagnose_reflection_file(_XDS_DEMO)
    assert d["format"] == FORMAT_XDS
    # The old getMerged() stub always returned True; real detection sees MERGE=FALSE.
    assert d["merged"] is False
    assert d["anomalous"] is True       # FRIEDEL'S_LAW=FALSE


# --- StarAniso detection (CCP4-free; the signal the Qt GUI had and React dropped) ---


def _write_mtz(path, extra_cols):
    """Minimal single-reflection MTZ with H,K,L,F,SIGF plus `extra_cols`
    (list of (label, type)); values default to 0. gemmi only, no CCP4."""
    import gemmi
    import numpy as np

    mtz = gemmi.Mtz(with_base=True)
    mtz.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    mtz.add_dataset("data")
    cols = [("F", "F"), ("SIGF", "Q")] + list(extra_cols)
    for label, ctype in cols:
        mtz.add_column(label, ctype)
    ncol = 3 + len(cols)  # H K L + the rest
    row = [0, 0, 1] + [0.0] * len(cols)
    mtz.set_data(np.array([row], dtype=np.float32).reshape(1, ncol))
    mtz.write_to_file(str(path))


def test_diagnose_mtz_staraniso_column(tmp_path):
    p = tmp_path / "sa.mtz"
    _write_mtz(p, [("SA_flag", "I")])
    d = diagnose_reflection_file(p)
    assert d["format"] == FORMAT_MTZ
    assert d["staraniso"] is True


def test_diagnose_mtz_no_staraniso(tmp_path):
    p = tmp_path / "plain.mtz"
    _write_mtz(p, [("FreeR_flag", "I")])
    d = diagnose_reflection_file(p)
    assert d["format"] == FORMAT_MTZ
    assert d["staraniso"] is False


MMCIF_STARANISO = (
    "data_r1abcsf\n"
    "_software.name STARANISO\n"
    "loop_\n"
    "_refln.index_h\n_refln.index_k\n_refln.index_l\n_refln.F_meas_au\n"
    "0 0 2 123.4\n"
)


def test_diagnose_mmcif_staraniso(tmp_path):
    p = tmp_path / "sa.cif"
    p.write_text(MMCIF_STARANISO)
    d = diagnose_reflection_file(p)
    assert d["format"] == FORMAT_MMCIF
    assert d["staraniso"] is True


def test_diagnose_mmcif_no_staraniso(tmp_path):
    p = tmp_path / "plain.cif"
    p.write_text(MMCIF_SF)   # the plain sfCIF fixture above, no _software.name
    d = diagnose_reflection_file(p)
    assert d["format"] == FORMAT_MMCIF
    assert d["staraniso"] is False
