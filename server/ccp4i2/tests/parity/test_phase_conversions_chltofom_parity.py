"""
Parity test: gemmi/numpy-native HL <-> PHIFOM vs the CCP4 chltofom binary.

PhaseDataConverter (ccp4i2.core.conversions.phase_data_converter) converts between
Hendrickson-Lattman coefficients and phase+figure-of-merit using pure gemmi+numpy:
- HL -> PHIFOM: centroid of the HL phase-probability distribution.
- PHIFOM -> HL: unimodal representation with X solving FOM = I1(X)/I0(X).

chltofom is clipper under the hood and computes the same quantities, so it is the
reference. This test runs both on the demo gamma data and asserts close numerical
agreement. It runs only when the chltofom binary is on PATH (i.e. under CCP4);
on a slim CCP4-free interpreter it is skipped.
"""

import shutil
import subprocess

import numpy as np
import pytest

from ccp4i2 import I2_TOP
from ccp4i2.core.conversions.phase_data_converter import PhaseDataConverter

CHLTOFOM = shutil.which("chltofom")
pytestmark = pytest.mark.skipif(CHLTOFOM is None, reason="chltofom binary not on PATH")

HL_INPUT = I2_TOP / "demo_data" / "gamma" / "initial_phases.mtz"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _read(path):
    import gemmi
    return gemmi.read_mtz_file(str(path))


def _col_by_type(mtz, ctype):
    for c in mtz.columns:
        if c.type == ctype:
            return np.array(c, dtype=np.float64)
    raise AssertionError(f"no column of type {ctype} in {[(c.label, c.type) for c in mtz.columns]}")


def _hkl_key(mtz):
    h = np.array(mtz.column_with_label("H"), dtype=np.int64)
    k = np.array(mtz.column_with_label("K"), dtype=np.int64)
    ll = np.array(mtz.column_with_label("L"), dtype=np.int64)
    return list(zip(h.tolist(), k.tolist(), ll.tolist()))


def _aligned(mtz_a, mtz_b, ctype):
    """Return (col_a, col_b) for type `ctype`, aligned on common HKL."""
    ka, kb = _hkl_key(mtz_a), _hkl_key(mtz_b)
    ca, cb = _col_by_type(mtz_a, ctype), _col_by_type(mtz_b, ctype)
    idx_b = {hkl: i for i, hkl in enumerate(kb)}
    rows = [(ca[i], cb[idx_b[hkl]]) for i, hkl in enumerate(ka) if hkl in idx_b]
    assert len(rows) > 0.9 * len(ka), "too few common reflections between outputs"
    a, b = np.array([r[0] for r in rows]), np.array([r[1] for r in rows])
    return a, b


def _circular_diff_deg(a, b):
    d = np.abs(a - b) % 360.0
    return np.minimum(d, 360.0 - d)


def _run_chltofom(mtz_in, mtz_out, colin_flag, colin_cols, colout_cols):
    subprocess.run(
        [CHLTOFOM, "-mtzin", str(mtz_in),
         colin_flag, f"/*/*/[{colin_cols}]",
         "-colout", f"/*/*/[{colout_cols}]",
         "-mtzout", str(mtz_out)],
        check=True, capture_output=True, text=True,
    )


# ---------------------------------------------------------------------------
# HL -> PHIFOM
# ---------------------------------------------------------------------------

def test_hl_to_phifom_matches_chltofom(tmp_path):
    gemmi_out = tmp_path / "gemmi_phifom.mtz"
    chlt_out = tmp_path / "chlt_phifom.mtz"

    PhaseDataConverter._hl_to_phifom(str(HL_INPUT), str(gemmi_out))
    _run_chltofom(HL_INPUT, chlt_out, "-colin-hl", "HLA,HLB,HLC,HLD", "PHI,FOM")

    ga, ca = _read(gemmi_out), _read(chlt_out)

    # FOM agreement (type W) — same algorithm as clipper, so machine-exact
    # (float32 storage ~1e-7).
    fom_g, fom_c = _aligned(ga, ca, "W")
    assert np.max(np.abs(fom_g - fom_c)) < 1e-4, f"max |dFOM|={np.max(np.abs(fom_g - fom_c)):.2e}"

    # Phase agreement (type P), circular, weighted by FOM. Phases of near-zero-FOM
    # reflections are physically undefined, so weight by FOM.
    phi_g, phi_c = _aligned(ga, ca, "P")
    dphi = _circular_diff_deg(phi_g, phi_c)
    w = fom_c
    weighted_mean = float(np.sum(dphi * w) / np.sum(w))
    assert weighted_mean < 0.5, f"FOM-weighted mean |dphi|={weighted_mean:.3f} deg"


# ---------------------------------------------------------------------------
# PHIFOM -> HL  (use a clean chltofom-produced PHIFOM as the common input)
# ---------------------------------------------------------------------------

def test_phifom_to_hl_matches_chltofom(tmp_path):
    phifom = tmp_path / "ref_phifom.mtz"
    _run_chltofom(HL_INPUT, phifom, "-colin-hl", "HLA,HLB,HLC,HLD", "PHI,FOM")

    gemmi_hl = tmp_path / "gemmi_hl.mtz"
    chlt_hl = tmp_path / "chlt_hl.mtz"
    PhaseDataConverter._phifom_to_hl(str(phifom), str(gemmi_hl))
    _run_chltofom(phifom, chlt_hl, "-colin-phifom", "PHI,FOM", "HLA,HLB,HLC,HLD")

    gh, ch = _read(gemmi_hl), _read(chlt_hl)

    # Same algorithm as clipper (atanh/invsim, 0.9999 clamp), so agreement is at
    # the float32 storage level (~1e-4). HLC, HLD are exactly zero.
    kg, kc = _hkl_key(gh), _hkl_key(ch)
    idx_c = {hkl: i for i, hkl in enumerate(kc)}
    common = [(i, idx_c[hkl]) for i, hkl in enumerate(kg) if hkl in idx_c]
    gi = [r[0] for r in common]
    ci = [r[1] for r in common]
    for label in ("HLA", "HLB"):
        g = np.array(gh.column_with_label(label), dtype=np.float64)[gi]
        c = np.array(ch.column_with_label(label), dtype=np.float64)[ci]
        assert np.max(np.abs(g - c)) < 1e-3, f"{label} max|d|={np.max(np.abs(g - c)):.2e}"

    for label in ("HLC", "HLD"):
        assert np.allclose(np.array(gh.column_with_label(label), dtype=np.float64), 0.0, atol=1e-6)
