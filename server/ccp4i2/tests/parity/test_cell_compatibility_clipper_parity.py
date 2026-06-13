"""
Parity test: pure-numpy cells_are_compatible() vs real clipper Cell::equals().

cells_are_compatible() (ccp4i2.core.CCP4XtalData) reimplements clipper's
Cell::equals(other, tol) — "the Frobenius norm of the difference of the two
orthogonalisation matrices is below tol (Angstroms)" — in pure numpy, so that
unit-cell compatibility checks (MTZ merge, sameCrystalAs runtime validity) need
no clipper / CCP4 binary.

This test locks that reimplementation against the genuine clipper result across a
panel of cell pairs and tolerances. It runs only where clipper is importable
(i.e. under ccp4-python); on a slim CCP4-free interpreter it is skipped.

Regression guard: a previous reimplementation used an unrelated per-reciprocal-
axis "mis-indexing resolution" heuristic that was far too lenient — it declared a
doubled cell (160 vs 80 Angstrom) and a 2-degree angle error "compatible". The
cases below would have caught that.
"""

import pytest

clipper = pytest.importorskip("clipper")

from ccp4i2.core.CCP4XtalData import cells_are_compatible


def _clipper_equals(params1, params2, tol):
    """Reference verdict from clipper. Cell_descr takes angles in degrees."""
    c1 = clipper.Cell(clipper.Cell_descr(*params1))
    c2 = clipper.Cell(clipper.Cell_descr(*params2))
    return c1.equals(c2, tol)


# (name, cell1, cell2) — cells as (a, b, c, alpha, beta, gamma)
_BASE = (80.0, 90.0, 100.0, 90.0, 90.0, 90.0)
_TRI = (50.0, 60.0, 70.0, 80.0, 85.0, 95.0)

CELL_PAIRS = [
    ("identical",          _BASE, _BASE),
    ("tiny_length_0.05",   _BASE, (80.05, 90.0, 100.0, 90.0, 90.0, 90.0)),
    ("small_length_0.3",   _BASE, (80.3, 90.3, 100.3, 90.0, 90.0, 90.0)),
    ("medium_length_1.0",  _BASE, (81.0, 91.0, 101.0, 90.0, 90.0, 90.0)),
    ("big_length_5.0",     _BASE, (85.0, 95.0, 105.0, 90.0, 90.0, 90.0)),
    ("doubled_a",          _BASE, (160.0, 90.0, 100.0, 90.0, 90.0, 90.0)),
    ("angle_2deg",         _BASE, (80.0, 90.0, 100.0, 92.0, 90.0, 90.0)),
    ("angle_0.3deg",       _BASE, (80.0, 90.0, 100.0, 90.3, 90.0, 90.0)),
    ("triclinic_same",     _TRI,  _TRI),
    ("triclinic_pert_0.5", _TRI,  (50.5, 60.5, 70.5, 80.0, 85.0, 95.0)),
    ("triclinic_angle",    _TRI,  (50.0, 60.0, 70.0, 81.0, 85.0, 95.0)),
]

TOLERANCES = [0.5, 1.0, 2.0, 3.0]


@pytest.mark.parametrize("tol", TOLERANCES)
@pytest.mark.parametrize("name,p1,p2", CELL_PAIRS, ids=[c[0] for c in CELL_PAIRS])
def test_validity_matches_clipper(name, p1, p2, tol):
    """cells_are_compatible()['validity'] must equal clipper Cell::equals()."""
    reimpl = cells_are_compatible(p1, p2, tolerance=tol)["validity"]
    reference = _clipper_equals(p1, p2, tol)
    assert reimpl == reference, (
        f"{name} @ tol={tol}: reimpl={reimpl} clipper={reference} "
        f"(diff={cells_are_compatible(p1, p2, tolerance=tol)['difference']:.4f})"
    )


def test_symmetry_and_self_identity():
    """Identical cells are always compatible; verdict is order-independent."""
    for _, p1, p2 in CELL_PAIRS:
        assert cells_are_compatible(p1, p1, tolerance=1.0)["validity"] is True
        fwd = cells_are_compatible(p1, p2, tolerance=1.0)["validity"]
        rev = cells_are_compatible(p2, p1, tolerance=1.0)["validity"]
        assert fwd == rev


def test_doubled_cell_is_incompatible():
    """Explicit regression: a doubled cell must NOT be judged compatible."""
    result = cells_are_compatible(_BASE, (160.0, 90.0, 100.0, 90.0, 90.0, 90.0), tolerance=1.0)
    assert result["validity"] is False
    assert _clipper_equals(_BASE, (160.0, 90.0, 100.0, 90.0, 90.0, 90.0), 1.0) is False
