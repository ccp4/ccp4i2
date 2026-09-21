"""
Unit tests for ``lib.superposition`` -- the CA fit behind ``superpose:
method: matrix``.

Pure gemmi, no Django, no CCP4. Structures are synthetic: a chain of CAs on
a slightly irregular helix, so the fit is well conditioned and the
expected transform is known exactly.
"""
import math

import gemmi
import pytest

from ccp4i2.lib import superposition
from ccp4i2.lib.superposition import (
    FIT_RADIUS_MAX,
    FIT_RADIUS_START,
    MIN_FIT_CAS,
    ca_positions,
    fit_files,
    fit_structures,
)


# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------

def _helix_cas(n, start=1):
    """[(seqnum, (x, y, z))] on a helix with a deterministic wobble."""
    out = []
    for i in range(n):
        theta = math.radians(100.0 * i)
        wobble = 0.3 * math.sin(2.1 * i)
        out.append((
            start + i,
            (2.3 * math.cos(theta) + wobble, 2.3 * math.sin(theta), 1.5 * i),
        ))
    return out


def _structure(cas, chain="A", extra_residues=()):
    """A structure whose ALA residues have only a CA each (all a fit needs)."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(100, 100, 100, 90, 90, 90)
    model = gemmi.Model(1)
    ch = gemmi.Chain(chain)
    for seqnum, (x, y, z) in cas:
        res = gemmi.Residue()
        res.name = "ALA"
        res.seqid = gemmi.SeqId(str(seqnum))
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(x, y, z)
        res.add_atom(atom)
        ch.add_residue(res)
    for res in extra_residues:
        ch.add_residue(res)
    model.add_chain(ch)
    st.add_model(model)
    return st


def _rotation(axis, degrees):
    """Row-major 3x3 rotation about ``axis`` (Rodrigues)."""
    ax, ay, az = axis
    norm = math.sqrt(ax * ax + ay * ay + az * az)
    ax, ay, az = ax / norm, ay / norm, az / norm
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    t = 1 - c
    return [
        [t * ax * ax + c, t * ax * ay - s * az, t * ax * az + s * ay],
        [t * ax * ay + s * az, t * ay * ay + c, t * ay * az - s * ax],
        [t * ax * az - s * ay, t * ay * az + s * ax, t * az * az + c],
    ]


def _apply(rot, trans, pos):
    x, y, z = pos
    return tuple(
        rot[r][0] * x + rot[r][1] * y + rot[r][2] * z + trans[r] for r in range(3)
    )


ROT = _rotation((1.0, 2.0, 3.0), 37.0)
TRANS = (3.0, -2.0, 5.0)


def _moved(cas, rot=ROT, trans=TRANS):
    return [(n, _apply(rot, trans, p)) for n, p in cas]


def _transpose(m):
    return [[m[c][r] for c in range(3)] for r in range(3)]


def _assert_inverts(result, rot=ROT, trans=TRANS):
    """The recovered transform is the inverse of the one applied."""
    inv = _transpose(rot)
    expected_vec = _apply(inv, (0, 0, 0), tuple(-t for t in trans))
    for got, want in zip(result.mat, [v for row in inv for v in row]):
        assert got == pytest.approx(want, abs=1e-6)
    for got, want in zip(result.vec, expected_vec):
        assert got == pytest.approx(want, abs=1e-6)


# --------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------

def test_round_trip_recovers_the_inverse():
    # The test that catches the superpose_positions argument-order sign
    # error: a backwards call recovers the forward transform instead.
    cas = _helix_cas(30)
    result = fit_structures(_structure(cas), _structure(_moved(cas)))
    assert result.ok, result.reason
    assert result.rmsd == pytest.approx(0.0, abs=1e-6)
    assert result.atoms == 30
    assert result.radius is None
    _assert_inverts(result)


def test_local_fit_reports_the_radius_it_used():
    cas = _helix_cas(30)
    centre = (0.0, 0.0, 20.0)  # mid-helix; 15 A holds well over MIN_FIT_CAS
    result = fit_structures(_structure(cas), _structure(_moved(cas)), centre)
    assert result.ok, result.reason
    assert result.radius == FIT_RADIUS_START
    assert MIN_FIT_CAS <= result.atoms < 30
    _assert_inverts(result)


def test_disordered_loop_is_tolerated_on_the_intersection():
    cas = _helix_cas(30)
    moving = [(n, p) for n, p in _moved(cas) if not 10 <= n <= 15]
    result = fit_structures(_structure(cas), _structure(moving))
    assert result.ok, result.reason
    assert result.atoms == 24
    assert result.rmsd == pytest.approx(0.0, abs=1e-6)
    _assert_inverts(result)


def test_radius_grows_until_the_count_gate_is_met():
    # Eight CAs inside 15 A of the centre, ten more in the 15-20 A shell:
    # the first sphere is short, the second is not.
    near = [(i + 1, (8.0 * math.cos(i), 8.0 * math.sin(i), 2.0 * i - 7.0)) for i in range(8)]
    shell = [(i + 20, (18.0 * math.cos(0.7 * i), 18.0 * math.sin(0.7 * i), 1.0 * i - 4.5))
             for i in range(10)]
    cas = near + shell
    result = fit_structures(_structure(cas), _structure(_moved(cas)), (0.0, 0.0, 0.0))
    assert result.ok, result.reason
    assert result.radius == FIT_RADIUS_START + 5.0
    assert result.atoms == 18
    _assert_inverts(result)


def test_cap_reached_without_enough_cas_is_a_failure_not_a_bad_fit():
    near = [(i + 1, (8.0 * math.cos(i), 8.0 * math.sin(i), 2.0 * i - 7.0)) for i in range(8)]
    result = fit_structures(_structure(near), _structure(_moved(near)), (0.0, 0.0, 0.0))
    assert not result.ok
    assert result.mat is None and result.vec is None
    assert result.radius == FIT_RADIUS_MAX
    assert result.atoms == 8
    assert "8 shared CA" in result.reason


def test_numbering_mismatch_fails_safely():
    cas = _helix_cas(30)
    renumbered = [(n + 100, p) for n, p in _moved(cas)]
    for centre in (None, (0.0, 0.0, 20.0)):
        result = fit_structures(_structure(cas), _structure(renumbered), centre)
        assert not result.ok
        assert result.atoms == 0
        assert result.mat is None


def test_chain_mismatch_fails_safely():
    cas = _helix_cas(30)
    result = fit_structures(_structure(cas), _structure(_moved(cas), chain="B"))
    assert not result.ok
    assert result.atoms == 0


def test_outliers_are_dropped_and_the_retained_fit_is_tight():
    cas = _helix_cas(30)
    moving = []
    for n, (x, y, z) in _moved(cas):
        if n in (5, 17, 26):
            moving.append((n, (x + 6.0, y - 3.0, z + 2.0)))
        else:
            moving.append((n, (x, y, z)))
    result = fit_structures(_structure(cas), _structure(moving))
    assert result.ok, result.reason
    assert result.atoms == 27
    assert result.rmsd < 0.05
    _assert_inverts(result)


def test_rejection_never_takes_the_set_below_the_gate():
    # Twelve CAs, four of them displaced: dropping them would leave eight,
    # so the first fit stands and its count is reported honestly.
    cas = _helix_cas(MIN_FIT_CAS)
    moving = [
        (n, (x + 6.0, y, z)) if n <= 4 else (n, (x, y, z))
        for n, (x, y, z) in _moved(cas)
    ]
    result = fit_structures(_structure(cas), _structure(moving))
    assert result.ok
    assert result.atoms == MIN_FIT_CAS
    assert result.rmsd > 1.0


# --------------------------------------------------------------------------
# CA selection
# --------------------------------------------------------------------------

def test_calcium_ion_is_not_a_c_alpha():
    ion = gemmi.Residue()
    ion.name = "CA"
    ion.seqid = gemmi.SeqId("501")
    ion.het_flag = "H"
    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("Ca")
    atom.pos = gemmi.Position(1, 1, 1)
    ion.add_atom(atom)
    st = _structure(_helix_cas(3), extra_residues=[ion])
    assert set(ca_positions(st)) == {("A", "1"), ("A", "2"), ("A", "3")}


def test_ligand_without_ca_is_ignored():
    lig = gemmi.Residue()
    lig.name = "DRG"
    lig.seqid = gemmi.SeqId("900")
    lig.het_flag = "H"
    atom = gemmi.Atom()
    atom.name = "C1"
    atom.element = gemmi.Element("C")
    lig.add_atom(atom)
    st = _structure(_helix_cas(2), extra_residues=[lig])
    assert len(ca_positions(st)) == 2


# --------------------------------------------------------------------------
# Scene entry + file wrapper
# --------------------------------------------------------------------------

def test_scene_entry_shape_global_and_local():
    cas = _helix_cas(30)
    ref, mov = _structure(cas), _structure(_moved(cas))

    entry = fit_structures(ref, mov).scene_entry("x0104", "reference")
    assert entry["method"] == "matrix"
    assert entry["move"] == "x0104"
    assert len(entry["mat"]) == 9 and len(entry["vec"]) == 3
    assert entry["fitted"]["onto"] == "reference"
    assert entry["fitted"]["atoms"] == 30
    assert "radius" not in entry["fitted"]
    assert entry["fitted"]["rmsd"] == pytest.approx(0.0, abs=1e-6)

    local = fit_structures(ref, mov, (0.0, 0.0, 20.0)).scene_entry("x0104", "reference")
    assert local["fitted"]["radius"] == FIT_RADIUS_START


def test_failed_fit_has_no_scene_entry():
    cas = _helix_cas(5)
    result = fit_structures(_structure(cas), _structure(_moved(cas)))
    assert not result.ok
    with pytest.raises(ValueError):
        result.scene_entry("a", "b")


def test_fit_files_round_trips_through_pdb(tmp_path):
    cas = _helix_cas(30)
    ref_path = tmp_path / "ref.pdb"
    mov_path = tmp_path / "mov.pdb"
    ref_path.write_text(_structure(cas).make_pdb_string())
    mov_path.write_text(_structure(_moved(cas)).make_pdb_string())
    result = fit_files(ref_path, mov_path)
    assert result.ok, result.reason
    # PDB coordinates carry three decimals, so the recovery is to ~1e-3.
    assert result.rmsd < 5e-3
    assert result.atoms == 30


def test_fit_files_unreadable_is_a_failure_not_an_exception(tmp_path):
    result = fit_files(tmp_path / "missing.pdb", tmp_path / "also_missing.pdb")
    assert not result.ok
    assert "could not read" in result.reason


def test_constants_are_as_designed():
    # The design doc names these; a change should be deliberate.
    assert superposition.FIT_RADIUS_START == 15.0
    assert superposition.FIT_RADIUS_STEP == 5.0
    assert superposition.FIT_RADIUS_MAX == 30.0
    assert superposition.MIN_FIT_CAS == 12
