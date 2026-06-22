"""Unit tests for dm_multidomain's gemmi-only NCS helpers (no CCP4 binaries).

Guards the mask/operator frame-consistency bug: the AHIR reference DARPIN sits
OUTSIDE the primary unit cell (fractional y<0), so a naive set_points_around
PBC-wraps the averaging mask ~118 A away from where the operators point, and dm
then reports spurious ~0 NCS correlation. write_domain_mask and
operator_ref_to_copy must bring the reference into the cell with the SAME shift.
"""
import os

import gemmi
import numpy as np
import pytest

import ccp4i2
from ccp4i2.wrappers.dm_multidomain.script import dm_ncs_lib as L

_MODEL = os.path.join(os.path.dirname(ccp4i2.__file__),
                      "demo_data", "ahir", "ahir_model.cif")


@pytest.fixture(scope="module")
def model():
    if not os.path.exists(_MODEL):
        pytest.skip("AHIR demo model not present")
    st = gemmi.read_structure(_MODEL)
    st.setup_entities()
    return st


def test_reference_darpin_is_outside_cell(model):
    """Precondition: the thing that triggered the bug actually holds here."""
    cell = model.cell
    fr = [cell.fractionalize(p)
          for p in L.ca_positions(model[0], "A", 13, 139).values()]
    assert min(f.y for f in fr) < 0.0   # DARPIN reaches outside [0,1)


def test_mask_lands_on_domain_not_lattice_image(model, tmp_path):
    m = model[0]
    cell, sg = model.cell, model.find_spacegroup()
    mask = str(tmp_path / "darpin.msk")
    n = L.write_domain_mask(m, "A", 13, 139, cell, sg, mask, radius=2.5)
    assert n > 1000

    # reference atoms shifted into the cell with the SAME shift the mask used
    s = L._lattice_shift(list(L.ca_positions(m, "A", 13, 139).values()), cell)
    ref = np.array([list(L._shift(p, s, cell))
                    for p in L.all_atom_positions(m, "A", 13, 139)])

    grid = gemmi.read_ccp4_mask(mask).grid
    idx = np.argwhere(np.array(grid, copy=False) > 0)[:3000]
    pos = np.array([list(grid.get_position(int(i), int(j), int(k)))
                    for i, j, k in idx])
    # every mask point must sit within ~ (radius + a grid step) of the domain,
    # NOT a lattice vector (~100 A) away
    from scipy.spatial import cKDTree
    d = cKDTree(ref).query(pos)[0]
    assert d.mean() < 3.0, f"mask displaced from domain (mean {d.mean():.1f} A)"


def test_operator_targets_intended_copy(model):
    """Each A->copy operator must superpose A's domain onto THAT copy."""
    m = model[0]
    cell = model.cell
    ref_ca = L.ca_positions(m, "A", 13, 139)
    for c in ["B", "C", "D", "E", "F"]:
        T, rmsd = L.operator_ref_to_copy(m, "A", c, 13, 139, cell=cell)
        assert rmsd < 1.5   # genuine superposition onto that copy


# -- role / segment parsing (pure, no model) --------------------------------

def test_parse_segments_implicit_and_roled():
    assert L.parse_segments("340-485") == [("_", 340, 485)]
    assert L.parse_segments("cyclin:10-95, CDK:45-60") == [
        ("cyclin", 10, 95), ("CDK", 45, 60)]


def test_parse_assembly_rows_bare_roled_and_partial():
    # homomer: bare chains, implicit role
    assert L.parse_assembly_rows(["A", "B", "C"]) == [
        {"_": "A"}, {"_": "B"}, {"_": "C"}]
    # heteromer with a partial (AmBn) final instance
    assert L.parse_assembly_rows(
        ["CDK=A cyclin=B", "CDK=C cyclin=D", "CDK=E"]) == [
        {"CDK": "A", "cyclin": "B"}, {"CDK": "C", "cyclin": "D"}, {"CDK": "E"}]


def test_body_name_and_instance_label():
    assert L.body_name([("_", 340, 485)]) == "340_485"
    assert L.body_name([("cyclin", 10, 95), ("CDK", 45, 60)]) == \
        "cyclin_10_95__CDK_45_60"
    assert L.instance_label({"CDK": "A", "cyclin": "B"}) == "A+B"


# -- multi-segment / cross-role correspondence (uses the model) -------------

def test_multisegment_body_concatenates_across_ranges(model):
    """A body spanning two halves of ONE rigid domain (340-410 + 411-485)
    superposes onto B as a single rigid unit -- the concatenated correspondence
    must fit as well as superposing the whole 340-485 in one go.

    (Fusing two genuinely DIFFERENT domains, e.g. 13-139 + 340-485, would fit
    poorly -- that they move independently is exactly why this task exists.)
    """
    m = model[0]
    cell = model.cell
    segs = [("_", 340, 410), ("_", 411, 485)]
    T, rmsd = L.operator_ref_to_copy_body(
        m, {"_": "A"}, {"_": "B"}, segs, cell=cell)
    _, rmsd_whole = L.operator_ref_to_copy_body(
        m, {"_": "A"}, {"_": "B"}, [("_", 340, 485)], cell=cell)
    assert rmsd < 1.5
    assert abs(rmsd - rmsd_whole) < 1e-6   # split == whole for a rigid domain


def test_cross_role_keying_matches_multisegment(model):
    """Routing the same atoms through two distinct roles (X<-A, Y<-A) must give
    the same operator as the single-role multi-segment body -- proving the
    (role, seqid) key gathers cross-chain atoms without colliding on seqid."""
    m = model[0]
    cell = model.cell
    ref = {"X": "A", "Y": "A"}
    copy = {"X": "B", "Y": "B"}
    segs = [("X", 340, 410), ("Y", 411, 485)]
    T, rmsd = L.operator_ref_to_copy_body(m, ref, copy, segs, cell=cell)
    _, rmsd_single = L.operator_ref_to_copy_body(
        m, {"_": "A"}, {"_": "B"}, [("_", 340, 410), ("_", 411, 485)], cell=cell)
    assert abs(rmsd - rmsd_single) < 1e-6


def test_partitioned_masks_are_disjoint(model, tmp_path):
    """Two abutting bodies (140-339 and 340-485) overlap when masked
    independently, but write_partitioned_masks resolves them into disjoint masks
    (dm requires non-overlapping NCS masks) without losing any claimed volume."""
    m = model[0]
    cell, sg = model.cell, model.find_spacegroup()
    segs_a = [("_", 140, 339)]
    segs_b = [("_", 340, 485)]

    # precondition: independent masks DO overlap at the shared boundary
    ma, mb = str(tmp_path / "a.msk"), str(tmp_path / "b.msk")
    L.write_body_mask(m, {"_": "A"}, segs_a, cell, sg, ma, radius=2.5)
    L.write_body_mask(m, {"_": "A"}, segs_b, cell, sg, mb, radius=2.5)
    ga = np.array(gemmi.read_ccp4_mask(ma).grid, copy=False) > 0
    gb = np.array(gemmi.read_ccp4_mask(mb).grid, copy=False) > 0
    assert int(np.count_nonzero(ga & gb)) > 0

    # partitioned masks are disjoint, and their union == the independent union
    pa, pb = str(tmp_path / "pa.msk"), str(tmp_path / "pb.msk")
    nsets, n_contested = L.write_partitioned_masks(
        m, {"_": "A"}, [segs_a, segs_b], cell, sg, [pa, pb], radius=2.5)
    assert n_contested > 0
    qa = np.array(gemmi.read_ccp4_mask(pa).grid, copy=False) > 0
    qb = np.array(gemmi.read_ccp4_mask(pb).grid, copy=False) > 0
    assert int(np.count_nonzero(qa & qb)) == 0
    assert int(np.count_nonzero(qa | qb)) == int(np.count_nonzero(ga | gb))


def test_body_operators_skips_partial_instance(model):
    """A copy instance missing one of the body's roles is skipped (AmBn): only
    identity remains when the sole copy lacks a required role."""
    m = model[0]
    cell = model.cell
    instances = [{"X": "A", "Y": "A"}, {"X": "B"}]   # copy lacks role Y
    segs = [("X", 13, 139), ("Y", 340, 485)]
    ops, rmsds = L.body_operators(m, instances, segs, cell=cell)
    assert ops == [L.IDENTITY_DM]   # identity only; partial copy skipped
    assert rmsds == {}
