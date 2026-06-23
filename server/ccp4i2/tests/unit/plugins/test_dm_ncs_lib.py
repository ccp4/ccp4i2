"""Unit tests for dm_multidomain's gemmi-only NCS helpers (no CCP4 binaries).

Guards mask/operator frame consistency: the AHIR reference DARPIN sits OUTSIDE
the primary unit cell (fractional y<0). The averaging mask carries a single
CONTIGUOUS copy of the domain at its real position as a signed-NSTART sub-volume
(write_domain_mask / write_partitioned_masks) -- not a PBC-wrapped image -- and
operator_ref_to_copy maps the real reference onto each copy. Mask and operator
share the crystal frame; if they ever diverged by a lattice vector dm would
report spurious ~0 NCS correlation.
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


def _subvolume_world_points(mask_path, cell):
    """Real-space (crystal-frame) positions of a sub-volume mask's set voxels.

    The mask is a contiguous NSTART sub-volume: voxel (i,j,k) of the written
    block sits at parent grid index (NSTART + (i,j,k)), i.e. fractional
    ((NSTART_x+i)/MX, ...) in the crystal cell. We reconstruct that explicitly
    (gemmi's get_position would use the block's box-local cell, not the crystal
    frame the NSTART/MX header encodes)."""
    mp = gemmi.read_ccp4_map(mask_path)
    arr = np.array(mp.grid, copy=False)
    NS = [mp.header_i32(5), mp.header_i32(6), mp.header_i32(7)]
    M = [mp.header_i32(8), mp.header_i32(9), mp.header_i32(10)]
    idx = np.argwhere(arr > 0)
    out = []
    for i, j, k in idx:
        f = gemmi.Fractional((NS[0] + i) / M[0], (NS[1] + j) / M[1],
                             (NS[2] + k) / M[2])
        p = cell.orthogonalize(f)
        out.append([p.x, p.y, p.z])
    return np.array(out), idx


def test_mask_lands_on_domain_not_lattice_image(model, tmp_path):
    m = model[0]
    cell, sg = model.cell, model.find_spacegroup()
    mask = str(tmp_path / "darpin.msk")
    n = L.write_domain_mask(m, "A", 13, 139, cell, sg, mask, radius=2.5)
    assert n > 1000

    # reference atoms at their REAL position (no lattice shift) -- the mask now
    # carries a contiguous copy there, via signed NSTART
    ref = np.array([[p.x, p.y, p.z]
                    for p in L.all_atom_positions(m, "A", 13, 139)])

    pos, _ = _subvolume_world_points(mask, cell)
    # every mask point must sit within ~ (radius + a grid step) of the domain at
    # its real position, NOT a lattice vector (~100 A) away
    from scipy.spatial import cKDTree
    d = cKDTree(ref).query(pos[:3000])[0]
    assert d.mean() < 3.0, f"mask displaced from domain (mean {d.mean():.1f} A)"

    # the header places it at the domain's REAL cell, which reaches y<0: a
    # negative NSTART proves it is NOT wrapped into the primary cell
    mp = gemmi.read_ccp4_map(mask)
    assert mp.header_i32(6) < 0   # NSTART_y negative (DARPIN reaches outside [0,1))


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

    # Each body is written as its own NSTART sub-volume, so disjointness/union
    # are checked in the shared PARENT grid frame (NSTART + local index), not by
    # overlaying differently-sized blocks.
    def parent_voxels(mask_path):
        mp = gemmi.read_ccp4_map(mask_path)
        arr = np.array(mp.grid, copy=False)
        NS = (mp.header_i32(5), mp.header_i32(6), mp.header_i32(7))
        return {(NS[0] + i, NS[1] + j, NS[2] + k)
                for i, j, k in np.argwhere(arr > 0)}

    # precondition: independent masks DO overlap at the shared boundary
    ma, mb = str(tmp_path / "a.msk"), str(tmp_path / "b.msk")
    L.write_body_mask(m, {"_": "A"}, segs_a, cell, sg, ma, radius=2.5)
    L.write_body_mask(m, {"_": "A"}, segs_b, cell, sg, mb, radius=2.5)
    ga, gb = parent_voxels(ma), parent_voxels(mb)
    assert len(ga & gb) > 0

    # partitioned masks are disjoint, and their union == the independent union
    pa, pb = str(tmp_path / "pa.msk"), str(tmp_path / "pb.msk")
    nsets, n_contested = L.write_partitioned_masks(
        m, {"_": "A"}, [segs_a, segs_b], cell, sg, [pa, pb], radius=2.5)
    assert n_contested > 0
    qa, qb = parent_voxels(pa), parent_voxels(pb)
    assert len(qa & qb) == 0
    assert (qa | qb) == (ga | gb)


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
