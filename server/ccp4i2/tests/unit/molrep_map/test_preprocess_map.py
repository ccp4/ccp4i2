"""Unit tests for the gemmi map surgery behind ``molrep_map``.

CCP4-free: needs only gemmi (a slim-server dependency), no CCP4 binaries and no
coot/chapi. These lock the invariants the design leans on -- an origin-inverting
hand flip, a same-cell down-sample, and a trim that keeps the full-cell sampling
and encodes its position via ``nxstart`` (what Moorhen/clipper read).
"""

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.wrappers.molrep_map.script import preprocess_map as pp


# --- fixtures ---------------------------------------------------------------

def _blob_map(n=72, cell_len=72.0, centre=(36, 36, 40), sigma2=60.0):
    cell = gemmi.UnitCell(cell_len, cell_len, cell_len, 90, 90, 90)
    grid = gemmi.FloatGrid(n, n, n)
    grid.set_unit_cell(cell)
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    arr = np.array(grid, copy=False)
    zz, yy, xx = np.mgrid[0:n, 0:n, 0:n]
    cx, cy, cz = centre
    arr[:] = np.exp(-(((xx - cx) ** 2 + (yy - cy) ** 2
                       + (zz - cz) ** 2) / sigma2)).astype("float32")
    m = gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header()
    m.setup(float("nan"), gemmi.MapSetup.Full)
    return m


def _three_atom_pdb(tmp_path, coords):
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(72, 72, 72, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    res = gemmi.Residue()
    res.name = "ALA"
    res.seqid = gemmi.SeqId("1")
    for i, (x, y, z) in enumerate(coords):
        at = gemmi.Atom()
        at.name = f"C{i}"
        at.element = gemmi.Element("C")
        at.pos = gemmi.Position(x, y, z)
        res.add_atom(at)
    chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    path = tmp_path / "model.pdb"
    st.write_pdb(str(path))
    return path


def _reread(m, tmp_path, name="out.map"):
    """Write and read back a map without re-expanding a partial map."""
    path = tmp_path / name
    m.write_ccp4_map(str(path))
    out = gemmi.read_ccp4_map(str(path))
    out.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
    return out


# --- composite_size ---------------------------------------------------------

@pytest.mark.parametrize("n,expected", [(1, 1), (35, 35), (70, 70), (37, 40),
                                        (101, 105), (127, 128)])
def test_composite_size_is_small_prime_product(n, expected):
    assert pp.composite_size(n) == expected
    # and it factorises over {2,3,5,7}
    m = pp.composite_size(n)
    for p in (2, 3, 5, 7):
        while m % p == 0:
            m //= p
    assert m == 1


# --- flip_hand --------------------------------------------------------------

def test_flip_hand_preserves_shape_and_cell():
    m = _blob_map()
    f = pp.flip_hand(m)
    assert (f.grid.nu, f.grid.nv, f.grid.nw) == (m.grid.nu, m.grid.nv, m.grid.nw)
    assert f.grid.unit_cell.parameters == m.grid.unit_cell.parameters


def test_flip_hand_is_an_involution():
    """Two flips return the original density exactly (point inversion is self-inverse)."""
    m = _blob_map()
    twice = pp.flip_hand(pp.flip_hand(m))
    assert np.allclose(np.array(twice.grid, copy=False),
                       np.array(m.grid, copy=False))


def test_flip_hand_inverts_through_origin():
    """A blob off-centre lands at the inverted position V'[i]=V[(N-i)%N]."""
    m = _blob_map(n=64, cell_len=64.0, centre=(20, 24, 28))
    f = pp.flip_hand(m)
    a = np.array(m.grid, copy=False)
    b = np.array(f.grid, copy=False)
    i, j, k = np.unravel_index(np.argmax(a), a.shape)
    N = a.shape
    expect = ((N[0] - i) % N[0], (N[1] - j) % N[1], (N[2] - k) % N[2])
    got = np.unravel_index(np.argmax(b), b.shape)
    assert got == expect


# --- downsample -------------------------------------------------------------

def test_downsample_coarsens_grid_same_cell():
    m = _blob_map(n=72, cell_len=72.0)
    d = pp.downsample(m, 2)
    assert d.grid.nu < m.grid.nu
    assert d.grid.unit_cell.parameters == m.grid.unit_cell.parameters  # same cell
    # peak survives a mild down-sample
    assert np.array(d.grid, copy=False).max() > 0.5


def test_downsample_factor_one_is_noop():
    m = _blob_map()
    assert pp.downsample(m, 1) is m


# --- blur -------------------------------------------------------------------

def test_blur_preserves_dims_and_reduces_peak():
    m = _blob_map()
    b = pp.blur(m, 80)
    assert (b.grid.nu, b.grid.nv, b.grid.nw) == (m.grid.nu, m.grid.nv, m.grid.nw)
    # blurring spreads density: the peak drops, the integral is ~conserved
    assert np.array(b.grid, copy=False).max() < np.array(m.grid, copy=False).max()
    assert np.isclose(np.array(b.grid, copy=False).sum(),
                      np.array(m.grid, copy=False).sum(), rtol=1e-3)


# --- trim (crop) ------------------------------------------------------------

def test_crop_keeps_full_sampling_and_spacing(tmp_path):
    m = _blob_map()
    pdb = _three_atom_pdb(tmp_path, [(30, 30, 34), (42, 42, 46), (36, 38, 40)])
    box = pp.model_frac_box(pdb, m.grid.unit_cell, border_a=5.0)
    cropped = pp.crop(pp.read_map(_write(m, tmp_path)), box)
    out = _reread(cropped, tmp_path, "cropped.map")
    # stored sub-block is smaller than the full grid...
    assert out.grid.nu < 72
    # ...but the header still says the cell is sampled on the full 72 grid,
    # so the voxel spacing (and thus registration) is unchanged.
    assert out.header_i32(8) == 72                 # MX
    assert np.isclose(out.header_float(11) / out.header_i32(8), 1.0, atol=1e-4)
    # and it is offset into the cell by a non-zero nxstart
    assert out.header_i32(5) > 0


def test_mask_aligns_with_cropped_map(tmp_path):
    m = _blob_map()
    pdb = _three_atom_pdb(tmp_path, [(30, 30, 34), (42, 42, 46), (36, 38, 40)])
    box = pp.model_frac_box(pdb, m.grid.unit_cell, border_a=5.0)
    mask = pp.atom_mask(pp.read_map(_write(m, tmp_path)), pdb, radius=3.0, frac_box=box)
    density = pp.crop(pp.read_map(_write(m, tmp_path)), box)
    mo = _reread(mask, tmp_path, "mask.map")
    do = _reread(density, tmp_path, "dens.map")
    # same stored dimensions and same nxstart => voxel-for-voxel aligned
    assert (mo.grid.nu, mo.grid.nv, mo.grid.nw) == (do.grid.nu, do.grid.nv, do.grid.nw)
    assert (mo.header_i32(5), mo.header_i32(6), mo.header_i32(7)) == \
           (do.header_i32(5), do.header_i32(6), do.header_i32(7))
    assert np.array(mo.grid, copy=False).sum() > 0


def test_model_frac_box_clamped_to_cell(tmp_path):
    m = _blob_map()
    # atoms hard against an edge: a 5 A border must not push the box below 0
    pdb = _three_atom_pdb(tmp_path, [(1, 1, 1), (3, 3, 3), (2, 2, 2)])
    box = pp.model_frac_box(pdb, m.grid.unit_cell, border_a=5.0)
    assert box.minimum.x >= 0.0 and box.minimum.y >= 0.0 and box.minimum.z >= 0.0


def test_map_model_cc_discriminates_the_hand(tmp_path):
    """CC is high when the model matches the map, low against the flipped map."""
    # A map computed from a chiral model; the correct hand must win decisively.
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(60, 60, 60, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for i in range(30):
        res = gemmi.Residue()
        res.name = "ALA"
        res.seqid = gemmi.SeqId(str(i + 1))
        at = gemmi.Atom()
        at.name = "CA"
        at.element = gemmi.Element("C")
        t = i * 0.6
        at.pos = gemmi.Position(30 + 6 * np.cos(t), 30 + 6 * np.sin(t), 15 + i * 0.9)
        at.b_iso = 20
        res.add_atom(at)
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    pdb = tmp_path / "model.pdb"
    st.write_pdb(str(pdb))

    dc = gemmi.DensityCalculatorX()
    dc.d_min = 2.5
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    m = gemmi.Ccp4Map()
    m.grid = dc.grid
    m.update_ccp4_header()
    m.setup(float("nan"), gemmi.MapSetup.Full)

    cc_correct = pp.map_model_cc(m, str(pdb), d_min=3.0)
    cc_wrong = pp.map_model_cc(pp.flip_hand(m), str(pdb), d_min=3.0)
    assert cc_correct > 0.5           # the model fits its own map
    assert cc_correct > 3 * cc_wrong  # and clearly beats the wrong hand


def test_empty_model_raises(tmp_path):
    m = _blob_map()
    empty = tmp_path / "empty.pdb"
    empty.write_text("END\n")
    with pytest.raises(ValueError):
        pp.model_frac_box(empty, m.grid.unit_cell, border_a=5.0)


def _write(m, tmp_path):
    path = tmp_path / "full.map"
    m.write_ccp4_map(str(path))
    return str(path)
