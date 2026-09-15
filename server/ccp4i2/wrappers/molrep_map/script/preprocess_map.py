"""Origin-aware cryo-EM map surgery for the ``molrep_map`` placement task.

Pure ``gemmi`` + ``numpy`` -- no CCP4 binaries, no ``chapi``/``coot``. The module
is dependency-light on purpose so its unit tests run in the CCP4-free venv, even
though the plugin that calls it runs under ``ccp4-python`` at execution time.

The crystallographic rationale lives in ``docs/molrep-map-design.md``; the points
that shape this code:

* **Hand flip = inversion through the origin** (``V'[i] = V[(N-i) % N]``), which
  is what distinguishes the two cryo-EM hands. molrep re-places the model into
  each hand independently, so the flip only has to invert the density correctly;
  registration of the *output* comes from the placed model, not the flip origin.

* **The trimmed output keeps the full cell and encodes its position via
  ``nxstart``**, not the MRC-2000 Angstrom ``ORIGIN``. This is what coot/clipper
  -- hence Moorhen -- actually read (confirmed in the Moorhen source: its MRC
  parser stops before the Angstrom ORIGIN words, and density is placed by clipper
  in the grid/cell frame). ``gemmi.Ccp4Map.set_extent`` produces exactly this: a
  standard CCP4 partial map whose header retains the full-cell sampling
  (``MX/MY/MZ``) and voxel spacing, with a stored sub-block at a grid ``nxstart``.

* **Down-sampling is a disposable, placement-only speed lever** -- it coarsens the
  grid for molrep's search and is never applied to the deliverable map.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import gemmi

logger = logging.getLogger(__name__)

_P1 = "P 1"


def composite_size(n: int, primes=(2, 3, 5, 7)) -> int:
    """Smallest integer >= ``n`` whose prime factors are all in ``primes``.

    FFT-friendly grid dimensions (products of small primes) keep any downstream
    transform -- molrep's, servalcat's, or a test's -- efficient and exact.
    """
    n = max(1, int(n))
    while True:
        m = n
        for p in primes:
            while m % p == 0:
                m //= p
        if m == 1:
            return n
        n += 1


def _spacegroup(m: gemmi.Ccp4Map) -> gemmi.SpaceGroup:
    sg = m.grid.spacegroup
    return sg if sg is not None else gemmi.SpaceGroup(_P1)


def _ccp4_from_array(arr: np.ndarray, cell: gemmi.UnitCell,
                     spacegroup: gemmi.SpaceGroup) -> gemmi.Ccp4Map:
    """Wrap a ``(nu, nv, nw)`` numpy array as a fresh full-cell CCP4 map."""
    arr = np.ascontiguousarray(arr, dtype="float32")
    grid = gemmi.FloatGrid(*arr.shape)
    grid.set_unit_cell(cell)
    grid.spacegroup = spacegroup
    np.array(grid, copy=False)[:] = arr
    m = gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header()
    return m


def _frac_to_cart_matrix(cell: gemmi.UnitCell) -> np.ndarray:
    """3x3 matrix M such that ``cartesian = M @ fractional`` (any cell)."""
    o0 = cell.orthogonalize(gemmi.Fractional(0, 0, 0))
    cols = []
    for f in (gemmi.Fractional(1, 0, 0),
              gemmi.Fractional(0, 1, 0),
              gemmi.Fractional(0, 0, 1)):
        p = cell.orthogonalize(f)
        cols.append([p.x - o0.x, p.y - o0.y, p.z - o0.z])
    return np.array(cols, dtype="float64").T


def read_map(path, normalize: bool = True) -> gemmi.Ccp4Map:
    """Read a CCP4/MRC map. With ``normalize`` (the default) the axes are put in
    standard order and the map is expanded to the full cell -- the right footing
    for building the disposable molrep-input map. Read the *output* frame with
    ``normalize=False`` (or ``MapSetup.ReorderOnly``) so a partial map is not
    re-expanded.
    """
    m = gemmi.read_ccp4_map(str(path))
    if normalize:
        m.setup(float("nan"), gemmi.MapSetup.Full)
    else:
        m.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
    return m


def flip_hand(m: gemmi.Ccp4Map) -> gemmi.Ccp4Map:
    """Invert the map through the origin -- the cryo-EM hand flip.

    ``V'[i, j, k] = V[(N-i) % N, (N-j) % N, (N-k) % N]``, i.e. reverse each axis
    and roll by one so grid point 0 maps to itself. Parity against
    ``coot_headless_api.flip_hand`` is checked in ``tests/parity/``.
    """
    arr = np.array(m.grid, copy=True)
    flipped = np.roll(arr[::-1, ::-1, ::-1], 1, axis=(0, 1, 2))
    return _ccp4_from_array(flipped, m.grid.unit_cell, _spacegroup(m))


def downsample(m: gemmi.Ccp4Map, factor: float) -> gemmi.Ccp4Map:
    """Resample onto a coarser grid spanning the *same* cell (disposable, lossy).

    A factor of ``2`` roughly halves each grid dimension and cuts an FFT-based
    search by ~8x. The cell is unchanged, so the result is a valid full periodic
    cell at coarser sampling -- exactly what a fast MR search wants. Never call
    this on a deliverable map.
    """
    if factor is None or factor <= 1:
        return m
    g = m.grid
    cell = g.unit_cell
    new_dims = [max(1, composite_size(round(n / factor)))
                for n in (g.nu, g.nv, g.nw)]
    fu = np.arange(new_dims[0]) / new_dims[0]
    fv = np.arange(new_dims[1]) / new_dims[1]
    fw = np.arange(new_dims[2]) / new_dims[2]
    FU, FV, FW = np.meshgrid(fu, fv, fw, indexing="ij")
    frac = np.stack([FU.ravel(), FV.ravel(), FW.ravel()], axis=-1)  # (Npts, 3)
    M = _frac_to_cart_matrix(cell)
    cart = frac @ M.T
    vals = g.interpolate_position_array(np.ascontiguousarray(cart), order=1)
    arr = np.asarray(vals, dtype="float32").reshape(new_dims)
    logger.debug("downsample %s -> %s (factor %s)",
                 (g.nu, g.nv, g.nw), tuple(new_dims), factor)
    return _ccp4_from_array(arr, cell, _spacegroup(m))


def blur(m: gemmi.Ccp4Map, b_add: float) -> gemmi.Ccp4Map:
    """Apply an isotropic B-factor via an FFT round-trip (orthogonal cells).

    ``F(s) *= exp(-B s^2 / 4)`` with ``s^2 = 1/d^2``. A positive ``b_add`` blurs,
    emphasising the low-resolution envelope MR keys on. Placement-only.
    """
    if not b_add:
        return m
    g = m.grid
    cell = g.unit_cell
    a, b, c = cell.parameters[:3]
    if not (abs(cell.alpha - 90) < 1e-3 and abs(cell.beta - 90) < 1e-3
            and abs(cell.gamma - 90) < 1e-3):
        logger.warning("blur() assumes an orthogonal cell; angles are "
                       "%.2f/%.2f/%.2f -- skipping blur",
                       cell.alpha, cell.beta, cell.gamma)
        return m
    arr = np.array(g, copy=False)
    nu, nv, nw = arr.shape
    F = np.fft.rfftn(arr)
    ku = np.fft.fftfreq(nu, d=a / nu)
    kv = np.fft.fftfreq(nv, d=b / nv)
    kw = np.fft.rfftfreq(nw, d=c / nw)
    s2 = (ku[:, None, None] ** 2 + kv[None, :, None] ** 2
          + kw[None, None, :] ** 2)
    F *= np.exp(-float(b_add) * s2 / 4.0)
    out = np.fft.irfftn(F, s=(nu, nv, nw), axes=(0, 1, 2)).astype("float32")
    return _ccp4_from_array(out, cell, _spacegroup(m))


def prepare_for_search(m: gemmi.Ccp4Map, downsample_factor: float = 1.0,
                       b_add: float = 0.0) -> gemmi.Ccp4Map:
    """The disposable molrep-input map: optional down-sample then blur."""
    out = downsample(m, downsample_factor)
    out = blur(out, b_add)
    return out


def model_frac_box(structure_path, cell: gemmi.UnitCell,
                   border_a: float) -> gemmi.FractionalBox:
    """Fractional bounding box of the model plus a solvent border (Angstrom).

    The box is clamped to the unit cell so a molecule near an edge cannot ask for
    an extent outside it.
    """
    st = gemmi.read_structure(str(structure_path))
    box = gemmi.FractionalBox()
    n_atoms = 0
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    box.extend(cell.fractionalize(atom.pos))
                    n_atoms += 1
    if n_atoms == 0:
        raise ValueError(f"No atoms in {structure_path} to bound the map to")
    a, b, c = cell.parameters[:3]
    lo, hi = box.minimum, box.maximum
    padded = gemmi.FractionalBox()
    padded.extend(gemmi.Fractional(
        max(0.0, lo.x - border_a / a),
        max(0.0, lo.y - border_a / b),
        max(0.0, lo.z - border_a / c)))
    padded.extend(gemmi.Fractional(
        min(1.0, hi.x + border_a / a),
        min(1.0, hi.y + border_a / b),
        min(1.0, hi.z + border_a / c)))
    return padded


def model_ortho_box(structure_path, border_a: float):
    """Orthogonal-Angstrom bounding box ``(min_xyz, max_xyz)`` of the model plus a
    solvent border. Frame-independent -- the caller converts it into whichever
    map's fractional coordinates it needs (see :func:`frac_box_from_ortho`)."""
    st = gemmi.read_structure(str(structure_path))
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    n_atoms = 0
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    p = atom.pos
                    for i, v in enumerate((p.x, p.y, p.z)):
                        if v < lo[i]:
                            lo[i] = v
                        if v > hi[i]:
                            hi[i] = v
                    n_atoms += 1
    if n_atoms == 0:
        raise ValueError(f"No atoms in {structure_path} to bound the map to")
    lo = tuple(v - border_a for v in lo)
    hi = tuple(v + border_a for v in hi)
    return lo, hi


def frac_box_from_ortho(cell: gemmi.UnitCell, ortho_min, ortho_max,
                        offset=(0.0, 0.0, 0.0)) -> gemmi.FractionalBox:
    """Fractional box in ``cell``'s frame from an orthogonal-Angstrom box, after
    translating it by ``offset`` (Angstrom).

    ``offset`` reconciles frames that a re-box reset to a common origin: a map
    cropped to a smaller centred box sits at ``(big_cell - small_cell)/2`` in the
    larger map's frame, and that lost registration is not in either header. All
    eight corners are fractionalised (correct for any cell) and clamped to the
    unit interval.
    """
    box = gemmi.FractionalBox()
    for x in (ortho_min[0], ortho_max[0]):
        for y in (ortho_min[1], ortho_max[1]):
            for z in (ortho_min[2], ortho_max[2]):
                f = cell.fractionalize(
                    gemmi.Position(x + offset[0], y + offset[1], z + offset[2]))
                box.extend(gemmi.Fractional(
                    min(1.0, max(0.0, f.x)),
                    min(1.0, max(0.0, f.y)),
                    min(1.0, max(0.0, f.z))))
    return box


def centred_crop_offset(inner_cell: gemmi.UnitCell,
                        outer_cell: gemmi.UnitCell):
    """Angstrom shift taking a point in a centred-crop map's frame to the same
    physical point in the larger (outer) map's frame: ``(outer - inner)/2`` per
    axis. Zero when the cells match. Assumes the standard central re-box (both
    maps origin-reset), which is what EMDB half-map/primary pairs use."""
    ip = inner_cell.parameters
    op = outer_cell.parameters
    return tuple((op[i] - ip[i]) / 2.0 for i in range(3))


def crop(m: gemmi.Ccp4Map, frac_box: gemmi.FractionalBox) -> gemmi.Ccp4Map:
    """Crop a full-cell map to ``frac_box`` in place, returning it.

    ``set_extent`` keeps the full-cell sampling (``MX/MY/MZ``) and voxel spacing
    and stores only the sub-block, positioned by ``nxstart`` -- the Moorhen/clipper
    -readable representation. Do not call ``update_ccp4_header`` afterwards (it
    requires a re-``setup``); ``write_ccp4_map`` writes the correct header.
    """
    m.set_extent(frac_box)
    return m


def atom_mask(reference_full_map: gemmi.Ccp4Map, structure_path,
              radius: float, frac_box: Optional[gemmi.FractionalBox] = None,
              value: float = 1.0) -> gemmi.Ccp4Map:
    """A hard atom mask on the reference map's full grid, optionally cropped to
    ``frac_box`` so it aligns voxel-for-voxel with a map cropped to the same box.

    Built on the *full* grid (where sampling is correct) then cropped, because a
    fresh grid sized to the sub-block would misplace ``set_points_around``.
    servalcat can soften the edge itself (``--mask_soft_edge``); a hard mask in.
    """
    g = reference_full_map.grid
    mgrid = gemmi.FloatGrid(g.nu, g.nv, g.nw)
    mgrid.set_unit_cell(g.unit_cell)
    mgrid.spacegroup = _spacegroup(reference_full_map)
    st = gemmi.read_structure(str(structure_path))
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    mgrid.set_points_around(atom.pos, radius=radius, value=value)
    mask = gemmi.Ccp4Map()
    mask.grid = mgrid
    mask.update_ccp4_header()
    if frac_box is not None:
        mask.set_extent(frac_box)
    return mask


def map_model_cc(m: gemmi.Ccp4Map, structure_path, d_min: float,
                 mask_radius: float = 2.5) -> float:
    """Real-space correlation between a (full-cell) map and the density calculated
    from the placed model, over a mask around the model.

    This is the honest hand-discriminator: it uses the map's *phases* (real space),
    so unlike an amplitude score it actually distinguishes the two hands -- the
    correct hand's model overlays the density and correlates; the wrong hand does
    not. Returns Pearson r in ``[-1, 1]`` (NaN if it cannot be computed). The map
    must span its full cell (read with ``MapSetup.Full``); pass the placement map
    *before* it is cropped.
    """
    grid = m.grid
    cell = grid.unit_cell
    st = gemmi.read_structure(str(structure_path))
    st.cell = cell
    st.spacegroup_hm = "P 1"
    dc = gemmi.DensityCalculatorX()
    dc.d_min = d_min
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    model_grid = dc.grid

    mask = gemmi.FloatGrid(grid.nu, grid.nv, grid.nw)
    mask.set_unit_cell(cell)
    mask.spacegroup = _spacegroup(m)
    for model in st:
        for chain in model:
            for res in chain:
                for atom in res:
                    mask.set_points_around(atom.pos, radius=mask_radius, value=1.0)
    mk = np.array(mask, copy=False) > 0.5
    if not mk.any():
        return float("nan")

    map_vals = np.array(grid, copy=False)[mk]
    idx = np.argwhere(mk)
    frac = idx / np.array([grid.nu, grid.nv, grid.nw])
    cart = np.ascontiguousarray((frac @ _frac_to_cart_matrix(cell).T).astype("float64"))
    model_vals = np.asarray(model_grid.interpolate_position_array(cart, order=1))

    good = np.isfinite(map_vals) & np.isfinite(model_vals)
    if good.sum() < 10:
        return float("nan")
    mv, dv = map_vals[good], model_vals[good]
    if mv.std() == 0 or dv.std() == 0:
        return float("nan")
    return float(np.corrcoef(mv, dv)[0, 1])


def write_map(m: gemmi.Ccp4Map, path) -> None:
    """Write a CCP4 map to ``path`` (str/Path)."""
    m.write_ccp4_map(str(path))
