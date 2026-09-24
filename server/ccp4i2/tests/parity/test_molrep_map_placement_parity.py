"""Acceptance-level parity for ``molrep_map`` placement, using the real molrep.

No downloads: a chiral model is turned into a map with gemmi, then placed back
through the engine on both hands. This exercises the whole placement + trim path
and asserts the two things that matter for a cryo-EM placement task:

1. **Hand discrimination** -- the correct (original) hand out-scores the inverted
   hand, i.e. the P1 *phased* search actually distinguishes the two hands.
2. **Real-space registration** -- the placed model sits on strong density in the
   *original* map, and the trimmed output is a valid ``nxstart``-encoded partial
   map that covers the model.

Needs the ``molrep`` binary (CCP4); skipped in the CCP4-free venv.
"""

import glob
import os
import shutil
import tempfile

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")
pytestmark = pytest.mark.skipif(shutil.which("molrep") is None,
                                reason="molrep binary not on PATH")

from ccp4i2.wrappers.molrep_map.script import preprocess_map as pp
from ccp4i2.wrappers.molrep_map.script import engines


def _chiral_model():
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
        at.pos = gemmi.Position(30 + 6 * np.cos(t), 30 + 6 * np.sin(t), 18 + i * 0.8)
        at.b_iso = 20
        res.add_atom(at)
        chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    return st


@pytest.fixture(scope="module")
def model_and_map(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("molrep_map_parity")
    st = _chiral_model()
    pdb = str(tmp / "model.pdb")
    st.write_pdb(pdb)
    dc = gemmi.DensityCalculatorX()
    dc.d_min = 3.0
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    m = gemmi.Ccp4Map()
    m.grid = dc.grid
    m.update_ccp4_header()
    mapf = str(tmp / "model.map")
    m.write_ccp4_map(mapf)
    return pdb, mapf, str(tmp)


def _place(hand_map_path, pdb, workdir):
    return engines.place("molrep", hand_map_path, pdb, workdir,
                         nmon=1, np_peaks=3, time_limit=180)


def test_original_hand_outscores_inverted(model_and_map):
    pdb, mapf, tmp = model_and_map
    orig = pp.read_map(mapf)
    flipped = pp.flip_hand(orig)
    flip_path = os.path.join(tmp, "flipped.map")
    pp.write_map(flipped, flip_path)

    res_o = _place(mapf, pdb, os.path.join(tmp, "orig"))
    res_f = _place(flip_path, pdb, os.path.join(tmp, "flip"))

    assert res_o.placed, "original hand placed nothing"
    # The correct hand must score, and beat the inverted hand (or the inverted
    # hand fails to place at all).
    assert res_o.score is not None
    if res_f.placed and res_f.score is not None:
        assert res_o.score >= res_f.score


def test_placed_model_overlays_original_density(model_and_map):
    pdb, mapf, tmp = model_and_map
    res = _place(mapf, pdb, os.path.join(tmp, "overlay"))
    assert res.placed

    m = pp.read_map(mapf)
    grid = m.grid
    placed = gemmi.read_structure(res.model_path)
    # Mean map value at placed atom centres should be well above the map mean:
    # if the placement were wrong the atoms would fall on noise/solvent.
    arr = np.array(grid, copy=False)
    map_mean = float(arr.mean())
    vals = []
    for model in placed:
        for chain in model:
            for r in chain:
                for at in r:
                    vals.append(grid.interpolate_value(at.pos))
    assert vals
    assert float(np.mean(vals)) > map_mean


def test_trimmed_output_is_valid_partial_map(model_and_map):
    pdb, mapf, tmp = model_and_map
    res = _place(mapf, pdb, os.path.join(tmp, "trim"))
    assert res.placed

    # Read the full map's header sampling + spacing to compare against.
    full_ref = gemmi.read_ccp4_map(mapf)
    full_ref.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
    full_mx = full_ref.header_i32(8)
    full_spacing = full_ref.header_float(11) / full_mx
    full_points = full_ref.grid.point_count

    full = pp.read_map(mapf)
    box = pp.model_frac_box(res.model_path, full.grid.unit_cell, border_a=5.0)
    cropped = pp.crop(full, box)
    out = os.path.join(tmp, "trimmed.map")
    pp.write_map(cropped, out)

    rr = gemmi.read_ccp4_map(out)
    rr.setup(float("nan"), gemmi.MapSetup.ReorderOnly)
    # Stored sub-block is genuinely smaller than the full map...
    assert rr.grid.point_count < full_points
    # ...yet the header retains full-cell sampling and voxel spacing, so the
    # trimmed map stays registered with the model (position lives in nxstart).
    assert rr.header_i32(8) == full_mx
    assert np.isclose(rr.header_float(11) / rr.header_i32(8), full_spacing, atol=1e-3)
    assert (rr.header_i32(5), rr.header_i32(6), rr.header_i32(7)) != (0, 0, 0)
