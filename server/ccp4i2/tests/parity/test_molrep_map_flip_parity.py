"""Parity: the gemmi hand-flip in ``molrep_map`` must equal coot's ``flip_hand``.

The Django port replaces the Qt task's ``chapi.flip_hand`` with a gemmi point
inversion (no coot dependency on the request path). coot/Moorhen is the eventual
consumer of the flipped density, so the gemmi flip has to reproduce coot's flip
exactly -- not merely "an inversion". This test reads the same synthetic map
through both and asserts they are voxel-for-voxel identical.

Needs ``coot_headless_api`` (CCP4); skipped in the CCP4-free venv.
"""

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")
chapi = pytest.importorskip("coot_headless_api",
                            reason="needs coot_headless_api (CCP4)")

from ccp4i2.wrappers.molrep_map.script import preprocess_map as pp


def _asymmetric_map(tmp_path, n=48):
    cell = gemmi.UnitCell(n, n, n, 90, 90, 90)
    grid = gemmi.FloatGrid(n, n, n)
    grid.set_unit_cell(cell)
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    arr = np.array(grid, copy=False)
    zz, yy, xx = np.mgrid[0:n, 0:n, 0:n]
    arr[:] = (np.exp(-(((xx - 14) ** 2 + (yy - 18) ** 2 + (zz - 22) ** 2) / 30.))
              + 0.5 * np.exp(-(((xx - 30) ** 2 + (yy - 10) ** 2
                                + (zz - 26) ** 2) / 20.))).astype("float32")
    m = gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header()
    m.setup(float("nan"), gemmi.MapSetup.Full)
    path = tmp_path / "input.map"
    m.write_ccp4_map(str(path))
    return m, str(path)


def test_gemmi_flip_matches_coot_flip_hand(tmp_path):
    m, path = _asymmetric_map(tmp_path)

    # coot's reference flip
    mc = chapi.molecules_container_t(False)
    imol = mc.read_ccp4_map(path, False)
    iflip = mc.flip_hand(imol)
    coot_out = str(tmp_path / "coot_flipped.map")
    mc.write_map(iflip, coot_out)
    coot = gemmi.read_ccp4_map(coot_out)
    coot.setup(float("nan"), gemmi.MapSetup.ReorderOnly)

    # the port's flip
    mine = pp.flip_hand(m)

    a_coot = np.array(coot.grid, copy=False)
    a_mine = np.array(mine.grid, copy=False)
    assert a_coot.shape == a_mine.shape
    assert np.allclose(a_coot, a_mine, atol=1e-4)
