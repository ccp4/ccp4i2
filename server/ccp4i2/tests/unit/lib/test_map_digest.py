"""Test the CMapDataFile digest.

A real-space CCP4/MRC map used to fall through the digest dispatch to the generic
handler and return ``null`` -- the file preview showed no header at all. The
``digest_cmapdatafile_file_object`` handler reports the map header (grid, cell,
spacing, start, spacegroup, stats, subType) instead. CCP4-free: gemmi only.
"""

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.core.CCP4XtalData import CMapDataFile
from ccp4i2.lib.utils.files.digest import (
    digest_cmapdatafile_file_object,
    digest_file_object,
)


@pytest.fixture
def map_file(tmp_path):
    """A synthetic non-cubic P1 map on disk, with a chosen subType."""
    cell = gemmi.UnitCell(72, 72, 90, 90, 90, 90)
    grid = gemmi.FloatGrid(48, 48, 60)
    grid.set_unit_cell(cell)
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    arr = np.array(grid, copy=False)
    arr[:] = np.random.default_rng(0).random(arr.shape).astype("float32")
    m = gemmi.Ccp4Map()
    m.grid = grid
    m.update_ccp4_header()
    path = tmp_path / "map.mrc"
    m.write_ccp4_map(str(path))
    return str(path)


def _digest(path, sub_type=None):
    f = CMapDataFile()
    f.setFullPath(path)
    if sub_type is not None:
        f.subType.set(sub_type)
    return f


def test_digest_reports_grid_cell_and_spacing(map_file):
    d = digest_cmapdatafile_file_object(_digest(map_file))
    assert d["grid_sampling"] == [48, 48, 60]
    assert d["cell"]["a"] == 72.0 and d["cell"]["c"] == 90.0
    assert d["spacing"] == [1.5, 1.5, 1.5]      # 72/48, 90/60
    assert d["spacegroup"] == 1
    assert d["likely_em"] is True               # P1 + orthogonal
    assert set(("min", "max", "mean", "rms")).issubset(d["statistics"])


def test_digest_reports_subtype_label(map_file):
    d = digest_cmapdatafile_file_object(_digest(map_file, CMapDataFile.SUBTYPE_MASK))
    assert d["sub_type"] == 4
    assert d["sub_type_label"] == "Mask"


def test_dispatch_routes_map_to_the_map_handler(map_file):
    # The generic dispatcher must reach the map handler (not the null fallback).
    d = digest_file_object(_digest(map_file, CMapDataFile.SUBTYPE_NORMAL))
    assert d.get("format") == "CCP4/MRC map"
    assert d["sub_type_label"] == "Normal (electron density)"


def test_digest_of_unset_file_fails_cleanly():
    d = digest_cmapdatafile_file_object(CMapDataFile())
    assert d["status"] == "Failed"
