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
