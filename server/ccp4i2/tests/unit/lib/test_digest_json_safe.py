"""A digest must be JSON, whatever the file header holds.

gemmi's CifToMtz writes a NaN dataset wavelength when the structure-factor
mmCIF records none. JSON has no NaN, so the digest of such an MTZ made the
renderer raise and the ``digest`` endpoint answer 500 -- seen on every
observation file of the demo campaign, whose data are imported from PDB
``-sf.cif`` files. CCP4-free: gemmi only.
"""

import json

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.core.CCP4XtalData import CObsDataFile
from ccp4i2.lib.utils.files.digest import digest_file_object, json_safe


def test_json_safe_replaces_non_finite_floats():
    raw = {"a": float("nan"), "b": [1.5, float("inf")], "c": (float("-inf"),), "d": "x"}
    assert json_safe(raw) == {"a": None, "b": [1.5, None], "c": [None], "d": "x"}


@pytest.fixture
def nan_wavelength_mtz(tmp_path):
    """An F/SIGF MTZ whose dataset wavelength is NaN, as CifToMtz writes it."""
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup("C 2 2 21")
    mtz.set_cell_for_all(gemmi.UnitCell(80.9, 96.2, 57.8, 90, 90, 90))
    mtz.add_dataset("unknown").wavelength = float("nan")
    mtz.add_column("F", "F")
    mtz.add_column("SIGF", "Q")
    mtz.set_data(np.array([[1, 1, 1, 100.0, 5.0], [2, 0, 2, 50.0, 3.0]], dtype="float32"))
    path = tmp_path / "nan_wavelength.mtz"
    mtz.write_to_file(str(path))
    return str(path)


def test_digest_of_mtz_with_nan_wavelength_is_json(nan_wavelength_mtz):
    f = CObsDataFile()
    f.setFullPath(nan_wavelength_mtz)
    digest = digest_file_object(f)
    assert digest.get("status") != "Failed", digest
    json.dumps(digest, allow_nan=False)
