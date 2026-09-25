"""``prepare_mtz_for_pandda`` finds the free-R set by MTZ column *type* when
no label matches. CCP4i2's own complete MTZs name the column after the
parameter that supplied it (``FREERFLAG_FREER``, dimple's COMPLETE_MTZ),
which no label list anticipates, and PanDDA accepts three exact labels."""
import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi")

from ccp4i2.lib.pandda_export import prepare_mtz_for_pandda


def _mtz(path, columns):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(40.0, 50.0, 60.0, 90.0, 90.0, 90.0))
    mtz.add_dataset("test")
    for label, ctype in columns:
        mtz.add_column(label, ctype)
    data = np.array([[1, 0, 0, 100.0, 5.0, 3], [2, 0, 0, 90.0, 4.0, 0]], dtype=np.float32)
    mtz.set_data(data[:, :3 + len(columns)])
    mtz.write_to_file(str(path))
    return path


def test_a_recognised_label_is_left_alone(tmp_path):
    src = _mtz(tmp_path / "a.mtz", [("F", "F"), ("SIGF", "Q"), ("FreeR_flag", "I")])
    assert prepare_mtz_for_pandda(src, tmp_path) == src


def test_ccp4i2s_own_naming_is_found_by_column_type(tmp_path):
    src = _mtz(tmp_path / "b.mtz", [("F_SIGF_F", "F"), ("F_SIGF_SIGF", "Q"), ("FREERFLAG_FREER", "I")])
    out = prepare_mtz_for_pandda(src, tmp_path)
    assert out != src
    labels = [c.label for c in gemmi.read_mtz_file(str(out)).columns]
    assert "FreeR_flag" in labels and "FREERFLAG_FREER" not in labels


def test_an_mtz_with_no_flag_column_is_returned_unchanged(tmp_path):
    src = _mtz(tmp_path / "c.mtz", [("F", "F"), ("SIGF", "Q"), ("ANOM", "D")])
    assert prepare_mtz_for_pandda(src, tmp_path) == src
