"""Extending a FreeR set across the crystals of a campaign.

In COMPLETE mode the freerflag wrapper joins the observed data with the input
FreeR set by reflection index (``makeHklin`` -> ``merge_mtz_files``, gemmi).
That join refused any pair of cells outside Clipper's 1 A test, which a
campaign's shared free set fails on most of its crystals: the CDK4/CyclinD1
pre-screen of September 2026 lost 23 of 26 datasets to it, after the
aimless_pipe pre-check had already been overridden.

These tests cover the permissive path: ``merge_mtz_files`` keeps its default
check, skips it when told, and keeps the first file's cell; the freerflag
wrapper exposes the override; aimless_pipe hands its own override down.
"""
import gemmi
import numpy as np
import pytest

from ccp4i2.core.CCP4Utils import MtzMergeError, merge_mtz_files
from ccp4i2.core.tasks import get_plugin_class


def _write_mtz(path, cell, columns):
    """A small P 1 21 1 MTZ with the given cell and (label, type, values)."""
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 1 21 1")
    mtz.set_cell_for_all(gemmi.UnitCell(*cell))
    mtz.add_dataset("data")
    for label, ctype, _ in columns:
        mtz.add_column(label, ctype)
    hkl = [(h, k, l) for h in range(0, 3) for k in range(0, 3) for l in range(1, 4)]
    rows = []
    for i, (h, k, l) in enumerate(hkl):
        rows.append([h, k, l] + [vals[i % len(vals)] for _, _, vals in columns])
    mtz.set_data(np.array(rows, dtype=np.float32))
    mtz.write_to_file(str(path))
    return path


def _data_and_freer(tmp_path, freer_cell):
    data = _write_mtz(tmp_path / "data.mtz", (58.0, 64.4, 186.2, 90.0, 90.7, 90.0),
                      [("F", "F", [10.0, 20.0]), ("SIGF", "Q", [1.0, 2.0])])
    freer = _write_mtz(tmp_path / "freer.mtz", freer_cell,
                       [("FreeR_flag", "I", [0, 1, 2, 3, 4])])
    return data, freer


def _specs(data, freer):
    return [
        {"path": data, "column_mapping": {"F": "F", "SIGF": "SIGF"}},
        {"path": freer, "column_mapping": {"FreeR_flag": "FreeR_flag"}},
    ]


CAMPAIGN_FREER_CELL = (57.6, 64.4, 184.7, 90.0, 91.7, 90.0)  # 1.5 A off on c


def test_the_default_check_refuses_a_campaign_freer_cell(tmp_path):
    data, freer = _data_and_freer(tmp_path, CAMPAIGN_FREER_CELL)
    with pytest.raises(MtzMergeError, match="Incompatible unit cells"):
        merge_mtz_files(_specs(data, freer), tmp_path / "out.mtz")


def test_no_tolerance_joins_by_index_and_keeps_the_first_cell(tmp_path):
    data, freer = _data_and_freer(tmp_path, CAMPAIGN_FREER_CELL)
    out = merge_mtz_files(_specs(data, freer), tmp_path / "out.mtz", cell_tolerance=None)
    mtz = gemmi.read_mtz_file(str(out))
    assert [c.label for c in mtz.columns][-3:] == ["F", "SIGF", "FreeR_flag"]
    assert mtz.cell.c == pytest.approx(186.2, abs=0.01)
    assert mtz.cell.beta == pytest.approx(90.7, abs=0.01)


def test_a_wider_tolerance_also_admits_it(tmp_path):
    data, freer = _data_and_freer(tmp_path, CAMPAIGN_FREER_CELL)
    merge_mtz_files(_specs(data, freer), tmp_path / "out.mtz", cell_tolerance=5.0)


def test_space_groups_must_still_agree_without_a_cell_check(tmp_path):
    data, freer = _data_and_freer(tmp_path, CAMPAIGN_FREER_CELL)
    mtz = gemmi.read_mtz_file(str(freer))
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.write_to_file(str(freer))
    with pytest.raises(MtzMergeError, match="Incompatible space groups"):
        merge_mtz_files(_specs(data, freer), tmp_path / "out.mtz", cell_tolerance=None)


def _plugin(tmp_path, name):
    plugin_class = get_plugin_class(name)
    assert plugin_class is not None, f"{name} is not in the registry"
    directory = tmp_path / name
    directory.mkdir(exist_ok=True)
    return plugin_class(workDirectory=str(directory), name=name)


def test_the_freerflag_wrapper_exposes_the_override_off_by_default(tmp_path):
    plugin = _plugin(tmp_path, "freerflag")
    assert not plugin.container.controlParameters.OVERRIDE_CELL_DIFFERENCE


def test_aimless_pipe_hands_its_override_to_freerflag(tmp_path):
    pipe = _plugin(tmp_path, "aimless_pipe")
    pipe.container.outputData.HKLOUT.addItem()
    pipe.container.controlParameters.OVERRIDE_CELL_DIFFERENCE.set(True)
    freerflag = pipe.makePluginObject("freerflag")
    pipe._configureFreerflag(freerflag, complete=True)
    assert freerflag.container.controlParameters.OVERRIDE_CELL_DIFFERENCE
    assert str(freerflag.container.controlParameters.GEN_MODE) == "COMPLETE"


def test_aimless_pipe_leaves_freerflag_strict_unless_asked(tmp_path):
    pipe = _plugin(tmp_path, "aimless_pipe")
    pipe.container.outputData.HKLOUT.addItem()
    freerflag = pipe.makePluginObject("freerflag")
    pipe._configureFreerflag(freerflag, complete=False)
    assert not freerflag.container.controlParameters.OVERRIDE_CELL_DIFFERENCE
    assert str(freerflag.container.controlParameters.GEN_MODE) == "GEN_NEW"
