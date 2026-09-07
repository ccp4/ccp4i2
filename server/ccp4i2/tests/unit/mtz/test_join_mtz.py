"""CPluginScript.joinMtz joins by HKL, as cmtzjoin did.

Observations and a Free-R set written by the same aimless run differ by a
reflection; the first gemmi port stacked columns as arrays and fell over
on exactly that (pointless_reindexToMatch, in SubstituteLigand's Phaser
route). Needs gemmi only.
"""
import math

import numpy as np
import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.CCP4PluginScript import CPluginScript


class _Joiner(CPluginScript):
    TASKNAME = "join_mtz_test"


def _write(path, columns, rows):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    mtz.add_dataset("d")
    for label, kind in columns:
        mtz.add_column(label, kind)
    mtz.set_data(np.array(rows, dtype=np.float32))
    mtz.update_reso()
    mtz.write_to_file(str(path))
    return path


@pytest.fixture
def obs_and_free(tmp_path):
    # Reflections allowed in P 21 21 21 (odd axial ones are absent)
    obs = _write(tmp_path / "obs.mtz", [("F", "F"), ("SIGF", "Q")], [
        [2, 0, 0, 100.0, 5.0], [0, 2, 0, 150.0, 7.0], [5, 6, 7, 200.0, 10.0], [8, 10, 11, 180.0, 9.0]])
    # The Free-R set lacks (8,10,11) and has (2,3,4) the observations lack
    free = _write(tmp_path / "free.mtz", [("FreeR_flag", "I")], [
        [2, 0, 0, 0], [0, 2, 0, 1], [5, 6, 7, 0], [2, 3, 4, 1]])
    return obs, free


def _rows(path):
    mtz = gemmi.read_mtz_file(str(path))
    labels = [c.label for c in mtz.columns]
    return labels, {tuple(int(x) for x in row[:3]): dict(zip(labels[3:], row[3:]))
                    for row in np.array(mtz, copy=False)}


def test_differing_reflection_lists_join_by_hkl(tmp_path, obs_and_free):
    obs, free = obs_and_free
    plugin = _Joiner(workDirectory=str(tmp_path))
    out = tmp_path / "joined.mtz"
    rv = plugin.joinMtz(str(out), [(str(obs), "F,SIGF", "F,SIGF"), (str(free), "FreeR_flag", "FreeR_flag")])
    assert rv == CPluginScript.SUCCEEDED
    labels, rows = _rows(out)
    assert labels[3:] == ["F", "SIGF", "FreeR_flag"]
    assert rows[(2, 0, 0)]["F"] == 100.0 and rows[(2, 0, 0)]["FreeR_flag"] == 0
    assert rows[(5, 6, 7)]["SIGF"] == 10.0 and rows[(5, 6, 7)]["FreeR_flag"] == 0
    # Present in one file only: the value it has, missing for the other
    assert rows[(8, 10, 11)]["F"] == 180.0 and math.isnan(rows[(8, 10, 11)]["FreeR_flag"])
    assert rows[(2, 3, 4)]["FreeR_flag"] == 1 and math.isnan(rows[(2, 3, 4)]["F"])


def test_columns_can_be_renamed(tmp_path, obs_and_free):
    obs, free = obs_and_free
    plugin = _Joiner(workDirectory=str(tmp_path))
    out = tmp_path / "renamed.mtz"
    assert plugin.joinMtz(str(out), [(str(obs), "F,SIGF", "FP,SIGFP"), (str(free), "FreeR_flag")]) \
        == CPluginScript.SUCCEEDED
    labels, rows = _rows(out)
    assert labels[3:] == ["FP", "SIGFP", "FreeR_flag"]
    assert rows[(0, 2, 0)]["FP"] == 150.0


def test_a_missing_file_fails_without_raising(tmp_path, obs_and_free):
    obs, _ = obs_and_free
    plugin = _Joiner(workDirectory=str(tmp_path))
    assert plugin.joinMtz(str(tmp_path / "x.mtz"), [(str(obs), "F,SIGF"), (str(tmp_path / "gone.mtz"), "FreeR_flag")]) \
        == CPluginScript.FAILED
