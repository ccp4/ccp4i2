"""A structure given to Phaser as already placed must be in THIS crystal.

Haiku, solving CDK4/cyclin D1, gave the cyclin chain of a homologous entry
(6p8e, cell 62.4 67.5 187.3) as the fixed structure for data of cell
57.8 64.7 186.1; Phaser held it in that frame and placed CDK4 against it
(LLG -892). The structure's cell and point group are now checked against
the data's before the job runs.
"""
from pathlib import Path

import gemmi
import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _model(path, cell, space_group):
    atoms = [l for l in (DEMO / "beta.pdb").read_text().splitlines() if l.startswith("ATOM")]
    path.write_text(f"CRYST1{cell[0]:9.3f}{cell[1]:9.3f}{cell[2]:9.3f}"
                    f"{cell[3]:7.2f}{cell[4]:7.2f}{cell[5]:7.2f} {space_group:<11s}\n"
                    + "\n".join(atoms) + "\nEND\n")
    return path


def _plugin(tmp_path, fixed, input_fixed=True):
    plugin = get_plugin_class("phaser_simple_phil")(workDirectory=str(tmp_path), name="mr")
    inp = plugin.container.inputData
    inp.F_SIGF.setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    inp.XYZIN.setFullPath(str(DEMO / "blip.pdb"))
    inp.INPUT_FIXED.set(input_fixed)
    inp.XYZIN_FIXED.setFullPath(str(fixed))
    return plugin


def _blocking_on_fixed(error):
    return [r for r in error._reports
            if str(r.get("name", "")).endswith("XYZIN_FIXED")
            and r["severity"] >= CCP4ErrorHandling.SEVERITY_ERROR]


@pytest.fixture
def data_cell():
    mtz = gemmi.read_mtz_file(str(DEMO / "beta_blip_P3221.mtz"))
    c = mtz.cell
    return (c.a, c.b, c.c, c.alpha, c.beta, c.gamma)


def test_a_structure_from_another_crystal_is_refused(tmp_path, data_cell):
    other = (data_cell[0] * 1.08, data_cell[1] * 1.08, data_cell[2]) + data_cell[3:]
    fixed = _model(tmp_path / "homologue.pdb", other, "P 32 2 1")
    assert _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_predicted_model_is_refused(tmp_path):
    # AlphaFold models carry a 1 A cell in P 1: placed in no crystal
    fixed = _model(tmp_path / "AF.pdb", (1, 1, 1, 90, 90, 90), "P 1")
    assert _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_structure_placed_in_this_crystal_is_accepted(tmp_path, data_cell):
    fixed = _model(tmp_path / "placed.pdb", data_cell, "P 32 2 1")
    assert not _blocking_on_fixed(_plugin(tmp_path, fixed).runTimeValidity())


def test_a_file_left_set_without_input_fixed_is_ignored(tmp_path):
    fixed = _model(tmp_path / "AF.pdb", (1, 1, 1, 90, 90, 90), "P 1")
    assert not _blocking_on_fixed(_plugin(tmp_path, fixed, input_fixed=False).runTimeValidity())
