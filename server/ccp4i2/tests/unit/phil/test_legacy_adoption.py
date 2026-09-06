"""A PHIL task adopts the front page of a classic job: typed inputs by
name (or a declared rename), and the values that became PHIL parameters.
Everything else takes the PHIL defaults."""
from pathlib import Path

import pytest

pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

import ccp4i2
from ccp4i2.core.tasks import get_plugin_class

GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma"


def make(task, tmp_path):
    return get_plugin_class(task)(workDirectory=str(tmp_path / task), name=task)


def test_phaser_simple_front_page(tmp_path):
    old = make("phaser_simple", tmp_path)
    inp = old.container.inputData
    inp.F_SIGF.setFullPath(str(GAMMA / "merged_intensities_Xe.mtz"))
    inp.XYZIN.setFullPath(str(GAMMA / "gamma_model.pdb"))
    inp.NCOPIES.set(2)
    inp.ID_RMS.set("RMS")
    inp.SEARCHRMS.set(1.2)
    inp.RESOLUTION_HIGH.set(2.5)
    inp.RUNREFMAC.set(False)
    old.container.keywords.TOPF.set(99)   # not front page: left behind

    new = make("phaser_simple_phil", tmp_path)
    adopted = new.adopt_legacy_container(old.container)
    ninp = new.container.inputData
    assert str(ninp.XYZIN.fullPath) == str(GAMMA / "gamma_model.pdb")
    assert str(ninp.F_SIGF.fullPath) == str(GAMMA / "merged_intensities_Xe.mtz")
    assert int(ninp.NCOPIES) == 2 and str(ninp.ID_RMS) == "RMS" and float(ninp.SEARCHRMS) == 1.2
    assert not bool(ninp.RUNREFMAC)
    assert float(new.find_phil("phaser.keywords.resolution.high")) == 2.5
    assert new.find_phil("phaser.keywords.resolution.low").isSet(allowDefault=False) is False
    assert "XYZIN" in adopted and "RESOLUTION_HIGH -> phaser.keywords.resolution.high" in adopted
    assert not any("TOPF" in a for a in adopted)


def test_ensembles_come_across_as_a_list(tmp_path):
    old = make("phaser_MR_AUTO", tmp_path)
    ensembles = old.container.inputData.ENSEMBLES
    ensembles.append(ensembles.makeItem())
    ensembles[-1].label.set("mine")
    ensembles[-1].number.set(3)
    item = ensembles[-1].pdbItemList.makeItem()
    ensembles[-1].pdbItemList.append(item)
    item.structure.setFullPath(str(GAMMA / "gamma_model.pdb"))
    item.identity_to_target.set(0.8)

    new = make("phaser_mr_auto_phil", tmp_path)
    new.adopt_legacy_container(old.container)
    got = new.container.inputData.ENSEMBLES
    assert len(got) == 1 and str(got[0].label) == "mine" and int(got[0].number) == 3
    assert str(got[0].pdbItemList[0].structure.fullPath) == str(GAMMA / "gamma_model.pdb")
    assert float(got[0].pdbItemList[0].identity_to_target) == 0.8


@pytest.mark.parametrize("old_choice, expected", [("MODEL", "MODEL"), ("NONE", "NONE"), ("MAP", None)])
def test_ep_partial_model_choice_is_renamed_where_it_exists(tmp_path, old_choice, expected):
    old = make("phaser_EP", tmp_path)
    old.container.inputData.PARTIALMODELORMAP.set(old_choice)
    old.container.inputData.WAVELENGTH.set(1.54)
    new = make("phaser_ep_phil", tmp_path)
    adopted = new.adopt_legacy_container(old.container)
    if expected is None:
        assert "PARTIAL_BY" not in adopted
    else:
        assert str(new.container.inputData.PARTIAL_BY) == expected
    assert float(new.container.inputData.WAVELENGTH) == 1.54
