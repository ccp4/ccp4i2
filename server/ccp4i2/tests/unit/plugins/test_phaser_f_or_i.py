"""
F_OR_I is settled against the data, not the GUI.

F_OR_I defaults to I. Given a file of amplitudes --- from which intensities
cannot be made --- every phaser MR wrapper asked makeHklin for IMEAN and the
job failed. Whether the file was picked by hand or autopopulated from an
earlier job, the GUI need not have fired anything, so the pipelines and the
wrapper base resolve the flag themselves and record the result.
"""

import pytest

from ccp4i2 import I2_TOP
from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR.script import phaser_MR


@pytest.fixture(scope="module")
def fmean_mtz(tmp_path_factory):
    """beta_blip's amplitudes, split to a canonical F,SIGF mini-MTZ."""
    pytest.importorskip("gemmi")
    from ccp4i2.core.CCP4Utils import split_mtz_file

    out = tmp_path_factory.mktemp("obs") / "beta_blip_F_SIGF.mtz"
    split_mtz_file(
        I2_TOP / "demo_data" / "beta_blip" / "beta_blip_P3221.mtz",
        out,
        CCP4XtalData.CObsDataFile.canonical_column_mapping(4, ["Fobs", "Sigma"]),
    )
    return str(out)


IMEAN_MTZ = str(I2_TOP / "demo_data" / "gamma" / "merged_intensities_Xe.mtz")

TASKS = ["phaser_pipeline", "phaser_simple", "phaser_MR_AUTO", "phaser_MR_FRF",
         "phaser_MR_FTF", "phaser_MR_RNP", "phaser_MR_PAK"]


def _task(name, obs_path, f_or_i=None):
    task = get_plugin_class(name)()
    task.container.inputData.F_SIGF.set(obs_path)
    if f_or_i is not None:
        task.container.inputData.F_OR_I.set(f_or_i)
    return task


def test_the_default_is_i():
    task = get_plugin_class("phaser_pipeline")()
    assert str(task.container.inputData.F_OR_I) == "I"
    assert not task.container.inputData.F_OR_I.isSet(allowDefault=False)


@pytest.mark.parametrize("name", TASKS)
def test_amplitudes_only_forces_f(name, fmean_mtz):
    task = _task(name, fmean_mtz)
    assert phaser_MR.resolve_f_or_i(task.container.inputData) == "F"
    assert str(task.container.inputData.F_OR_I) == "F"
    # ...and it is recorded, so the params file says what was used
    assert task.container.inputData.F_OR_I.isSet(allowDefault=False)


@pytest.mark.parametrize("name", ["phaser_pipeline", "phaser_MR_AUTO"])
def test_intensities_keep_the_default_i(name):
    task = _task(name, IMEAN_MTZ)
    assert phaser_MR.resolve_f_or_i(task.container.inputData) == "I"
    assert str(task.container.inputData.F_OR_I) == "I"


def test_a_choice_of_f_for_intensities_is_respected():
    task = _task("phaser_pipeline", IMEAN_MTZ, f_or_i="F")
    assert phaser_MR.resolve_f_or_i(task.container.inputData) == "F"


def test_no_file_changes_nothing():
    task = get_plugin_class("phaser_pipeline")()
    assert phaser_MR.resolve_f_or_i(task.container.inputData) == "I"
    assert not task.container.inputData.F_OR_I.isSet(allowDefault=False)


def test_content_flag_is_introspected_when_not_recorded(fmean_mtz):
    task = _task("phaser_MR_AUTO", fmean_mtz)
    task.container.inputData.F_SIGF.contentFlag.unSet()
    assert phaser_MR.resolve_f_or_i(task.container.inputData) == "F"


# --- validity: advisory, on the field the user can change ---------------------

def _f_or_i_reports(task):
    return [e for e in task.validity().getErrors()
            if e.get("name", "").endswith("inputData.F_OR_I")]


@pytest.mark.parametrize("name", ["phaser_pipeline", "phaser_simple", "phaser_MR_AUTO"])
def test_explicit_i_on_amplitudes_is_a_warning_not_an_error(name, fmean_mtz):
    task = _task(name, fmean_mtz, f_or_i="I")
    (report,) = _f_or_i_reports(task)
    assert report["severity"] == CCP4ErrorHandling.SEVERITY_WARNING
    assert report["name"] == f"{name}.container.inputData.F_OR_I"
    assert "amplitudes" in report["details"]


def test_the_default_i_on_amplitudes_is_reported_too(fmean_mtz):
    """The field is visible beside the reflections; the user can act on it."""
    (report,) = _f_or_i_reports(_task("phaser_pipeline", fmean_mtz))
    assert report["severity"] == CCP4ErrorHandling.SEVERITY_WARNING


def test_choosing_f_clears_the_warning(fmean_mtz):
    assert _f_or_i_reports(_task("phaser_pipeline", fmean_mtz, f_or_i="F")) == []


def test_the_warning_reaches_run_time_validity(fmean_mtz):
    task = _task("phaser_pipeline", fmean_mtz)
    names = [e.get("name", "") for e in task.runTimeValidity().getErrors()]
    assert "phaser_pipeline.container.inputData.F_OR_I" in names


def test_explicit_i_on_intensities_is_not_reported():
    assert _f_or_i_reports(_task("phaser_pipeline", IMEAN_MTZ, f_or_i="I")) == []
