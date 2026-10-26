"""modelcraft warns when its phases and its model come from different jobs."""
from pathlib import Path

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _advice(error):
    return [r for r in error._reports if r["code"] == 202]


def test_phases_from_another_job_than_the_model_get_advice(tmp_path):
    plugin = get_plugin_class("modelcraft")(workDirectory=str(tmp_path), name="mc")
    inp = plugin.container.inputData
    plugin.container.controlParameters.USE_MODEL_PHASES.set(False)  # the default starts from the model alone
    inp.XYZIN.setFullPath(str(DEMO / "beta.pdb"))
    inp.PHASES.setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    inp.XYZIN.relPath.set("CCP4_JOBS/job_6")
    inp.PHASES.relPath.set("CCP4_JOBS/job_5")
    advice = _advice(plugin.validity())
    assert advice and advice[0]["severity"] == CCP4ErrorHandling.SEVERITY_WARNING
    assert str(advice[0]["name"]).endswith("inputData.PHASES") and "job_5" in advice[0]["details"]
    inp.PHASES.relPath.set("CCP4_JOBS/job_6")
    assert not _advice(plugin.validity())
    inp.PHASES.relPath.set("CCP4_JOBS/job_5")
    plugin.container.controlParameters.USE_MODEL_PHASES.set(True)
    assert not _advice(plugin.validity())


def _no_start(error):
    return [r for r in error._reports if r["code"] == 203]


def test_a_job_with_nothing_to_start_from_cannot_run(tmp_path):
    """Box unticked and no phases, or ticked and no model: ModelCraft has no
    starting phases, so the run is refused here rather than failing later."""
    plugin = get_plugin_class("modelcraft")(workDirectory=str(tmp_path), name="mc")
    inp = plugin.container.inputData
    params = plugin.container.controlParameters

    params.USE_MODEL_PHASES.set(False)
    blocked = _no_start(plugin.validity())
    assert blocked and blocked[0]["severity"] == CCP4ErrorHandling.SEVERITY_ERROR
    assert str(blocked[0]["name"]).endswith("inputData.PHASES")
    inp.PHASES.setFullPath(str(DEMO / "beta_blip_P3221.mtz"))
    assert not _no_start(plugin.validity())

    params.USE_MODEL_PHASES.set(True)
    blocked = _no_start(plugin.validity())
    assert blocked and str(blocked[0]["name"]).endswith("inputData.XYZIN")
    inp.XYZIN.setFullPath(str(DEMO / "beta.pdb"))
    assert not _no_start(plugin.validity())
