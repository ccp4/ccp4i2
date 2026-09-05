"""A Phaser failure reaches the user as an error, in Phaser's words.

appendErrorReport() records a warning unless told otherwise, and the job
framework then reports "the most recent thing it recorded was only a
warning" -- true, and useless. The two failure paths, a Result that says
no and a Sorry from Phaser's driver, both record code 202 at ERROR
severity with Phaser's own sentence, and leave program.xml behind so the
report can show how far the run got.
"""
import os
import tempfile
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_ERROR
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.wrappers.phaser_phil.script import phaser_run


class FailedResult:
    def Success(self):
        return False

    def ErrorName(self):
        return "INPUT"

    def ErrorMessage(self):
        return "Molecular scattering of beta.pdb deviates more than 20% from the mean."


@pytest.mark.parametrize("task", ["phaser_mr_auto_phil", "phaser_ep_auto_phil"])
class TestFailures:

    def plugin(self, tmp, task):
        plugin = get_plugin_class(task)(workDirectory=tmp, name=task)
        plugin._phil_path = os.path.join(tmp, "working.phil")
        return plugin

    def reports(self, plugin):
        return [(r["code"], r["severity"], str(r["details"])) for r in plugin.errorReport._reports]

    def test_a_result_that_says_no_is_an_error(self, task, monkeypatch):
        monkeypatch.setattr(phaser_run, "run_mode", lambda *a, **k: FailedResult())
        with tempfile.TemporaryDirectory() as tmp:
            plugin = self.plugin(tmp, task)
            assert plugin.startProcess() == CPluginScript.FAILED
            (code, severity, details), = self.reports(plugin)
            assert code == 202 and severity == SEVERITY_ERROR
            assert details.startswith("INPUT: Molecular scattering")
            assert ET.parse(os.path.join(tmp, "program.xml")).getroot().tag.startswith("Phaser")

    def test_a_sorry_from_the_driver_is_phasers_sentence(self, task, monkeypatch):
        def refuse(*a, **k):
            raise phaser_run.PhaserInputError("Couldn't find array F,SIGF in file x.mtz.")
        monkeypatch.setattr(phaser_run, "run_mode", refuse)
        with tempfile.TemporaryDirectory() as tmp:
            plugin = self.plugin(tmp, task)
            assert plugin.startProcess() == CPluginScript.FAILED
            (code, severity, details), = self.reports(plugin)
            assert code == 202 and severity == SEVERITY_ERROR
            assert details == "Couldn't find array F,SIGF in file x.mtz."
            assert "Traceback" not in details


@pytest.mark.parametrize("task", ["phaser_mr_auto_phil", "phaser_ep_auto_phil"])
def test_unpreparable_reflections_are_an_error(task):
    """F_SIGF that cannot be rewritten (here: not set at all) fails the job
    with code 204 at ERROR severity -- and does not raise on the way."""
    with tempfile.TemporaryDirectory() as tmp:
        plugin = get_plugin_class(task)(workDirectory=tmp, name=task)
        assert plugin.processInputFiles() == CPluginScript.FAILED
        codes = [(r["code"], r["severity"]) for r in plugin.errorReport._reports]
        assert (204, SEVERITY_ERROR) in codes


class _Known:
    def __init__(self, modlid):
        self.MODLID = modlid


class _Set:
    def __init__(self, *modlids):
        self.KNOWN = [_Known(m) for m in modlids]


def test_solutions_must_name_the_search_models():
    """A classic Phaser job names its ensembles Ensemble_0...; solutions from
    it given to a job whose search model is called SearchModel fail in Phaser
    with "No model for ensemble Ensemble_0". runTimeValidity says so first."""
    import pickle
    with tempfile.TemporaryDirectory() as tmp:
        plugin = get_plugin_class("phaser_mr_auto_phil")(workDirectory=tmp, name="t")
        inp = plugin.container.inputData
        inp.F_SIGF.setFullPath(os.path.join(tmp, "x.mtz")); open(os.path.join(tmp, "x.mtz"), "w").close()
        e = inp.ENSEMBLES.makeItem(); inp.ENSEMBLES.append(e); e.label.set("SearchModel"); e.number.set(1); e.use.set(True)
        item = e.pdbItemList.makeItem(); e.pdbItemList.append(item)
        item.structure.setFullPath(os.path.join(tmp, "m.pdb")); open(os.path.join(tmp, "m.pdb"), "w").close()
        item.identity_to_target.set(0.9)
        path = os.path.join(tmp, "SOLOUT.phaser_sol.pkl")
        with open(path, "wb") as f:
            pickle.dump([_Set("Ensemble_0")], f)
        inp.SOLIN.setFullPath(path)
        reports = [(r["code"], r["severity"], str(r["details"])) for r in plugin.runTimeValidity()._reports]
        assert any(code == 117 and sev == SEVERITY_ERROR and "Ensemble_0" in d and "SearchModel" in d
                   for code, sev, d in reports), reports
        with open(path, "wb") as f:
            pickle.dump([_Set("SearchModel")], f)
        assert not any(r["code"] == 117 for r in plugin.runTimeValidity()._reports)
