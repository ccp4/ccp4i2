"""mrbump_basic stops before MrBUMP when the AU contents give no sequence.

Given an empty sequence, MrBUMP does not stop: it spun at full CPU for as
long as it was left, writing nothing, the job Running, when an agent gave
ASUIN a .seq file, which does not parse as AU contents (2026-10-12).
"""
from pathlib import Path

import ccp4i2
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import CErrorReport
from ccp4i2.core.tasks import get_plugin_class

BETA_SEQ = Path(ccp4i2.__file__).parent / "demo_data" / "beta_blip" / "beta.seq"


def _plugin(tmp_path, monkeypatch):
    plugin = get_plugin_class("mrbump_basic")(workDirectory=str(tmp_path), name="m")
    monkeypatch.setattr(plugin, "makeHklin", lambda objects: (str(tmp_path / "hklin.mtz"), CErrorReport()))
    return plugin


def test_an_asu_file_that_is_not_one_stops_the_job(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, monkeypatch)
    plugin.container.inputData.ASUIN.setFullPath(str(BETA_SEQ))
    assert plugin.processInputFiles() == CPluginScript.FAILED
    assert 201 in [r["code"] for r in plugin.errorReport._reports]


def test_no_asu_file_at_all_stops_it_too(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, monkeypatch)
    assert plugin.processInputFiles() == CPluginScript.FAILED


def test_real_au_contents_go_through(tmp_path, monkeypatch):
    plugin = _plugin(tmp_path, monkeypatch)
    plugin.container.inputData.ASUIN.setFullPath(
        str(Path(ccp4i2.__file__).parent / "demo_data" / "gamma" / "gamma.asu.xml"))
    assert plugin.processInputFiles() == CPluginScript.SUCCEEDED
    assert ">" in Path(plugin.seqFile).read_text()
