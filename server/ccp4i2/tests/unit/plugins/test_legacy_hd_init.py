"""A refinement job from before 2025 says HD_INIT NO; the current def.xml
allows only DEUTERIUM and MIXTURE, with HD_INIT_TOGGLE saying whether. Cloning
such a job used to give "Value 'NO' not in allowed values" and no way to run
it (#860). Loading its parameters now turns the toggle off instead."""
from pathlib import Path

import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.wrappers.refmac.script.refmac import upgrade_legacy_hd_init

LEGACY_PARAMS = """<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <function>PARAMS</function>
    <pluginName>{task}</pluginName>
  </ccp4i2_header>
  <ccp4i2_body id="{task}">
    <controlParameters>
      <HD_INIT>NO</HD_INIT>
    </controlParameters>
  </ccp4i2_body>
</ccp4:ccp4i2>
"""


def _hd_errors(plugin):
    return [r for r in plugin.validity().getErrors()
            if r["severity"] >= CCP4ErrorHandling.SEVERITY_ERROR
            and str(r.get("name", "")).endswith("HD_INIT")]


@pytest.mark.parametrize("task", ["refmac", "prosmart_refmac"])
def test_a_2024_job_loads_with_the_toggle_off(tmp_path, task):
    params = Path(tmp_path) / "input_params.xml"
    params.write_text(LEGACY_PARAMS.format(task=task))
    plugin = get_plugin_class(task)(workDirectory=str(tmp_path), name=task)
    error = plugin.loadDataFromXml(str(params))
    assert error.maxSeverity() < CCP4ErrorHandling.SEVERITY_ERROR, str(error)
    control = plugin.container.controlParameters
    assert str(control.HD_INIT) == "MIXTURE"
    assert bool(control.HD_INIT_TOGGLE) is False
    assert _hd_errors(plugin) == []


def test_current_values_are_left_alone(tmp_path):
    plugin = get_plugin_class("refmac")(workDirectory=str(tmp_path), name="refmac")
    control = plugin.container.controlParameters
    control.HD_INIT.set("DEUTERIUM")
    control.HD_INIT_TOGGLE.set(True)
    assert upgrade_legacy_hd_init(control) is False
    assert str(control.HD_INIT) == "DEUTERIUM"
    assert bool(control.HD_INIT_TOGGLE) is True


def test_a_container_without_the_parameters_is_ignored():
    class Bare:
        pass
    assert upgrade_legacy_hd_init(Bare()) is False
