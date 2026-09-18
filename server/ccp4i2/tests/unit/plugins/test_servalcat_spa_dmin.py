"""SPA refinement requires a resolution (d_min).

servalcat's SPA path (`refine_spa_norefmac`) always puts `-d <RES_MIN>` on the
command line and builds its grid from it, so an unset RES_MIN is not a soft
default: the job dies at runtime with `initialize_grid(): d_min is not set`.
That is the trap the cryo-EM placement -> servalcat handoff walked people into
(CryoMapMR job 9). validity() now turns it into a blocking field error on
RES_MIN, before submission -- and only in SPA mode.
"""

import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class


def _res_min_errors(task):
    return [
        e for e in task.validity().getErrors()
        if e.get("name", "").endswith("controlParameters.RES_MIN")
    ]


def _spa_task(res_min=None):
    task = get_plugin_class("servalcat")()
    task.container.controlParameters.DATA_METHOD.set("spa")
    if res_min is not None:
        task.container.controlParameters.RES_MIN.set(res_min)
    return task


def test_spa_without_res_min_is_a_blocking_error():
    (report,) = _res_min_errors(_spa_task())
    assert report["severity"] == CCP4ErrorHandling.SEVERITY_ERROR
    assert report["name"] == "servalcat.container.controlParameters.RES_MIN"


def test_setting_res_min_clears_it():
    assert _res_min_errors(_spa_task(res_min=3.2)) == []


def test_xtal_mode_does_not_require_res_min():
    task = get_plugin_class("servalcat")()
    task.container.controlParameters.DATA_METHOD.set("xtal")
    assert _res_min_errors(task) == []
