"""servalcat_pipe validity is mode-aware.

The pipeline (what users run) validated the same way in both modes: it always
recommended a Free R set and never required a resolution. But in SPA mode Free R
is meaningless (cross-validation is via half maps) and a resolution is mandatory
(servalcat SPA builds its grid from d_min). So:

- SPA  -> no Free R recommendation; RES_MIN is a blocking error until set.
- xtal -> Free R recommendation as before; RES_MIN not required.
"""

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.tasks import get_plugin_class


def _reports(task, suffix):
    return [e for e in task.validity().getErrors()
            if e.get("name", "").endswith(suffix)]


def _pipe(method, res_min=None):
    t = get_plugin_class("servalcat_pipe")()
    t.container.controlParameters.DATA_METHOD.set(method)
    if res_min is not None:
        t.container.controlParameters.RES_MIN.set(res_min)
    return t


# --- SPA mode ---------------------------------------------------------------

def test_spa_does_not_recommend_free_r():
    assert _reports(_pipe("spa"), "inputData.FREERFLAG") == []


def test_spa_requires_resolution():
    (report,) = _reports(_pipe("spa"), "controlParameters.RES_MIN")
    assert report["severity"] == CCP4ErrorHandling.SEVERITY_ERROR


def test_spa_with_resolution_is_clean():
    assert _reports(_pipe("spa", res_min=3.2), "controlParameters.RES_MIN") == []


# --- xtal mode --------------------------------------------------------------

def test_xtal_still_recommends_free_r():
    (report,) = _reports(_pipe("xtal"), "inputData.FREERFLAG")
    assert report["severity"] == CCP4ErrorHandling.SEVERITY_WARNING


def test_xtal_does_not_require_resolution():
    assert _reports(_pipe("xtal"), "controlParameters.RES_MIN") == []
