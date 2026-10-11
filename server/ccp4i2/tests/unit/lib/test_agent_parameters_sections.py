"""A job's parameters include every section the task declares."""
import tempfile

from ccp4i2.agent.parameters import summarise
from ccp4i2.core.tasks import get_plugin_class


def test_a_tasks_own_section_is_listed():
    # slicendice keeps MIN_SPLITS in modelParameters; a fixed list of three
    # sections hid it, so set_parameter reported "set: false" though it had
    # taken and a query found nothing (an agent run, 2026-10-11).
    plugin = get_plugin_class("slicendice")(workDirectory=tempfile.mkdtemp(), name="s")
    paths = {e["path"] for e in summarise(plugin.container)}
    assert "modelParameters.MIN_SPLITS" in paths and "inputData.XYZIN" in paths
    assert {e["path"] for e in summarise(plugin.container, sections=("modelParameters",))} <= \
        {p for p in paths if p.startswith("modelParameters.")}
    assert [e["path"] for e in summarise(plugin.container, query="split")] == \
        ["modelParameters.MIN_SPLITS", "modelParameters.MAX_SPLITS"]
