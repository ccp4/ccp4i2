"""``refreshRemoteReport``: the report of a run happening somewhere else.

A dispatched run's plugin is not the process doing the work, so nothing
rewrote program.xml between submit and harvest and the Report page stayed on
"handed to a run target" for the hours the run took. The output tree is the
progress -- PanDDA makes a directory per dataset and fills each in turn -- so
the method counts the tree and rewrites the report from it. It is a plugin
method, so the generic object_method endpoint reaches it, and reconcile calls
it whenever the user asks how the run is doing.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.lib.utils.jobs import dispatch_record as dr


def _tree(work: Path, analysed: int, loaded: int) -> None:
    """An output tree partway through: ``loaded`` dataset directories, of
    which ``analysed`` have their Z-map (what summarise_output_tree counts)."""
    processed = work / "pandda2_out" / "processed_datasets"
    for i in range(loaded):
        name = f"xtal-{i:04d}"
        (processed / name).mkdir(parents=True)
        if i < analysed:
            (processed / name / f"{name}-z_map.native.ccp4").write_bytes(b"map")


@pytest.fixture
def plugin(tmp_path):
    work = tmp_path / "job_3"
    work.mkdir()
    dr.new_record(work, target="batch", handle="batch-42")
    dr.write_record(work, state="running")
    return get_plugin_class("pandda_campaign")(workDirectory=str(work))


def _report(plugin) -> ET.Element:
    return ET.parse(plugin.makeFileName("PROGRAMXML")).getroot()


def test_counts_the_tree_and_rewrites_the_report(plugin):
    work = Path(plugin.workDirectory)
    _tree(work, analysed=2, loaded=158)

    summary = plugin.refreshRemoteReport()

    assert (summary["processed"], summary["analysed"]) == (158, 2)
    root = _report(plugin)
    assert root.findtext("state") == "dispatched"
    assert root.findtext("n_processed") == "158"
    assert root.findtext("n_analysed") == "2"
    assert root.find("dispatch/handle").text == "batch-42"


def test_reports_progress_again_as_the_run_moves_on(plugin):
    work = Path(plugin.workDirectory)
    _tree(work, analysed=2, loaded=158)
    plugin.refreshRemoteReport()

    for i in range(2, 40):
        name = f"xtal-{i:04d}"
        (work / "pandda2_out" / "processed_datasets" / name / f"{name}-z_map.native.ccp4").write_bytes(b"map")
    plugin.refreshRemoteReport()

    assert _report(plugin).findtext("n_analysed") == "40"


def test_an_empty_tree_is_not_an_error(plugin):
    summary = plugin.refreshRemoteReport()
    assert summary["processed"] == 0
    assert _report(plugin).findtext("state") == "dispatched"


def test_an_unreadable_tree_never_raises(plugin, monkeypatch):
    monkeypatch.setattr(plugin, "_record_tree",
                        lambda: (_ for _ in ()).throw(OSError("share gone")))
    assert plugin.refreshRemoteReport() == {}


def test_the_method_is_reachable_by_name(plugin):
    """object_method calls a plugin method by name when the object path is
    just the task, which is how the interface would offer a refresh."""
    assert callable(getattr(plugin, "refreshRemoteReport", None))
