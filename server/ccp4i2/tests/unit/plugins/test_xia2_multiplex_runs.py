"""xia2_multiplex offered "Previous xia2 run directories" but never read them,
so a job given only xia2 runs ran with no data ("No Experiments found"). Each
run's integrated DIALS files (not its scaled ones) now reach the command."""
import pytest

from ccp4i2.core.tasks import get_plugin_class

pytest.importorskip("libtbx.phil", reason="the command builds a PHIL file (needs libtbx, CCP4/cctbx)")


def test_xia2_runs_reach_the_command(tmp_path):
    runs = []
    for name in ("job_3", "job_4"):
        data = tmp_path / name / "DataFiles"
        data.mkdir(parents=True)
        for f in ("SWEEP1.expt", "SWEEP1.refl", "AUTOMATIC_DEFAULT_scaled.expt",
                  "AUTOMATIC_DEFAULT_scaled.refl"):
            (data / f).write_text("")
        runs.append(tmp_path / name)
    work = tmp_path / "work"
    work.mkdir()
    plugin = get_plugin_class("xia2_multiplex")(workDirectory=str(work), name="mx")
    for run in runs:
        plugin.container.inputData.XIA2_RUN.append(plugin.container.inputData.XIA2_RUN.makeItem())
        plugin.container.inputData.XIA2_RUN[-1].setFullPath(str(run))
    plugin.makeCommandAndScript()
    args = [str(a) for a in plugin.commandLine]
    assert [a for a in args if a.startswith("experiments=")] == [
        f"experiments={run / 'DataFiles' / 'SWEEP1.expt'}" for run in runs]
