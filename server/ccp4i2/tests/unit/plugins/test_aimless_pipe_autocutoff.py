"""aimless_pipe's automatic resolution cutoff must do what it says.

The cutoff is estimated by phaser_analysis after a first Aimless run and
applied by a second. Several things went wrong with that, all silent, and
all ending in data scaled to a resolution nobody chose.

The root cause is one shape: this pipeline reports its verdict from deep
inside a chain of ordinary calls, and process_finish returns like any other
function. So "exit the pipeline" did not exit anything -- it unwound back
into the loop that called it, which carried on.

* A first run that had already finished the job was followed by a second
  run of Aimless, phaser_analysis and ctruncate, whose results replaced it.
* The branch for "the data already reach the edge, so no cutoff is needed"
  fell through into the code that sets the cutoff, overwriting both flags
  it had just set.
* An explicit resolution range reached the Aimless job only on the
  no-cutoff path: with a cutoff, the high-resolution end was set on a
  container nobody had copied into, so a low-resolution limit was dropped.

And asking for a cutoff with the Phaser analysis switched off cannot work
at all -- the pipeline dropped AUTOCUTOFF to a print(), leaving the user to
believe their data had been cut where the data stop.

These drive the real functions with only the subprocess-running parts
stubbed, so they fail if the fixes are reverted.

CCP4-free: no binary is run.
"""
import pytest
from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.error_reporting import Severity
from ccp4i2.core.tasks import get_plugin_class

# Enough of an Aimless XML for checkaimlessresult: healthy inner shell, so
# the data are not a "complete disaster" and the run proceeds.
HEALTHY_AIMLESS_XML = """<AIMLESS><Result><Dataset name="p/x/d">
  <CChalf><Inner>0.977</Inner></CChalf>
  <Multiplicity><Inner>3.4</Inner></Multiplicity>
</Dataset></Result></AIMLESS>"""


def _pipeline(tmp_path, **control):
    plugin = get_plugin_class("aimless_pipe")(
        workDirectory=str(tmp_path), name="aimless_pipe")
    for name, value in control.items():
        getattr(plugin.container.controlParameters, name).set(value)
    # What process() sets up before the first sub-job runs.
    plugin.fatalError = None
    plugin._pipelineFinished = False
    plugin.rootXML = etree.Element("AIMLESS_PIPE")
    plugin.aimless1xml = None
    plugin.phaser_analysisxml = None
    plugin.phaser_analysis1xml = None
    return plugin


def _subplugin(task, tmp_path, name):
    directory = tmp_path / name
    directory.mkdir(exist_ok=True)
    plugin = get_plugin_class(task)(workDirectory=str(directory), name=name)
    plugin.process = lambda: CPluginScript.SUCCEEDED
    return plugin


def _codes(report, severity=None):
    return [e.get("code") for e in report.getErrors()
            if severity is None or e["severity"] == severity]


# --------------------------------------------------------------------------
# Asking for a cutoff that cannot be estimated
# --------------------------------------------------------------------------

def test_cutoff_without_phaser_analysis_says_so(tmp_path):
    plugin = _pipeline(tmp_path, AUTOCUTOFF=True, DOPHASERANALYSIS=False)
    assert 203 in _codes(plugin.validity(), Severity.WARNING)


def test_cutoff_without_phaser_analysis_does_not_block_the_run(tmp_path):
    """Advisory: the run is still valid, it just will not cut the data."""
    plugin = _pipeline(tmp_path, AUTOCUTOFF=True, DOPHASERANALYSIS=False)
    assert 203 not in _codes(plugin.validity(), Severity.ERROR)


@pytest.mark.parametrize("autocutoff,phaser", [
    (True, True),     # asked for it, and it can be estimated
    (False, False),   # did not ask for it
    (False, True),    # the default
])
def test_no_complaint_about_any_other_combination(tmp_path, autocutoff, phaser):
    plugin = _pipeline(tmp_path, AUTOCUTOFF=autocutoff, DOPHASERANALYSIS=phaser)
    assert 203 not in _codes(plugin.validity())


# --------------------------------------------------------------------------
# Finishing the pipeline ends it
# --------------------------------------------------------------------------

def test_a_finished_pipeline_does_not_run_aimless_again(tmp_path):
    """The headline symptom: everything ran twice and the second pass won.

    Drives the real run loop in process_cycle_aimless, with process_aimless
    standing in for a first run that finished the pipeline -- what happens
    when the data already reach the edge, or are declared hopeless.
    """
    plugin = _pipeline(tmp_path, AUTOCUTOFF=True, DOPHASERANALYSIS=True)
    plugin.doPhaserAnalysis = True
    runs = []

    def finish_on_first_run():
        # What process_finish does, without its FreeR and report work --
        # that it sets the flag once and only once is a separate test.
        runs.append(plugin.aimlessruncount)
        plugin._pipelineFinished = True

    plugin.process_aimless = finish_on_first_run

    plugin.process_cycle_aimless(CPluginScript.SUCCEEDED)

    assert runs == [1], "Aimless ran again after the pipeline had finished"


def test_process_finish_reports_once(tmp_path):
    """The backstop under the returns: a second verdict is refused."""
    plugin = _pipeline(tmp_path)
    reported = []
    plugin.reportStatus = reported.append
    plugin._pipelineFinished = True

    plugin.process_finish(CPluginScript.SUCCEEDED)

    assert reported == []


def test_data_reaching_the_edge_are_not_then_cut(tmp_path):
    """analyseResolution says no cutoff is needed; the flags must survive.

    Falling through overwrote both of them, and the second Aimless run then
    applied a cutoff that this branch exists to say was unnecessary.
    """
    plugin = _pipeline(tmp_path, AUTOCUTOFF=True, DOPHASERANALYSIS=True)
    plugin.doPhaserAnalysis = True
    plugin.AUTOCUTOFF = True
    plugin.aimlessruncount = 1
    plugin.aimless = _subplugin("aimless", tmp_path, "aimless_1")
    with open(plugin.aimless.makeFileName("PROGRAMXML"), "w") as handle:
        handle.write(HEALTHY_AIMLESS_XML)

    # Data reach the edge: no cutoff wanted, and 2.0 A is what it would have
    # been cut to had one been wanted.
    plugin.analyseResolution = lambda: (True, 2.0)
    plugin.process_cycle_ctruncate = lambda: None

    plugin.process_post_aimless(CPluginScript.SUCCEEDED)

    assert plugin.highrescutoff == -1.0, "a cutoff was recorded anyway"
    assert plugin.cutoffdone is False


# --------------------------------------------------------------------------
# What the second Aimless run is told about resolution
# --------------------------------------------------------------------------

def _resolution_range_given_to_aimless(plugin, tmp_path, autocutoff, cutoff):
    """Run process_aimless's second-run path and report what Aimless got."""
    aimless = _subplugin("aimless", tmp_path, "aimless_2")
    plugin.pointless = _subplugin("pointless", tmp_path, "pointless")
    plugin.makePluginObject = lambda task: aimless
    plugin.process_post_aimless = lambda status: None
    plugin.AUTOCUTOFF = autocutoff
    plugin.highrescutoff = cutoff
    plugin.aimlessruncount = 2

    plugin.process_aimless()

    return aimless.container.controlParameters.RESOLUTION_RANGE


def test_an_explicit_low_resolution_limit_survives_the_cutoff(tmp_path):
    """Asking for a range and a cutoff together used to lose the range."""
    plugin = _pipeline(tmp_path, AUTOCUTOFF=True)
    plugin.container.controlParameters.RESOLUTION_RANGE.start.set(30.0)
    plugin.container.controlParameters.RESOLUTION_RANGE.end.set(1.8)

    resrange = _resolution_range_given_to_aimless(plugin, tmp_path, True, 2.4)

    assert float(resrange.start) == 30.0, "low-resolution limit dropped"
    assert float(resrange.end) == 2.4, "estimated cutoff not applied"


def test_without_a_cutoff_the_users_range_is_used_unchanged(tmp_path):
    plugin = _pipeline(tmp_path)
    plugin.container.controlParameters.RESOLUTION_RANGE.start.set(30.0)
    plugin.container.controlParameters.RESOLUTION_RANGE.end.set(1.8)

    resrange = _resolution_range_given_to_aimless(plugin, tmp_path, False, -1.0)

    assert float(resrange.start) == 30.0
    assert float(resrange.end) == 1.8
