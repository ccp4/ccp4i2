"""A campaign's shared FreeR set has to survive the cell check.

SubstituteLigand hands an input FreeR set to aimless_pipe, which extends it
onto the newly merged data --- but only when the two cells agree to within
Clipper's 1 A tolerance. On a fragment campaign the FreeR set comes from one
reference crystal and is extended onto every other crystal, whose cells drift
by a percent or so; on a 185 A axis that is already outside the tolerance.
The whole CDK4/CyclinD1 pre-screen of September 2026 failed this way, 23 of
26 crystals, each with 'Aimless did not produce FreeR output' and no hint that
the cell test was the reason.

aimless_pipe already has OVERRIDE_CELL_DIFFERENCE for exactly this. These
tests cover its exposure on SubstituteLigand: off by default so the desktop
behaviour is unchanged, and passed through to aimless_pipe when set.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class


def _plugin(tmp_path):
    plugin_class = get_plugin_class("SubstituteLigand")
    assert plugin_class is not None, "SubstituteLigand is not in the registry"
    directory = tmp_path / "SubstituteLigand"
    directory.mkdir(exist_ok=True)
    return plugin_class(workDirectory=str(directory), name="SubstituteLigand")


def _freer(tmp_path, plugin):
    path = tmp_path / "campaign_freer.mtz"
    path.write_bytes(b"MTZ not really, but it exists\n")
    plugin.container.inputData.FREERFLAG_IN.setFullPath(str(path))
    return path


def test_the_override_is_off_by_default(tmp_path):
    plugin = _plugin(tmp_path)
    assert not plugin.container.controlParameters.OVERRIDE_CELL_DIFFERENCE


def test_aimless_keeps_its_cell_check_unless_asked(tmp_path):
    plugin = _plugin(tmp_path)
    _freer(tmp_path, plugin)
    aimless = plugin.makePluginObject("aimless_pipe")
    plugin._configureAimless(aimless)
    assert aimless.container.inputData.FREERFLAG.isSet()
    assert not aimless.container.controlParameters.OVERRIDE_CELL_DIFFERENCE


def test_the_override_reaches_aimless(tmp_path):
    plugin = _plugin(tmp_path)
    _freer(tmp_path, plugin)
    plugin.container.controlParameters.OVERRIDE_CELL_DIFFERENCE.set(True)
    aimless = plugin.makePluginObject("aimless_pipe")
    plugin._configureAimless(aimless)
    assert aimless.container.controlParameters.OVERRIDE_CELL_DIFFERENCE


def test_no_free_r_set_means_nothing_to_override(tmp_path):
    plugin = _plugin(tmp_path)
    plugin.container.controlParameters.OVERRIDE_CELL_DIFFERENCE.set(True)
    aimless = plugin.makePluginObject("aimless_pipe")
    plugin._configureAimless(aimless)
    assert not aimless.container.inputData.FREERFLAG.isSet()
    assert not aimless.container.controlParameters.OVERRIDE_CELL_DIFFERENCE


def test_aimless_still_matches_to_the_reference_model(tmp_path):
    """The rest of the configuration is untouched by the extraction."""
    plugin = _plugin(tmp_path)
    aimless = plugin.makePluginObject("aimless_pipe")
    plugin._configureAimless(aimless)
    cp = aimless.container.controlParameters
    assert str(cp.MODE) == "MATCH"
    assert str(cp.REFERENCE_DATASET) == "XYZ"
    assert cp.AUTOCUTOFF
    assert float(cp.TOLERANCE) == pytest.approx(10.0)
