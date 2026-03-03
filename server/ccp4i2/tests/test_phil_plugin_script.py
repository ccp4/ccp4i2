"""
Tests for PhilPluginScript — base class for PHIL-based tool wrappers.

Uses a mock subclass with synthetic PHIL to avoid dependency on phasertng.
"""

import os
import pytest
import tempfile
from pathlib import Path

from libtbx.phil import parse

from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.core.base_object.base_classes import CContainer, ValueState
from ccp4i2.core.base_object.fundamental_types import CInt, CFloat, CBoolean, CString


# ---------------------------------------------------------------------------
# Mock PhilPluginScript subclass with synthetic PHIL
# ---------------------------------------------------------------------------

MOCK_PHIL = parse("""
    refinement {
      resolution = 2.0
        .short_caption = "Resolution limit"
        .type = float(value_min=0.5, value_max=20.0)
      cycles = 10
        .type = int
      use_solvent = True
        .type = bool
      mode = *auto manual
        .type = choice
    }
    advanced {
      weight = 0.5
        .type = float
      strategy = *fast thorough
        .type = choice
    }
""")


class MockPhilPlugin(PhilPluginScript):
    TASKNAME = "mock_phil_plugin"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return MOCK_PHIL

    def get_phil_exclude_scopes(self):
        return []

    def get_command_target(self):
        return "echo"


class MockPhilPluginWithExclusions(PhilPluginScript):
    TASKNAME = "mock_phil_excluded"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return MOCK_PHIL

    def get_phil_exclude_scopes(self):
        return ["advanced"]

    def get_command_target(self):
        return "echo"


# ---------------------------------------------------------------------------
# Instantiation and parameter merge
# ---------------------------------------------------------------------------

class TestInstantiation:

    def test_creates_standard_containers(self):
        plugin = MockPhilPlugin()
        assert plugin.container is not None
        assert plugin.container.inputData is not None
        assert plugin.container.outputData is not None
        assert plugin.container.controlParameters is not None

    def test_phil_parameters_merged_into_control_parameters(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters
        children = list(cp.children())
        # Should have at least the 2 scope containers (refinement, advanced)
        assert len(children) >= 2

    def test_phil_leaf_types(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters

        # Navigate into refinement scope
        ref = cp.refinement
        assert isinstance(ref, CContainer)

        resolution = ref.refinement__resolution
        assert isinstance(resolution, CFloat)
        assert resolution.get() == 2.0

        cycles = ref.refinement__cycles
        assert isinstance(cycles, CInt)
        assert cycles.get() == 10

        use_solvent = ref.refinement__use_solvent
        assert isinstance(use_solvent, CBoolean)
        assert use_solvent.get() is True

        mode = ref.refinement__mode
        assert isinstance(mode, CString)
        assert mode.get() == "auto"

    def test_phil_path_qualifier_preserved(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters
        resolution = cp.refinement.refinement__resolution
        assert resolution.get_qualifier("philPath") == "refinement.resolution"

    def test_exclusions_applied(self):
        plugin = MockPhilPluginWithExclusions()
        cp = plugin.container.controlParameters
        # refinement should be present
        assert hasattr(cp, "refinement")
        # advanced should NOT be present
        assert not hasattr(cp, "advanced")


# ---------------------------------------------------------------------------
# Extract parameters
# ---------------------------------------------------------------------------

class TestExtractParameters:

    def test_no_user_modified_params_initially(self):
        plugin = MockPhilPlugin()
        params = plugin.extract_phil_parameters()
        assert params == [], \
            f"Expected no user-modified params, got {params}"

    def test_user_change_is_extracted(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters

        # Simulate user changing resolution
        cp.refinement.refinement__resolution.value = 3.0

        params = plugin.extract_phil_parameters()
        assert len(params) == 1
        assert params[0][0] == "refinement.resolution"
        assert "3.0" in params[0][1]

    def test_multiple_user_changes(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters

        # Change multiple parameters
        cp.refinement.refinement__resolution.value = 3.0
        cp.refinement.refinement__cycles.value = 20
        cp.advanced.advanced__weight.value = 0.8

        params = plugin.extract_phil_parameters()
        paths = [p[0] for p in params]
        assert "refinement.resolution" in paths
        assert "refinement.cycles" in paths
        assert "advanced.weight" in paths
        assert len(params) == 3

    def test_unchanged_defaults_not_extracted(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters

        # Only change one parameter
        cp.refinement.refinement__resolution.value = 3.0

        params = plugin.extract_phil_parameters()
        paths = [p[0] for p in params]

        # These should NOT appear (still at defaults)
        assert "refinement.cycles" not in paths
        assert "refinement.use_solvent" not in paths
        assert "refinement.mode" not in paths
        assert "advanced.weight" not in paths


# ---------------------------------------------------------------------------
# build_working_phil
# ---------------------------------------------------------------------------

class TestBuildWorkingPhil:

    def test_generates_working_phil_file(self):
        plugin = MockPhilPlugin()

        with tempfile.TemporaryDirectory() as tmpdir:
            plugin._workDirectory = tmpdir
            phil_path = plugin.build_working_phil()

            assert os.path.exists(phil_path)
            assert phil_path.endswith("working.phil")

            with open(phil_path) as f:
                content = f.read()
            assert len(content) > 0

    def test_working_phil_contains_user_change(self):
        plugin = MockPhilPlugin()
        cp = plugin.container.controlParameters

        # Set a non-default value
        cp.refinement.refinement__resolution.value = 3.5

        with tempfile.TemporaryDirectory() as tmpdir:
            plugin._workDirectory = tmpdir
            phil_path = plugin.build_working_phil()

            with open(phil_path) as f:
                content = f.read()

            # The working.phil should contain the user's value
            assert "3.5" in content

    def test_working_phil_contains_defaults(self):
        plugin = MockPhilPlugin()

        with tempfile.TemporaryDirectory() as tmpdir:
            plugin._workDirectory = tmpdir
            phil_path = plugin.build_working_phil()

            with open(phil_path) as f:
                content = f.read()

            # Default values should still appear
            assert "cycles = 10" in content
            assert "use_solvent = True" in content


# ---------------------------------------------------------------------------
# Parameter round-trip (set → save → reload)
# ---------------------------------------------------------------------------

class TestRoundTrip:

    def test_set_save_reload_preserves_value(self):
        """Set a parameter, save to XML, reload, verify the value persists."""
        plugin1 = MockPhilPlugin()

        # Set a value
        plugin1.container.controlParameters.refinement.refinement__resolution.value = 4.2

        # Verify it's set
        assert plugin1.container.controlParameters.refinement.refinement__resolution.get() == 4.2

        # Save to XML
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "input_params.xml")
            plugin1.saveDataToXml(xml_path)

            assert os.path.exists(xml_path)

            # Reload into a fresh plugin
            plugin2 = MockPhilPlugin(xmlFile=xml_path)

            # Verify the value survived the round-trip
            val = plugin2.container.controlParameters.refinement.refinement__resolution.get()
            assert val == pytest.approx(4.2), \
                f"Expected 4.2 after round-trip, got {val}"


# ---------------------------------------------------------------------------
# Subclass contract
# ---------------------------------------------------------------------------

class TestSubclassContract:

    def test_get_master_phil_not_implemented(self):
        """PhilPluginScript.get_master_phil() should raise if not overridden."""

        class BarePlugin(PhilPluginScript):
            TASKNAME = "bare_plugin"

        with pytest.raises(NotImplementedError):
            BarePlugin().get_master_phil()

    def test_get_command_target_not_implemented(self):
        """PhilPluginScript.get_command_target() should raise if not overridden."""

        class BarePlugin(PhilPluginScript):
            TASKNAME = "bare_plugin"

            def get_master_phil(self):
                return parse("x = 1\n  .type = int")

        plugin = BarePlugin()
        with pytest.raises(NotImplementedError):
            plugin.get_command_target()

    def test_default_shim_definitions_empty(self):
        plugin = MockPhilPlugin()
        assert plugin.get_shim_definitions() == []

    def test_default_exclude_scopes_from_class_var(self):
        """PHIL_EXCLUDE_SCOPES class variable is used by default."""

        class ExcludePlugin(PhilPluginScript):
            TASKNAME = "exclude_plugin"
            PHIL_EXCLUDE_SCOPES = ["output", "debug"]

            def get_master_phil(self):
                return parse("x = 1\n  .type = int")

            def get_command_target(self):
                return "echo"

        plugin = ExcludePlugin()
        assert plugin.get_phil_exclude_scopes() == ["output", "debug"]
