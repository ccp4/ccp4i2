"""
Tests for PhilPluginScript — base class for PHIL-based tool wrappers.

Uses a mock subclass with synthetic PHIL to avoid dependency on phasertng.
"""

import os
import pytest
import tempfile
from pathlib import Path

# libtbx ships with cctbx inside the CCP4 suite; there is no pip wheel for it.
# Skip rather than error at collection so the suite still runs CCP4-free in CI.
parse = pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)").parse

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
            plugin.workDirectory = tmpdir
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
            plugin.workDirectory = tmpdir
            phil_path = plugin.build_working_phil()

            with open(phil_path) as f:
                content = f.read()

            # The working.phil should contain the user's value
            assert "3.5" in content

    def test_working_phil_contains_defaults(self):
        plugin = MockPhilPlugin()

        with tempfile.TemporaryDirectory() as tmpdir:
            plugin.workDirectory = tmpdir
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


# ---------------------------------------------------------------------------
# Repeated scopes and definitions in the working phil
# ---------------------------------------------------------------------------

MULTI_PHIL = parse("""
    composition {
      solvent = None
        .type = float
      chain
        .optional = True
        .multiple = True
      {
        chain_type = *protein na
          .type = choice
        nres = None
          .type = int
        num = 1
          .type = int
        dataset
          .multiple = True
        {
          label = None
            .type = str
        }
      }
    }
    model = None
      .type = path
      .multiple = True
""")


class MockMultiPlugin(PhilPluginScript):
    TASKNAME = "mock_multi_plugin"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return MULTI_PHIL

    def get_phil_exclude_scopes(self):
        return []

    def get_command_target(self):
        return "echo"


def _fetched(plugin):
    """What the tool would see: working.phil fetched against the master.

    A bare parse of the file carries no .multiple attributes, so it must be
    fetched against the master to read as the tool reads it."""
    with tempfile.TemporaryDirectory() as tmpdir:
        plugin.workDirectory = tmpdir
        phil_path = plugin.build_working_phil()
        return MULTI_PHIL.fetch(sources=[parse(file_name=phil_path)]).extract()


class TestRepeatedScopes:

    def test_no_items_means_no_instances(self):
        plugin = MockMultiPlugin()
        assert plugin.extract_phil_lines() == []
        assert len(_fetched(plugin).composition.chain) == 0
        assert _fetched(plugin).model == []

    def test_two_items_become_two_blocks(self):
        plugin = MockMultiPlugin()
        chain = plugin.container.controlParameters.composition.composition__chain
        chain.append(chain.makeItem())
        chain.append(chain.makeItem())
        chain[0].composition__chain__nres.value = 120
        chain[0].composition__chain__num.value = 2
        chain[1].composition__chain__chain_type.value = "na"

        lines = plugin.extract_phil_lines()
        assert lines[0] == "composition.chain {"
        assert "  nres = 120" in lines and "  num = 2" in lines
        assert lines.count("composition.chain {") == 2

        chains = _fetched(plugin).composition.chain
        assert [(c.nres, c.num, c.chain_type) for c in chains] == [
            (120, 2, "protein"), (None, 1, "na")]

    def test_an_item_identical_to_the_defaults_is_not_an_instance(self):
        # libtbx semantics, not ours: fetch() folds a repeated-scope instance
        # that equals the master template back into the template, so PHIL
        # cannot express "a repeat with all the defaults". The block is
        # written (membership is explicit on our side) and the tool sees none.
        plugin = MockMultiPlugin()
        chain = plugin.container.controlParameters.composition.composition__chain
        chain.append(chain.makeItem())
        assert plugin.extract_phil_lines() == ["composition.chain {", "}"]
        assert len(_fetched(plugin).composition.chain) == 0

    def test_nested_repeated_scope_renders_nested_blocks(self):
        plugin = MockMultiPlugin()
        chain = plugin.container.controlParameters.composition.composition__chain
        chain.append(chain.makeItem())
        chain[0].composition__chain__nres.value = 120
        datasets = chain[0].composition__chain__dataset
        datasets.append(datasets.makeItem())
        datasets.append(datasets.makeItem())
        datasets[0].composition__chain__dataset__label.value = "peak"
        datasets[1].composition__chain__dataset__label.value = "remote"
        assert plugin.extract_phil_lines() == [
            "composition.chain {", "  nres = 120",
            "  dataset {", "    label = peak", "  }",
            "  dataset {", "    label = remote", "  }",
            "}"]
        chains = _fetched(plugin).composition.chain
        assert [d.label for d in chains[0].dataset] == ["peak", "remote"]

    def test_repeated_definition_is_repeated_lines(self):
        plugin = MockMultiPlugin()
        model = plugin.container.controlParameters.model
        model.append("a.pdb")
        model.append("b.pdb")
        assert plugin.extract_phil_lines() == ["model = a.pdb", "model = b.pdb"]
        assert _fetched(plugin).model == ["a.pdb", "b.pdb"]

    def test_flat_extraction_keeps_its_contract(self):
        plugin = MockMultiPlugin()
        cp = plugin.container.controlParameters
        cp.composition.composition__solvent.value = 0.5
        chain = cp.composition.composition__chain
        chain.append(chain.makeItem())
        cp.model.append("a.pdb")
        # scalars and repeated leaves as pairs; blocks are for the lines API
        assert plugin.extract_phil_parameters() == [
            ("composition.solvent", "0.5"), ("model", "a.pdb")]

    def test_items_round_trip_through_params_xml(self):
        plugin = MockMultiPlugin()
        chain = plugin.container.controlParameters.composition.composition__chain
        chain.append(chain.makeItem())
        chain.append(chain.makeItem())
        chain[1].composition__chain__nres.value = 77
        with tempfile.TemporaryDirectory() as tmpdir:
            xml_path = os.path.join(tmpdir, "input_params.xml")
            plugin.saveDataToXml(xml_path)
            again = MockMultiPlugin()
            again.loadDataFromXml(xml_path)
        chain2 = again.container.controlParameters.composition.composition__chain
        assert len(chain2) == 2
        assert chain2[1].composition__chain__nres.value == 77
        assert chain2[0].composition__chain__nres.getValueState() == ValueState.NOT_SET

    def test_json_offers_a_container_shaped_template(self):
        import json
        from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder
        plugin = MockMultiPlugin()
        chain = plugin.container.controlParameters.composition.composition__chain
        encoded = json.loads(json.dumps(chain, cls=CCP4i2JsonEncoder))
        assert encoded["_baseClass"] == "CList"
        template = encoded["_subItem"]
        assert template["_baseClass"] == "CContainer"
        assert "composition__chain__num" in template["_value"]
        assert template["_value"]["composition__chain__num"]["_value"] == 1
        assert "[?]" in template["_objectPath"]


# ---------------------------------------------------------------------------
# Shims: block entries, and their targets leave the generic tree
# ---------------------------------------------------------------------------

from ccp4i2.utils.phil_shims import PhilShim, MtzFileShim, PdbFileShim


class ChainShim(PhilShim):
    """Pretends to map ASU contents to composition.chain blocks."""

    def __init__(self):
        self.phil_chain_path = "composition.chain"

    def convert(self, container, work_directory):
        return [
            ("composition.chain", [("nres", 120), ("num", 2)]),
            ("composition.chain", [("nres", 50)]),
            ("composition.solvent", 0.5),
        ]


class MockShimmedPlugin(MockMultiPlugin):
    TASKNAME = "mock_shimmed_plugin"

    def get_shim_definitions(self):
        return [ChainShim()]


class TestShimEntries:

    def test_phil_targets_read_from_the_phil_attributes(self):
        assert MtzFileShim("HKLIN", "picard.hklin").phil_targets() == ["picard.hklin"]
        assert PdbFileShim("XYZIN", ["a.model", "b.model"]).phil_targets() == ["a.model", "b.model"]
        assert ChainShim().phil_targets() == ["composition.chain"]

    def test_a_shim_target_is_not_offered_in_the_tree(self):
        assert hasattr(MockMultiPlugin().container.controlParameters.composition,
                       "composition__chain")
        assert not hasattr(MockShimmedPlugin().container.controlParameters.composition,
                           "composition__chain")

    def test_shim_blocks_reach_the_tool(self):
        plugin = MockShimmedPlugin()
        fetched = _fetched(plugin)
        assert [(c.nres, c.num) for c in fetched.composition.chain] == [(120, 2), (50, 1)]
        assert fetched.composition.solvent == 0.5
