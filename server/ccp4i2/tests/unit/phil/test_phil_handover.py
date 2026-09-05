"""Handing PHIL parameters between tasks that host the same master.

A pipeline hosts the tool's PHIL and gives the tree to the sub-job that
runs the tool; the two containers have the same shape, so the copy walks
them in parallel. compose_master_phil() lets the pipeline keep its own
parameters beside the tool's in one master.
"""
import pytest

parse = pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)").parse

from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.core.base_object.fundamental_types import CList

TOOL = parse("""
    tool {
      resolution = 2.0
        .type = float
      title = None
        .type = str
      chain
        .optional = True
        .multiple = True
      {
        nres = None
          .type = int
        num = 1
          .type = int
      }
      model = None
        .type = path
        .multiple = True
    }
""")

PIPELINE_EXTRA = """
pipeline {
  run_refmac = True
    .type = bool
  cycles = 10
    .type = int
}
"""


class Tool(PhilPluginScript):
    TASKNAME = "tool"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return TOOL

    def get_phil_exclude_scopes(self):
        return []

    def get_command_target(self):
        return "echo"


class Pipeline(Tool):
    TASKNAME = "pipeline"

    def get_master_phil(self):
        return self.compose_master_phil(TOOL, PIPELINE_EXTRA)

    def get_phil_exclude_scopes(self):
        return []


class TestFindAndSet:

    def test_find_and_set_a_scalar(self):
        t = Tool()
        t.set_phil("tool.resolution", 1.5)
        assert t.find_phil("tool.resolution").get() == 1.5
        assert t.extract_phil_parameters() == [("tool.resolution", "1.5")]

    def test_set_a_repeated_definition(self):
        t = Tool()
        t.set_phil("tool.model", ["a.pdb", "b.pdb"])
        assert [str(m) for m in t.container.controlParameters.tool.tool__model] == ["a.pdb", "b.pdb"]

    def test_unknown_path_is_an_error(self):
        with pytest.raises(AttributeError):
            Tool().set_phil("tool.nothing", 1)


class TestHandover:

    def test_set_values_lists_and_expert_level_cross_over(self):
        a, b = Tool(), Tool()
        a.set_phil("tool.title", "x")
        chain = a.container.controlParameters.tool.tool__chain
        chain.append(chain.makeItem()); chain[0].tool__chain__nres.value = 120
        chain.append(chain.makeItem()); chain[1].tool__chain__num.value = 3
        a.set_phil("tool.model", ["a.pdb"])
        a.container.controlParameters.PHIL_EXPERT_LEVEL.set(2)
        a.hand_phil_to(b)
        assert b.extract_phil_lines() == a.extract_phil_lines()
        assert b.extract_phil_lines() == [
            "tool.title = x", "tool.chain {", "  nres = 120", "}", "tool.chain {", "  num = 3", "}",
            "tool.model = a.pdb"]
        assert b.container.controlParameters.PHIL_EXPERT_LEVEL.get() == 2

    def test_defaults_are_not_copied_as_set(self):
        a, b = Tool(), Tool()
        a.hand_phil_to(b)
        assert b.extract_phil_parameters() == []


class TestComposedMaster:

    def test_the_pipeline_has_both_scopes(self):
        p = Pipeline()
        cp = p.container.controlParameters
        assert hasattr(cp, "tool") and hasattr(cp, "pipeline")
        assert cp.pipeline.pipeline__run_refmac.value is True
        assert cp.pipeline.pipeline__cycles.value == 10

    def test_the_base_master_is_not_mutated(self):
        Pipeline()
        assert [o.name for o in TOOL.active_objects()] == ["tool"]

    def test_handover_to_the_tool_ignores_the_pipeline_scope(self):
        p, t = Pipeline(), Tool()
        p.set_phil("tool.resolution", 1.7)
        p.set_phil("pipeline.cycles", 5)
        p.hand_phil_to(t)
        assert t.extract_phil_parameters() == [("tool.resolution", "1.7")]
        # and the tool's working phil is valid against its own master
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            t.workDirectory = tmp
            path = t.build_working_phil()
            assert TOOL.fetch(sources=[parse(file_name=path)]).extract().tool.resolution == 1.7
