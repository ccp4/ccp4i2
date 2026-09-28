"""One authority decides how a parameter is named on an i2run command line.

The renderer behind the job panel's "i2run command" button and i2run's own
argparse builder must agree, or the button prints commands the parser rejects.
They used to compute the answer separately and agreed only by coincidence --
both happened to minimise against every leaf in the container, including the
output sections. Nothing recorded that, so these tests do.
"""
import pytest

from ccp4i2.lib.utils.parameters.argument_names import (
    argument_names,
    compute_minimum_paths,
    leaf_paths,
)


class TestMinimumPaths:
    """The rule itself, on hand-built tables -- no CCP4, no plugins."""

    def test_a_unique_name_shortens_to_itself(self):
        kws = [{"path": "container.inputData.F_SIGF"}, {"path": "container.controlParameters.FRAC"}]
        compute_minimum_paths(kws)
        assert kws[0]["minimumPath"] == "F_SIGF"
        assert kws[1]["minimumPath"] == "FRAC"

    def test_a_shared_name_grows_until_it_is_unique(self):
        """servalcat_pipe's real shape: two XYZINs, so neither may be bare."""
        kws = [
            {"path": "container.inputData.XYZIN"},
            {"path": "container.metalCoordWrapper.inputData.XYZIN"},
        ]
        compute_minimum_paths(kws)
        assert kws[0]["minimumPath"] == "container.inputData.XYZIN"
        assert kws[1]["minimumPath"] == "metalCoordWrapper.inputData.XYZIN"
        assert kws[0]["isAmbiguousSimpleName"] is True

    def test_suffixes_are_compared_per_element_not_per_character(self):
        """A string endswith would match halfway into an element: 'B.X' is not
        a suffix of 'a.AB.X' by elements, though 'B.X' ends the string."""
        kws = [{"path": "root.AB.X"}, {"path": "root.B.X"}]
        compute_minimum_paths(kws)
        assert kws[0]["minimumPath"] == "AB.X"
        assert kws[1]["minimumPath"] == "B.X"

    def test_the_shortest_bearer_of_an_ambiguous_name_is_flagged(self):
        kws = [
            {"path": "container.inputData.XYZIN"},
            {"path": "container.deep.nested.inputData.XYZIN"},
        ]
        compute_minimum_paths(kws)
        assert kws[0]["isShortestForSimpleName"] is True
        assert kws[1]["isShortestForSimpleName"] is False


class TestAgainstRealTasks:
    """Needs the plugins, so it needs CCP4 -- skips rather than fails."""

    @pytest.fixture(autouse=True)
    def _needs_plugins(self):
        pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

    # Spread of shapes: one flat task, one with sub-wrappers and colliding
    # names, one pipeline, one PHIL task.
    TASKS = ("freerflag", "servalcat_pipe", "prosmart_refmac", "aimless_pipe")

    @pytest.mark.parametrize("task_name", TASKS)
    def test_the_renderer_and_the_parser_use_the_same_table(self, task_name):
        """The invariant that was previously only a coincidence: every name the
        renderer would emit is a name i2run's keyword table carries."""
        from ccp4i2.cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase
        from ccp4i2.core.tasks import get_plugin_class

        plugin = get_plugin_class(task_name)(parent=None)
        rendered = argument_names(plugin.container)

        parser_table = {
            kw["path"]: kw["minimumPath"]
            for kw in CCP4i2RunnerBase.keywordsOfTaskName(task_name)
            if "path" in kw and "minimumPath" in kw
        }

        assert parser_table, f"no keywords for {task_name}"
        for path, name in parser_table.items():
            assert rendered.get(path) == name, (
                f"{task_name}: parser calls {path} '{name}', "
                f"renderer calls it '{rendered.get(path)}'"
            )

    @pytest.mark.parametrize("task_name", TASKS)
    def test_every_parameter_gets_a_name(self, task_name):
        from ccp4i2.core.tasks import get_plugin_class

        plugin = get_plugin_class(task_name)(parent=None)
        leaves = leaf_paths(plugin.container)
        names = argument_names(plugin.container)

        assert len(names) == len(leaves)
        assert all(names.values()), "a parameter rendered with an empty name"

    def test_a_file_is_one_parameter_not_one_per_field(self):
        """Recursion stops at non-containers, so CDataFile's project/baseName/
        dbFileId are not separate arguments -- the command names the file."""
        from ccp4i2.core.tasks import get_plugin_class

        plugin = get_plugin_class("freerflag")(parent=None)
        names = argument_names(plugin.container)

        assert "container.inputData.F_SIGF" in names
        assert not [p for p in names if p.endswith(".F_SIGF.dbFileId")]
