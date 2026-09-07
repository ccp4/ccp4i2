"""What the RNP task and its pipeline accept as something to refine.

The task needs either a solution file or search models placed at the
origin of their coordinates; the pipeline builds the latter from atom
selections of a parent model, and checks the selections first.
"""
from pathlib import Path

import pytest

pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)")

import ccp4i2
from ccp4i2.core.tasks import get_plugin_class

GAMMA = Path(ccp4i2.__file__).parent / "demo_data" / "gamma" / "gamma_model.pdb"


def codes(plugin):
    return [r["code"] for r in plugin.validity()._reports]


def make(task, tmp_path):
    return get_plugin_class(task)(workDirectory=str(tmp_path), name=task)


class TestTask:

    def test_nothing_placed_is_an_error(self, tmp_path):
        assert 116 in codes(make("phaser_mr_rnp_phil", tmp_path))

    def test_fixed_ensembles_stand_in_for_a_solution_file(self, tmp_path):
        plugin = make("phaser_mr_rnp_phil", tmp_path)
        plugin.container.inputData.FIXENSEMBLES.append("Fragment_1")
        assert 116 not in codes(plugin)


class TestPipelineSelections:

    def plugin(self, tmp_path, *selections):
        plugin = make("phaser_rnp_pipeline_phil", tmp_path)
        inp = plugin.container.inputData
        inp.XYZIN_PARENT.setFullPath(str(GAMMA))
        for text in selections:
            inp.SELECTIONS.append(inp.SELECTIONS.makeItem())
            inp.SELECTIONS[-1].text.set(text)
        return plugin

    def test_no_selection_means_the_whole_model(self, tmp_path):
        assert not {210, 212} & set(codes(self.plugin(tmp_path)))

    def test_a_selection_of_nothing_is_an_error(self, tmp_path):
        assert 210 in codes(self.plugin(tmp_path, "A/", "Z/"))

    def test_overlapping_selections_are_an_error(self, tmp_path):
        assert 212 in codes(self.plugin(tmp_path, "A/", "A/"))

    def test_disjoint_selections_pass(self, tmp_path):
        assert not {210, 212} & set(codes(self.plugin(tmp_path, "A/703-750", "A/751-822")))

    def test_fragments_are_cut_and_placed(self, tmp_path):
        plugin = self.plugin(tmp_path, "A/703-750", "A/751-822")
        sub = plugin.makePluginObject("phaser_mr_rnp_phil")
        assert plugin.createEnsembles(sub)
        inp = sub.container.inputData
        assert [str(e.label) for e in inp.ENSEMBLES] == ["Fragment_1", "Fragment_2"]
        assert [str(x) for x in inp.FIXENSEMBLES] == ["Fragment_1", "Fragment_2"]
        sizes = [sum(1 for line in (tmp_path / f"Fragment_{i}.pdb").read_text().splitlines()
                     if line.startswith("ATOM")) for i in (1, 2)]
        assert all(sizes) and sum(sizes) == sum(1 for line in GAMMA.read_text().splitlines()
                                                if line.startswith("ATOM"))
        assert 116 not in codes(sub)
