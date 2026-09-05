"""The shims from CCP4i2's typed inputs to Phaser's PHIL entries."""
import os
import tempfile

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ModelData import (
    CAsuDataFile, CEnsembleList, CSeqDataFile, CSequenceWithCopiesList)
from ccp4i2.core.CCP4XtalData import CObsDataFile
from ccp4i2.core.base_object.fundamental_types import CString, CFloat
from ccp4i2.utils.phil_shims import FixedPhilShim
from ccp4i2.wrappers.phaser_phil.script.phaser_shims import (
    CompositionShim, EnsembleListShim, ObsDataShim)


def demo(*parts):
    from ccp4i2.core.CCP4Utils import getCCP4I2Dir
    return os.path.join(getCCP4I2Dir(), "demo_data", *parts)


@pytest.fixture(name="plugin")
def plugin_fixture():
    script = CPluginScript(name="shim_test")
    inp = script.container.inputData
    inp.addContent(CEnsembleList, "ENSEMBLES")
    inp.addContent(CString, "COMP_BY")
    inp.addContent(CAsuDataFile, "ASUFILE")
    inp.addContent(CSequenceWithCopiesList, "SEQUENCES")
    inp.addContent(CFloat, "SOLVENT_FRACTION")
    return script


class TestEnsembleListShim:

    def test_ensemble_and_search_blocks(self, plugin):
        ensembles = plugin.container.inputData.ENSEMBLES
        e = ensembles.makeItem()
        ensembles.append(e)
        e.label.set("model")
        e.number.set(2)
        e.use.set(True)
        item = e.pdbItemList.makeItem()
        e.pdbItemList.append(item)
        item.structure.setFullPath("/x/gamma_model.pdb")
        item.identity_to_target.set(0.9)
        entries = EnsembleListShim("ENSEMBLES", "phaser.ensemble", "phaser.search").convert(
            plugin.container, "/tmp")
        assert entries[0] == ("phaser.ensemble", [
            ("model_id", "model"),
            ("coordinates", [("pdb", "/x/gamma_model.pdb"), ("identity", 0.9)])])
        assert entries[1] == ("phaser.search", [("ensembles", "model"), ("copies", 2)])

    def test_rmsd_and_variance_from_remarks(self, plugin):
        ensembles = plugin.container.inputData.ENSEMBLES
        e = ensembles.makeItem(); ensembles.append(e); e.label.set("m"); e.use.set(False)
        a = e.pdbItemList.makeItem(); e.pdbItemList.append(a)
        a.structure.setFullPath("/x/a.pdb"); a.rms_to_target.set(1.2)
        b = e.pdbItemList.makeItem(); e.pdbItemList.append(b)
        b.structure.setFullPath("/x/b.cif")
        shim = EnsembleListShim("ENSEMBLES", "phaser.ensemble", "phaser.search",
                                path_map={"/x/b.cif": "/job/b.pdb"})
        entries = shim.convert(plugin.container, "/tmp")
        coords = [c for k, c in entries[0][1] if k == "coordinates"]
        assert coords[0] == [("pdb", "/x/a.pdb"), ("rmsd", 1.2)]
        assert coords[1] == [("pdb", "/job/b.pdb"), ("read_variance_from_pdb_remarks", True)]
        assert len(entries) == 1        # 'use' off: defined, not searched for

    def test_targets(self):
        assert EnsembleListShim("E", "phaser.ensemble", "phaser.search").phil_targets() == [
            "phaser.ensemble", "phaser.search"]


class TestCompositionShim:

    def test_default_is_average_solvent(self, plugin):
        plugin.container.inputData.COMP_BY.set("DEFAULT")
        assert CompositionShim().convert(plugin.container, "/tmp") == [("phaser.composition.solvent", 0.5)]

    def test_solvent_fraction(self, plugin):
        inp = plugin.container.inputData
        inp.COMP_BY.set("SOLVENT"); inp.SOLVENT_FRACTION.set(0.42)
        assert CompositionShim().convert(plugin.container, "/tmp") == [("phaser.composition.solvent", 0.42)]

    def test_sequence_files_with_copies(self, plugin):
        inp = plugin.container.inputData
        inp.COMP_BY.set("SEQUENCES")
        item = inp.SEQUENCES.makeItem(); inp.SEQUENCES.append(item)
        item.seqFile.setFullPath("/x/a.seq"); item.nCopies.set(3)
        item = inp.SEQUENCES.makeItem(); inp.SEQUENCES.append(item)
        item.seqFile.setFullPath("/x/b.seq")
        assert CompositionShim().convert(plugin.container, "/tmp") == [
            ("phaser.composition.chain", [("sequence_file", "/x/a.seq"), ("num", 3)]),
            ("phaser.composition.chain", [("sequence_file", "/x/b.seq"), ("num", 1)])]

    def test_asu_file_writes_one_fasta_per_sequence(self, plugin):
        inp = plugin.container.inputData
        inp.COMP_BY.set("ASU")
        inp.ASUFILE.setFullPath(demo("gamma", "gamma.asu.xml"))
        with tempfile.TemporaryDirectory() as tmp:
            entries = CompositionShim().convert(plugin.container, tmp)
            assert len(entries) >= 1
            path, fields = entries[0]
            assert path == "phaser.composition.chain"
            fields = dict(fields)
            assert fields["sequence_file"].startswith(tmp) and fields["num"] >= 1
            text = open(fields["sequence_file"]).read()
            assert text.startswith(">") and len(text.splitlines()) >= 2

    def test_targets(self):
        assert sorted(CompositionShim().phil_targets()) == [
            "phaser.composition.chain", "phaser.composition.solvent"]


class TestObsDataShim:

    def test_names_intensities_or_amplitudes(self):
        pytest.importorskip("gemmi", reason="makeHklin needs gemmi")
        import shutil
        if shutil.which("servalcat") is None:
            pytest.skip("the IPAIR -> IMEAN conversion runs servalcat (CCP4)")
        with tempfile.TemporaryDirectory() as tmp:
            plugin = CPluginScript(name="obs_shim_test", workDirectory=tmp)
            inp = plugin.container.inputData
            inp.addContent(CObsDataFile, "F_SIGF")
            inp.F_SIGF.setFullPath(demo("gamma", "merged_intensities_Xe.mtz"))
            inp.F_SIGF.contentFlag.set(CObsDataFile.CONTENT_FLAG_IPAIR)
            shim = ObsDataShim(plugin, "F_SIGF", "phaser.hklin", "phaser.labin")
            assert shim.convert(plugin.container, tmp) == []      # not prepared yet
            error = shim.prepare()
            assert error.maxSeverity() <= 2, str(error)
            entries = dict(shim.convert(plugin.container, tmp))
            assert entries["phaser.labin"] == "I,SIGI"
            assert os.path.exists(entries["phaser.hklin"])
            assert shim.phil_targets() == ["phaser.hklin", "phaser.labin"]


class TestFixedPhilShim:

    def test_values_and_targets(self):
        shim = FixedPhilShim({"phaser.keywords.general.root": "{work_directory}/PHASER", "a.b": 3})
        assert shim.convert(None, "/job") == [("phaser.keywords.general.root", "/job/PHASER"), ("a.b", 3)]
        assert shim.phil_targets() == ["phaser.keywords.general.root", "a.b"]


from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.base_object.fundamental_types import CInt, CBoolean, CList
from ccp4i2.wrappers.phaser_phil.script.phaser_shims import (
    EpCrystalShim, LlgCompletionShim, PartialModelShim)


@pytest.fixture(name="ep_plugin")
def ep_plugin_fixture():
    script = CPluginScript(name="ep_shim_test")
    inp = script.container.inputData
    inp.addContent(CPdbDataFile, "XYZIN_HA")
    inp.addContent(CPdbDataFile, "XYZIN_PARTIAL")
    inp.addContent(CString, "PARTIAL_BY")
    inp.addContent(CList, "ELEMENTS")
    inp.addContent(CInt, "LLGC_CYCLES")
    inp.addContent(CBoolean, "PURE_ANOMALOUS")
    inp.addContent(CFloat, "WAVELENGTH")
    return script


class TestLlgCompletionShim:

    def test_elements_cycles_and_method(self, ep_plugin):
        inp = ep_plugin.container.inputData
        inp.ELEMENTS.append("Xe"); inp.ELEMENTS.append("Se")
        inp.LLGC_CYCLES.set(20); inp.PURE_ANOMALOUS.set(True)
        assert LlgCompletionShim().convert(ep_plugin.container, "/tmp") == [
            ("phaser.keywords.llgcompletion", [("complete", True), ("scatterer", "Xe Se"),
                                               ("ncycle", 20), ("method", "imaginary")])]

    def test_nothing_without_elements(self, ep_plugin):
        assert LlgCompletionShim().convert(ep_plugin.container, "/tmp") == []


class TestPartialModelShim:

    def test_model(self, ep_plugin):
        inp = ep_plugin.container.inputData
        inp.PARTIAL_BY.set("MODEL"); inp.XYZIN_PARTIAL.setFullPath("/x/partial.pdb")
        assert PartialModelShim().convert(ep_plugin.container, "/tmp") == [
            ("phaser.keywords.partial", [("mode", "model"), ("pdb", "/x/partial.pdb")])]

    def test_none(self, ep_plugin):
        ep_plugin.container.inputData.PARTIAL_BY.set("NONE")
        assert PartialModelShim().convert(ep_plugin.container, "/tmp") == []


class TestEpCrystalShim:

    def test_crystal_block(self):
        pytest.importorskip("iotbx", reason="needs iotbx (CCP4/cctbx)")
        import shutil
        if shutil.which("servalcat") is None:
            pytest.skip("the IPAIR -> FPAIR conversion runs servalcat (CCP4)")
        with tempfile.TemporaryDirectory() as tmp:
            plugin = CPluginScript(name="ep_crystal_test", workDirectory=tmp)
            inp = plugin.container.inputData
            inp.addContent(CObsDataFile, "F_SIGF"); inp.addContent(CPdbDataFile, "XYZIN_HA")
            inp.addContent(CFloat, "WAVELENGTH")
            inp.F_SIGF.setFullPath(demo("gamma", "merged_intensities_Xe.mtz"))
            inp.F_SIGF.contentFlag.set(CObsDataFile.CONTENT_FLAG_IPAIR)
            inp.XYZIN_HA.setFullPath(demo("gamma", "heavy_atoms.pdb")); inp.WAVELENGTH.set(1.542)
            shim = EpCrystalShim(plugin)
            assert shim.prepare().maxSeverity() <= 2
            (path, fields), = shim.convert(plugin.container, tmp)
            assert path == "phaser.crystal"
            fields = dict(fields)
            assert fields["pdb_file"].endswith("heavy_atoms.pdb")
            dataset = dict(fields["dataset"])
            assert dataset["labin"] == "Fplus,SIGFplus,Fminus,SIGFminus,merged"
            assert dataset["wavelength"] == 1.542 and os.path.exists(dataset["hklin"])
            assert shim.phil_targets() == ["phaser.crystal"]


class TestMrValidity:
    """The task warns before Phaser refuses."""

    def test_two_molecules_in_one_ensemble_is_advised_against(self):
        pytest.importorskip("lxml")
        from ccp4i2.core.tasks import get_plugin_class
        with tempfile.TemporaryDirectory() as tmp:
            plugin = get_plugin_class("phaser_mr_auto_phil")(workDirectory=tmp, name="t")
            inp = plugin.container.inputData
            e = inp.ENSEMBLES.makeItem(); inp.ENSEMBLES.append(e); e.label.set("both"); e.number.set(1); e.use.set(True)
            for name in ("beta.pdb", "blip.pdb"):
                item = e.pdbItemList.makeItem(); e.pdbItemList.append(item)
                item.structure.setFullPath(f"/x/{name}"); item.identity_to_target.set(0.9)
            codes = [r["code"] for r in plugin.validity()._reports]
            assert 115 in codes
            assert 112 not in codes
