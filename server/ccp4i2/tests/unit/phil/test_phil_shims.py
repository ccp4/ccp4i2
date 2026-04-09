"""
Tests for PHIL shims — converting CCP4i2 rich file types to PHIL-compatible values.
"""

import os
import pytest
import tempfile
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CMtzDataFile
from ccp4i2.core.CCP4ModelData import CPdbDataFile, CAsuDataFile, CDictDataFile
from ccp4i2.utils.phil_shims import (
    MtzFileShim,
    PdbFileShim,
    AsuContentShim,
    DictFileShim,
)

# Try to find demo data path
try:
    from ccp4i2.core.CCP4Utils import getCCP4I2Dir
except ImportError:
    def getCCP4I2Dir():
        return str(Path(__file__).parent.parent)


def demoData(*paths):
    return os.path.join(getCCP4I2Dir(), "demo_data", *paths)


# ---------------------------------------------------------------------------
# Helper to create a plugin with file objects in inputData
# ---------------------------------------------------------------------------

def _make_plugin_with_file(file_class, file_name, file_path=None):
    """Create a minimal plugin and attach a file object to inputData."""
    script = CPluginScript(name="shim_test")
    file_obj = file_class(parent=script.container.inputData, name=file_name)
    if file_path:
        file_obj.setFullPath(file_path)
    return script


# ---------------------------------------------------------------------------
# MtzFileShim
# ---------------------------------------------------------------------------

class TestMtzFileShim:

    def test_converts_set_mtz_file(self):
        mtz_path = demoData("gamma", "merged_intensities_Xe.mtz")
        if not os.path.exists(mtz_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CMtzDataFile, "HKLIN", mtz_path)
        shim = MtzFileShim("HKLIN", "picard.hklin")
        result = shim.convert(script.container, "/tmp")

        assert len(result) == 1
        assert result[0][0] == "picard.hklin"
        assert result[0][1] == mtz_path

    def test_returns_empty_when_not_set(self):
        script = _make_plugin_with_file(CMtzDataFile, "HKLIN")
        shim = MtzFileShim("HKLIN", "picard.hklin")
        result = shim.convert(script.container, "/tmp")
        assert result == []

    def test_returns_empty_when_attribute_missing(self):
        script = CPluginScript(name="shim_test")
        shim = MtzFileShim("NONEXISTENT", "picard.hklin")
        result = shim.convert(script.container, "/tmp")
        assert result == []


# ---------------------------------------------------------------------------
# PdbFileShim
# ---------------------------------------------------------------------------

class TestPdbFileShim:

    def test_converts_set_pdb_file(self):
        pdb_path = demoData("gamma", "gamma_model.pdb")
        if not os.path.exists(pdb_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CPdbDataFile, "XYZIN", pdb_path)
        shim = PdbFileShim("XYZIN", "picard.xyzin")
        result = shim.convert(script.container, "/tmp")

        assert len(result) == 1
        assert result[0][0] == "picard.xyzin"
        assert result[0][1] == pdb_path

    def test_multi_path_output(self):
        pdb_path = demoData("gamma", "gamma_model.pdb")
        if not os.path.exists(pdb_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CPdbDataFile, "XYZIN", pdb_path)
        shim = PdbFileShim("XYZIN", ["picard.xyzin", "phasertng.model.filename"])
        result = shim.convert(script.container, "/tmp")

        assert len(result) == 2
        assert result[0] == ("picard.xyzin", pdb_path)
        assert result[1] == ("phasertng.model.filename", pdb_path)

    def test_returns_empty_when_not_set(self):
        script = _make_plugin_with_file(CPdbDataFile, "XYZIN")
        shim = PdbFileShim("XYZIN", "picard.xyzin")
        result = shim.convert(script.container, "/tmp")
        assert result == []

    def test_returns_empty_when_attribute_missing(self):
        script = CPluginScript(name="shim_test")
        shim = PdbFileShim("NONEXISTENT", "picard.xyzin")
        result = shim.convert(script.container, "/tmp")
        assert result == []


# ---------------------------------------------------------------------------
# AsuContentShim
# ---------------------------------------------------------------------------

class TestAsuContentShim:

    def test_converts_asu_to_fasta(self):
        asu_path = demoData("gamma", "gamma.asu.xml")
        if not os.path.exists(asu_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CAsuDataFile, "ASUIN", asu_path)
        shim = AsuContentShim("ASUIN", "picard.seqin")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = shim.convert(script.container, tmpdir)

            assert len(result) == 1
            assert result[0][0] == "picard.seqin"

            fasta_path = result[0][1]
            assert os.path.exists(fasta_path)

            # Verify FASTA format
            with open(fasta_path) as f:
                content = f.read()

            assert content.startswith(">")
            lines = content.strip().split("\n")
            # First line is header
            assert lines[0].startswith(">")
            # Sequence lines should be <= 80 chars
            for line in lines[1:]:
                assert len(line) <= 80

    def test_fasta_contains_sequence(self):
        asu_path = demoData("gamma", "gamma.asu.xml")
        if not os.path.exists(asu_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CAsuDataFile, "ASUIN", asu_path)
        shim = AsuContentShim("ASUIN", "picard.seqin")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = shim.convert(script.container, tmpdir)
            with open(result[0][1]) as f:
                content = f.read()

            # gamma.asu.xml contains the Gamma sequence starting with MHHHHHH
            assert "MHHHHHH" in content

    def test_returns_empty_when_not_set(self):
        script = _make_plugin_with_file(CAsuDataFile, "ASUIN")
        shim = AsuContentShim("ASUIN", "picard.seqin")
        result = shim.convert(script.container, "/tmp")
        assert result == []


# ---------------------------------------------------------------------------
# DictFileShim
# ---------------------------------------------------------------------------

class TestDictFileShim:

    def test_converts_set_dict_file(self):
        dict_path = demoData("mdm2", "4hg7.cif")
        if not os.path.exists(dict_path):
            pytest.skip("Demo data not available")

        script = _make_plugin_with_file(CDictDataFile, "DICT", dict_path)
        shim = DictFileShim("DICT", "phasertng.cluster_compound.filename")
        result = shim.convert(script.container, "/tmp")

        assert len(result) == 1
        assert result[0][0] == "phasertng.cluster_compound.filename"
        assert result[0][1] == dict_path

    def test_returns_empty_when_not_set(self):
        script = _make_plugin_with_file(CDictDataFile, "DICT")
        shim = DictFileShim("DICT", "phasertng.cluster_compound.filename")
        result = shim.convert(script.container, "/tmp")
        assert result == []

    def test_returns_empty_when_attribute_missing(self):
        script = CPluginScript(name="shim_test")
        shim = DictFileShim("NONEXISTENT", "some.phil.path")
        result = shim.convert(script.container, "/tmp")
        assert result == []
