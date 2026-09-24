"""
Tests for CUnmergedDataContent.loadFile() implementation.

Tests loading of unmerged reflection data from multiple formats:
- MTZ (merged and unmerged)
- Scalepack (.sca - merged and unmerged formats)
- XDS (INTEGRATE.HKL, XDS_ASCII.HKL)
"""

import pytest
from pathlib import Path
from ccp4i2.core.CCP4XtalData import CUnmergedDataContent, CUnmergedDataFile
from ccp4i2.core.base_object.error_reporting import CErrorReport
from ccp4i2.core import CCP4Utils

# Get the path to test data
CCP4I2_ROOT = Path(CCP4Utils.getCCP4I2Dir())
TEST_MTZ = CCP4I2_ROOT / "wrappers/pointless/test_data/brap_pk_6A.mtz"
TEST_SCA_MERGED = CCP4I2_ROOT / "demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled.sca"
TEST_SCA_UNMERGED = CCP4I2_ROOT / "demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled_unmerged.sca"
# XDS writes .HKL, and so does scalepack in some pipelines, so these exercise
# the content sniff rather than the extension.
TEST_XDS_ASCII = CCP4I2_ROOT / "demo_data/ceue/apo-ceue-sad-sweep1.hkl"
TEST_XDS_INTEGRATE = (
    CCP4I2_ROOT
    / "demo_data/baz2b/BAZ2BA_x828.xia2/3daii-run/DataFiles/Integrate"
    / "nt5073v16_xBAZ2BAx8281_SAD_SWEEP1_INTEGRATE.HKL"
)
TEST_XDS_CORRECT = (
    CCP4I2_ROOT
    / "demo_data/baz2b/BAZ2BA_x828.xia2/3daii-run/DataFiles/Integrate"
    / "nt5073v16_xBAZ2BAx8281_SAD_SWEEP1_CORRECT.HKL"
)


class TestCUnmergedDataContent:
    """Test CUnmergedDataContent.loadFile() method for various formats."""

    def test_loadfile_nonexistent_file(self):
        """Test loading a non-existent file returns error."""
        data = CUnmergedDataContent()
        error = data.loadFile('/nonexistent/file.mtz')

        assert isinstance(error, CErrorReport)
        assert error.count() > 0
        assert error.maxSeverity() > 0

    def test_loadfile_empty_path(self):
        """Test loading with empty path returns empty error report."""
        data = CUnmergedDataContent()
        error = data.loadFile('')

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    def test_loadfile_none_path(self):
        """Test loading with None path returns empty error report."""
        data = CUnmergedDataContent()
        error = data.loadFile(None)

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    @pytest.mark.skipif(
        not TEST_MTZ.exists(),
        reason="Test MTZ file not available"
    )
    def test_loadfile_unmerged_mtz(self):
        """Test loading an unmerged MTZ file."""
        test_mtz = str(TEST_MTZ)

        data = CUnmergedDataContent()
        error = data.loadFile(test_mtz)

        # Should load without errors
        if error.count() > 0:
            print(f"\nError loading MTZ: {error}")
        assert error.count() == 0

        # Check format detection
        assert data.format.value == 'mtz'
        assert data.merged.value == 'unmerged'

        # Check metadata extraction
        assert data.cell is not None
        assert data.cell.a.value > 0
        assert data.spaceGroup.value is not None
        assert len(data.spaceGroup.value) > 0
        assert data.knowncell.value == True

        # Print summary
        print(f"\nLoaded unmerged MTZ:")
        print(f"  Format: {data.format.value}")
        print(f"  Merged: {data.merged.value}")
        print(f"  Space group: {data.spaceGroup.value}")
        print(f"  Cell: a={data.cell.a.value:.2f}")
        print(f"  Known cell: {data.knowncell.value}")

    @pytest.mark.skipif(
        not TEST_SCA_MERGED.exists(),
        reason="Test Scalepack merged file not available"
    )
    def test_loadfile_scalepack_merged(self):
        """Test loading a merged Scalepack file."""
        test_sca = str(TEST_SCA_MERGED)

        data = CUnmergedDataContent()
        error = data.loadFile(test_sca)

        # Should load without errors
        if error.count() > 0:
            print(f"\nError loading Scalepack: {error}")
        assert error.count() == 0

        # Check format detection
        assert data.format.value == 'sca'
        assert data.merged.value == 'merged'

        # Check metadata extraction
        assert data.cell is not None
        assert data.cell.a.value > 0
        assert data.spaceGroup.value is not None
        assert len(data.spaceGroup.value) > 0
        assert data.knowncell.value == True
        assert data.knownwavelength.value == False  # Scalepack doesn't have wavelength

        # Print summary
        print(f"\nLoaded merged Scalepack:")
        print(f"  Format: {data.format.value}")
        print(f"  Merged: {data.merged.value}")
        print(f"  Space group: {data.spaceGroup.value}")
        print(f"  Cell: a={data.cell.a.value:.2f}, b={data.cell.b.value:.2f}, c={data.cell.c.value:.2f}")
        print(f"  Known cell: {data.knowncell.value}")
        print(f"  Known wavelength: {data.knownwavelength.value}")

    @pytest.mark.skipif(
        not TEST_SCA_UNMERGED.exists(),
        reason="Test Scalepack unmerged file not available"
    )
    def test_loadfile_scalepack_unmerged(self):
        """Test loading an unmerged Scalepack file."""
        test_sca = str(TEST_SCA_UNMERGED)

        data = CUnmergedDataContent()
        error = data.loadFile(test_sca)

        # Should load without errors
        if error.count() > 0:
            print(f"\nError loading Scalepack unmerged: {error}")
        assert error.count() == 0

        # Check format detection
        assert data.format.value == 'sca'
        assert data.merged.value == 'unmerged'

        # Check that cell is NOT known (unmerged Scalepack doesn't have cell)
        assert data.knowncell.value == False
        assert data.knownwavelength.value == False

        # Space group should be detected
        assert data.spaceGroup.value is not None
        assert len(data.spaceGroup.value) > 0

        # Print summary
        print(f"\nLoaded unmerged Scalepack:")
        print(f"  Format: {data.format.value}")
        print(f"  Merged: {data.merged.value}")
        print(f"  Space group: {data.spaceGroup.value}")
        print(f"  Known cell: {data.knowncell.value}")
        print(f"  Known wavelength: {data.knownwavelength.value}")

    def test_overwrite_on_load(self):
        """Test that loadFile() can be called multiple times."""
        data = CUnmergedDataContent()

        # Call with empty path first
        error1 = data.loadFile('')
        assert error1.count() == 0

        # Call again with empty path
        error2 = data.loadFile('')
        assert error2.count() == 0

        # No errors should occur
        assert True


class TestXdsVersusScalepackDetection:
    """An XDS file named .hkl must be read as XDS, not as scalepack.

    Both formats use the .hkl/.HKL extension, so dispatching on the extension
    alone sent every XDS file down the scalepack reader. That reader takes
    line 1 as "<nsyms> <space group>" without checking, so XDS's
    "!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=FALSE" was reported as
    a space group, with format 'sca' and no cell or wavelength at all. It
    raised nothing, so the try/except fallback around it could never fire.
    """

    def test_xds_ascii_hkl_is_detected_as_xds(self):
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_XDS_ASCII))

        assert str(data.format) == 'xds'
        assert str(data.merged) == 'unmerged'

    def test_xds_header_is_not_reported_as_a_space_group(self):
        """The specific regression: the header line leaking into spaceGroup."""
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_XDS_ASCII))

        space_group = str(data.spaceGroup)
        assert 'MERGE' not in space_group, space_group
        assert 'FRIEDEL' not in space_group, space_group
        assert '!' not in space_group, space_group
        # P1 at this stage of processing: unindexed, but a real answer.
        assert space_group.replace(' ', '') == 'P1', space_group

    def test_xds_cell_and_wavelength_are_recovered(self):
        """The scalepack path discarded both; gemmi reads them fine."""
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_XDS_ASCII))

        assert bool(data.knowncell)
        assert float(data.cell.a) == pytest.approx(57.0, abs=1.0)
        assert bool(data.knownwavelength)
        assert float(data.wavelength) == pytest.approx(0.9763, abs=0.001)

    def test_integrate_hkl_is_detected_as_xds(self):
        """INTEGRATE.HKL opens with !OUTPUT_FILE=, not !FORMAT=, so the sniff
        keys on the leading '!' rather than on any single keyword."""
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_XDS_INTEGRATE))

        assert str(data.format) == 'xds'
        assert bool(data.knowncell)
        # !SPACE_GROUP_NUMBER=20 in the file itself.
        assert str(data.spaceGroup).replace(' ', '') == 'C2221', str(data.spaceGroup)

    def test_correct_hkl_is_detected_as_xds(self):
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_XDS_CORRECT))

        assert str(data.format) == 'xds'
        assert bool(data.knownwavelength)

    def test_scalepack_still_reads_as_scalepack(self):
        """The sniff must not drag genuine scalepack files into the XDS path."""
        for path in (TEST_SCA_MERGED, TEST_SCA_UNMERGED):
            data = CUnmergedDataContent()
            data.loadFile(str(path))
            assert str(data.format) == 'sca', path

    def test_scalepack_unmerged_space_group_is_unharmed(self):
        data = CUnmergedDataContent()
        data.loadFile(str(TEST_SCA_UNMERGED))

        assert str(data.merged) == 'unmerged'
        assert str(data.spaceGroup).replace(' ', '') == 'C2221', str(data.spaceGroup)



if __name__ == '__main__':
    pytest.main([__file__, '-v'])
