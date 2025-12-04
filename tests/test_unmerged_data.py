"""
Tests for CUnmergedDataContent.loadFile() implementation.

Tests loading of unmerged reflection data from multiple formats:
- MTZ (merged and unmerged)
- Scalepack (.sca - merged and unmerged formats)
- XDS (INTEGRATE.HKL, XDS_ASCII.HKL)
"""

import pytest
from pathlib import Path
from core.CCP4XtalData import CUnmergedDataContent, CUnmergedDataFile
from core.base_object.error_reporting import CErrorReport


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
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz').exists(),
        reason="Test MTZ file not available"
    )
    def test_loadfile_unmerged_mtz(self):
        """Test loading an unmerged MTZ file."""
        test_mtz = '/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz'

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
        not Path('/Users/nmemn/Developer/ccp4i2/demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled.sca').exists(),
        reason="Test Scalepack merged file not available"
    )
    def test_loadfile_scalepack_merged(self):
        """Test loading a merged Scalepack file."""
        test_sca = '/Users/nmemn/Developer/ccp4i2/demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled.sca'

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
        not Path('/Users/nmemn/Developer/ccp4i2/demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled_unmerged.sca').exists(),
        reason="Test Scalepack unmerged file not available"
    )
    def test_loadfile_scalepack_unmerged(self):
        """Test loading an unmerged Scalepack file."""
        test_sca = '/Users/nmemn/Developer/ccp4i2/demo_data/baz2b/BAZ2BA_x839.xia2/3daii-run/DataFiles/nt5073v16_xBAZ2BAx8392_scaled_unmerged.sca'

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

    @pytest.mark.skipif(
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz').exists(),
        reason="Test MTZ file not available"
    )
    def test_cdatafile_loadfile_integration(self):
        """Test CDataFile.loadFile() correctly instantiates CUnmergedDataContent and loads."""
        test_mtz = '/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz'

        # Create a CUnmergedDataFile
        data_file = CUnmergedDataFile()
        data_file.setFullPath(test_mtz)

        # Load file using base CDataFile.loadFile()
        error = data_file.loadFile()

        # Should load without errors
        assert error.count() == 0

        # Content should be instantiated
        assert data_file.content is not None
        assert isinstance(data_file.content, CUnmergedDataContent)

        # Content should be populated
        assert data_file.content.cell is not None
        assert data_file.content.spaceGroup is not None

        # Verify data was actually loaded
        cell = data_file.content.cell
        assert cell.a.value > 0

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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
