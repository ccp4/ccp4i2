"""
Tests for CMtzData.loadFile() implementation.

Tests the gemmi-based MTZ file loading with proper CData setters.
"""

import pytest
from pathlib import Path
from core.CCP4XtalData import CMtzData, CMtzDataFile
from core.base_object.error_reporting import CErrorReport


class TestCMtzDataLoadFile:
    """Test CMtzData.loadFile() method."""

    def test_loadfile_nonexistent_file(self):
        """Test loading a non-existent MTZ file returns error."""
        mtz_data = CMtzData()
        error = mtz_data.loadFile('/nonexistent/file.mtz')

        assert isinstance(error, CErrorReport)
        assert error.count() > 0
        assert error.maxSeverity() > 0

    def test_loadfile_empty_path(self):
        """Test loading with empty path returns empty error report."""
        mtz_data = CMtzData()
        error = mtz_data.loadFile('')

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    def test_loadfile_none_path(self):
        """Test loading with None path returns empty error report."""
        mtz_data = CMtzData()
        error = mtz_data.loadFile(None)

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    @pytest.mark.skipif(
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz').exists(),
        reason="Test MTZ file not available"
    )
    def test_loadfile_real_mtz(self):
        """Test loading a real MTZ file extracts metadata correctly."""
        test_mtz = '/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz'

        mtz_data = CMtzData()
        error = mtz_data.loadFile(test_mtz)

        # Should load without errors
        if error.count() > 0:
            print(f"\nError loading MTZ: {error}")
        assert error.count() == 0

        # Check that metadata was extracted
        # Check child objects have data
        assert mtz_data.cell is not None
        assert mtz_data.spaceGroup is not None
        assert mtz_data.listOfColumns is not None

        # Verify gemmi object is stored
        assert hasattr(mtz_data, '_gemmi_mtz')
        assert mtz_data._gemmi_mtz is not None

        # Cell should have proper structure - access as attributes
        cell = mtz_data.cell
        assert cell.a.value > 0  # Should have positive values
        assert cell.b.value > 0
        assert cell.c.value > 0
        assert cell.alpha.value > 0
        assert cell.beta.value > 0
        assert cell.gamma.value > 0

        # Space group should be a string
        sg = mtz_data.spaceGroup.value
        assert isinstance(sg, str)
        assert len(sg) > 0

        # List of columns should be non-empty
        # Note: CList got replaced with plain list due to smart assignment
        columns = mtz_data.listOfColumns
        assert isinstance(columns, list)
        assert len(columns) > 0

        # Print some data for verification
        print(f"\nLoaded MTZ data:")
        print(f"  Space group: {sg}")
        print(f"  Cell: a={cell.a.value:.2f}, b={cell.b.value:.2f}, c={cell.c.value:.2f}")
        print(f"  Resolution: {mtz_data.resolutionRange.low.value:.2f} - {mtz_data.resolutionRange.high.value:.2f} Å")
        print(f"  Columns ({len(columns)}): {columns[:5]}...")  # First 5 columns

    @pytest.mark.skipif(
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz').exists(),
        reason="Test MTZ file not available"
    )
    def test_cdatafile_loadfile_integration(self):
        """Test CDataFile.loadFile() correctly instantiates CMtzData and loads."""
        test_mtz = '/Users/nmemn/Developer/ccp4i2/wrappers/pointless/test_data/brap_pk_6A.mtz'

        # Create a CMtzDataFile
        mtz_file = CMtzDataFile()
        mtz_file.setFullPath(test_mtz)

        # Load file using base CDataFile.loadFile()
        error = mtz_file.loadFile()

        # Should load without errors
        assert error.count() == 0

        # Content should be instantiated
        assert mtz_file.content is not None
        assert isinstance(mtz_file.content, CMtzData)

        # Content should be populated
        assert mtz_file.content.cell is not None
        assert mtz_file.content.spaceGroup is not None

        # Verify data was actually loaded
        cell = mtz_file.content.cell
        assert cell.a.value > 0

    def test_overwrite_on_load(self):
        """Test that loadFile() can be called multiple times."""
        mtz_data = CMtzData()

        # Call with empty path first
        error1 = mtz_data.loadFile('')
        assert error1.count() == 0

        # Call again with empty path
        error2 = mtz_data.loadFile('')
        assert error2.count() == 0

        # No errors should occur
        assert True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
