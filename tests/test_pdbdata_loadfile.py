"""
Tests for CPdbData.loadFile() implementation.

Tests the gemmi-based PDB/mmCIF file loading.
"""

import pytest
from pathlib import Path
from core.CCP4ModelData import CPdbData, CPdbDataFile
from core.base_object.error_reporting import CErrorReport


class TestCPdbDataLoadFile:
    """Test CPdbData.loadFile() method."""

    def test_loadfile_nonexistent_file(self):
        """Test loading a non-existent PDB file returns error."""
        pdb_data = CPdbData()
        error = pdb_data.loadFile('/nonexistent/file.pdb')

        assert isinstance(error, CErrorReport)
        assert error.count() > 0
        assert error.maxSeverity() > 0

    def test_loadfile_empty_path(self):
        """Test loading with empty path returns empty error report."""
        pdb_data = CPdbData()
        error = pdb_data.loadFile('')

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    def test_loadfile_none_path(self):
        """Test loading with None path returns empty error report."""
        pdb_data = CPdbData()
        error = pdb_data.loadFile(None)

        assert isinstance(error, CErrorReport)
        assert error.count() == 0

    @pytest.mark.skipif(
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/chainsaw/test_data/1a80_A.pdb').exists(),
        reason="Test PDB file not available"
    )
    def test_loadfile_real_pdb(self):
        """Test loading a real PDB file."""
        test_pdb = '/Users/nmemn/Developer/ccp4i2/wrappers/chainsaw/test_data/1a80_A.pdb'

        pdb_data = CPdbData()
        error = pdb_data.loadFile(test_pdb)

        # Should load without errors
        if error.count() > 0:
            print(f"\nError loading PDB: {error}")
        assert error.count() == 0

        # Verify gemmi object is stored
        assert hasattr(pdb_data, '_gemmi_structure')
        assert pdb_data._gemmi_structure is not None

        # Structure should have models
        structure = pdb_data._gemmi_structure
        assert len(structure) > 0  # At least one model

        # Print some data for verification
        print(f"\nLoaded PDB data:")
        print(f"  Name: {structure.name}")
        print(f"  Models: {len(structure)}")
        if len(structure) > 0:
            model = structure[0]
            print(f"  Chains in first model: {len(model)}")
            if len(model) > 0:
                chain = model[0]
                print(f"  First chain ID: {chain.name}")
                print(f"  Residues in first chain: {len(chain)}")

    @pytest.mark.skipif(
        not Path('/Users/nmemn/Developer/ccp4i2/wrappers/chainsaw/test_data/1a80_A.pdb').exists(),
        reason="Test PDB file not available"
    )
    def test_cdatafile_loadfile_integration(self):
        """Test CDataFile.loadFile() correctly instantiates CPdbData and loads."""
        test_pdb = '/Users/nmemn/Developer/ccp4i2/wrappers/chainsaw/test_data/1a80_A.pdb'

        # Create a CPdbDataFile
        pdb_file = CPdbDataFile()
        pdb_file.setFullPath(test_pdb)

        # Load file using base CDataFile.loadFile()
        error = pdb_file.loadFile()

        # Should load without errors
        assert error.count() == 0

        # Content should be instantiated
        assert pdb_file.content is not None
        assert isinstance(pdb_file.content, CPdbData)

        # Content should have structure loaded
        assert hasattr(pdb_file.content, '_gemmi_structure')
        assert pdb_file.content._gemmi_structure is not None

        # Verify structure has data
        structure = pdb_file.content._gemmi_structure
        assert len(structure) > 0

    def test_overwrite_on_load(self):
        """Test that loadFile() can be called multiple times."""
        pdb_data = CPdbData()

        # Call with empty path first
        error1 = pdb_data.loadFile('')
        assert error1.count() == 0

        # Call again with empty path
        error2 = pdb_data.loadFile('')
        assert error2.count() == 0

        # No errors should occur
        assert True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
