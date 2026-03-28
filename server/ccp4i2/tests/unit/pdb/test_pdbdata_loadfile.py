# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Tests for CPdbData.loadFile() implementation.

Tests the gemmi-based PDB/mmCIF file loading.
"""

import pytest
from pathlib import Path
from ccp4i2.core.CCP4ModelData import CPdbData, CPdbDataFile
from ccp4i2.core.base_object.error_reporting import CErrorReport
from ccp4i2.core import CCP4Utils

# Get the path to test data
CCP4I2_ROOT = Path(CCP4Utils.getCCP4I2Dir())
TEST_PDB = CCP4I2_ROOT / "wrappers/chainsaw/test_data/1a80_A.pdb"


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
        not TEST_PDB.exists(),
        reason="Test PDB file not available"
    )
    def test_loadfile_real_pdb(self):
        """Test loading a real PDB file."""
        test_pdb = str(TEST_PDB)

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
