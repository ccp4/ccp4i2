"""
Test CAsuDataFile digest functionality.

Tests the digest_casudatafile_file_object function which provides
sequence selection information for the frontend.
"""

import json
import logging
from pathlib import Path

import pytest

from ccp4i2.core.CCP4ModelData import CAsuDataFile

logger = logging.getLogger(f"ccp4x::{__name__}")


class TestCAsuDataFileDigest:
    """Test CAsuDataFile digest functionality."""

    @pytest.fixture
    def gamma_asu_file(self):
        """Path to gamma demo ASU file."""
        asu_path = Path(__file__).parent.parent.parent.parent.parent / "demo_data" / "gamma" / "gamma.asu.xml"
        if not asu_path.exists():
            pytest.skip(f"Demo data not found: {asu_path}")
        return str(asu_path)

    def test_casudatafile_loads_correctly(self, gamma_asu_file):
        """Test that CAsuDataFile can load the gamma ASU file."""
        asu_file = CAsuDataFile()
        asu_file.setFullPath(gamma_asu_file)

        # Verify file is set
        assert asu_file.isSet(), "File should be set after setFullPath"

        # Load file content
        asu_file.loadFile()
        contents = asu_file.getFileContent()

        # Verify content loaded
        assert contents is not None, "File content should not be None after loadFile"
        assert hasattr(contents, 'seqList'), "Content should have seqList attribute"
        assert len(contents.seqList) == 1, "Gamma ASU should have 1 sequence"

    def test_digest_casudatafile_returns_sequences(self, gamma_asu_file):
        """Test that digest function returns sequence information."""
        from ccp4x.lib.utils.files.digest import digest_casudatafile_file_object

        asu_file = CAsuDataFile()
        asu_file.setFullPath(gamma_asu_file)

        result = digest_casudatafile_file_object(asu_file)

        # Check structure
        assert "sequences" in result, f"Digest should contain 'sequences' key, got: {result.keys()}"
        assert "sequenceCount" in result, "Digest should contain 'sequenceCount' key"

        # Check sequence count
        assert result["sequenceCount"] == 1, f"Expected 1 sequence, got {result['sequenceCount']}"
        assert len(result["sequences"]) == 1, f"Expected 1 sequence in list"

        # Check sequence details
        seq = result["sequences"][0]
        assert seq["index"] == 0, "First sequence should have index 0"
        assert seq["name"] == "Gamma", f"Expected name 'Gamma', got '{seq['name']}'"
        assert seq["polymerType"] == "PROTEIN", f"Expected polymerType 'PROTEIN', got '{seq['polymerType']}'"
        assert seq["nCopies"] == 1, f"Expected nCopies 1, got {seq['nCopies']}"
        assert seq["sequenceLength"] == 135, f"Expected sequenceLength 135, got {seq['sequenceLength']}"
        assert seq["selected"] is True, "Sequence should be selected by default"

        logger.info(f"Digest result: {json.dumps(result, indent=2)}")

    def test_digest_casudatafile_with_empty_file(self):
        """Test that digest handles unset file gracefully."""
        from ccp4x.lib.utils.files.digest import digest_casudatafile_file_object

        asu_file = CAsuDataFile()
        # Don't set any path

        result = digest_casudatafile_file_object(asu_file)

        # Should return failure status
        assert result.get("status") == "Failed", f"Expected 'Failed' status for unset file, got: {result}"
        assert "reason" in result, "Should include failure reason"

    def test_digest_casudatafile_sequence_selection(self, gamma_asu_file):
        """Test that digest respects selection CDict values."""
        from ccp4x.lib.utils.files.digest import digest_casudatafile_file_object

        asu_file = CAsuDataFile()
        asu_file.setFullPath(gamma_asu_file)

        # Set a sequence as not selected
        asu_file.selection["Gamma"] = False

        result = digest_casudatafile_file_object(asu_file)

        # Check that selection is reflected
        assert len(result["sequences"]) == 1
        seq = result["sequences"][0]
        assert seq["selected"] is False, "Sequence should respect selection CDict value"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
