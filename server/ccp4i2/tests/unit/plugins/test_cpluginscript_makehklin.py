"""Integration tests for CPluginScript makeHklin methods."""

import pytest
from pathlib import Path
import tempfile
import shutil

try:
    import gemmi
    import numpy as np
    GEMMI_AVAILABLE = True
except ImportError:
    GEMMI_AVAILABLE = False

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.base_object.base_classes import CContainer
from ccp4i2.core.base_object.fundamental_types import CInt
from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_OK

# Skip all tests if gemmi is not installed
pytestmark = pytest.mark.skipif(not GEMMI_AVAILABLE, reason="gemmi not installed")


# ============================================================================
# Mock File Classes (Minimal implementations for testing)
# ============================================================================

class MockMiniMtzDataFile:
    """Mock CMiniMtzDataFile for testing."""

    CONTENT_SIGNATURE_LIST = []  # Override in subclasses

    def __init__(self, path: Path, content_flag: int = 1):
        self._path = path
        self.contentFlag = CInt(content_flag)

    def getFullPath(self) -> Path:
        return self._path


class MockCObsDataFile(MockMiniMtzDataFile):
    """Mock CObsDataFile with CONTENT_SIGNATURE_LIST."""

    # Content flag constants (from stub)
    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    CONTENT_ANNOTATION = [
        'Anomalous Is',
        'Anomalous SFs',
        'Mean Is',
        'Mean SFs'
    ]

    CONTENT_SIGNATURE_LIST = [
        ['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'],
        ['Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'],
        ['I', 'SIGI'],
        ['F', 'SIGF']
    ]


class MockCFreeRDataFile(MockMiniMtzDataFile):
    """Mock CFreeRDataFile with CONTENT_SIGNATURE_LIST."""

    CONTENT_FLAG_FREER = 1

    CONTENT_ANNOTATION = ['FreeR']

    CONTENT_SIGNATURE_LIST = [['FREER']]


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    tmp = tempfile.mkdtemp()
    yield Path(tmp)
    shutil.rmtree(tmp)


@pytest.fixture
def obs_mtz_file(temp_dir):
    """Create a sample MTZ file with F, SIGF columns (FMEAN format)."""
    if not GEMMI_AVAILABLE:
        return None

    mtz_path = temp_dir / "obs.mtz"

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    mtz.add_dataset("test_dataset")
    mtz.add_column("F", "F")
    mtz.add_column("SIGF", "Q")

    # Add realistic resolution data (100-6 Å range)
    data = np.array([
        [1, 0, 0, 100.0, 5.0],    # d ≈ 50 Å
        [0, 1, 0, 150.0, 7.0],    # d ≈ 60 Å
        [5, 6, 7, 200.0, 10.0],   # d ≈ 8 Å
        [8, 10, 11, 180.0, 9.0],  # d ≈ 6 Å
    ], dtype=np.float32)
    mtz.set_data(data)

    mtz.update_reso()
    mtz.write_to_file(str(mtz_path))

    return mtz_path


@pytest.fixture
def freer_mtz_file(temp_dir):
    """Create a sample MTZ file with FREER column."""
    if not GEMMI_AVAILABLE:
        return None

    mtz_path = temp_dir / "freer.mtz"

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
    mtz.add_dataset("test_dataset")
    mtz.add_column("FREER", "I")

    # Same reflections as obs_mtz_file
    data = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 1],
        [5, 6, 7, 0],
        [8, 10, 11, 1],
    ], dtype=np.float32)
    mtz.set_data(data)

    mtz.update_reso()
    mtz.write_to_file(str(mtz_path))

    return mtz_path


@pytest.fixture
def plugin_script(temp_dir, obs_mtz_file, freer_mtz_file):
    """Create a CPluginScript instance with mock file objects."""
    script = CPluginScript()
    script.workDirectory = temp_dir

    # Create input container
    script.container.inputData = CContainer(name="inputData")

    # Add mock file objects to inputData
    obs_file = MockCObsDataFile(obs_mtz_file, content_flag=4)  # FMEAN
    script.container.inputData.HKLIN1 = obs_file

    freer_file = MockCFreeRDataFile(freer_mtz_file, content_flag=1)
    script.container.inputData.FREERFLAG = freer_file

    # Create output container
    script.container.outputData = CContainer(name="outputData")

    return script


# ============================================================================
# Tests for makeHklinGemmi
# ============================================================================

class TestMakeHklinGemmi:
    """Tests for makeHklinGemmi method (new Pythonic API)."""

    def test_simple_merge_with_string_names(self, plugin_script, temp_dir):
        """Test basic merge using simple string names."""
        result = plugin_script.makeHklinGemmi(['HKLIN1', 'FREERFLAG'])

        # Verify output file exists
        assert result.exists()
        assert result == temp_dir / "hklin.mtz"

        # Read and verify merged file
        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Should have H, K, L + prefixed columns from each file
        # HKLIN1_F, HKLIN1_SIGF (from HKLIN1 contentFlag=4)
        # FREERFLAG_FREER (from FREERFLAG contentFlag=1)
        assert 'H' in column_labels
        assert 'K' in column_labels
        assert 'L' in column_labels
        assert 'HKLIN1_F' in column_labels
        assert 'HKLIN1_SIGF' in column_labels
        assert 'FREERFLAG_FREER' in column_labels

        # Check reflections count (merge creates complete reflection set)
        assert merged_mtz.nreflections >= 4  # At least the 4 we input

    def test_merge_with_rename(self, plugin_script, temp_dir):
        """Test merge with column renaming."""
        result = plugin_script.makeHklinGemmi([
            {
                'name': 'HKLIN1',
                'rename': {'F': 'F_native', 'SIGF': 'SIGF_native'}
            },
            'FREERFLAG'
        ])

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Should have renamed columns from HKLIN1 (explicit rename)
        # and prefixed column from FREERFLAG (default naming)
        assert 'F_native' in column_labels
        assert 'SIGF_native' in column_labels
        assert 'F' not in column_labels  # Original name should not exist
        assert 'SIGF' not in column_labels
        assert 'FREERFLAG_FREER' in column_labels  # Default prefix

    def test_custom_output_name(self, plugin_script, temp_dir):
        """Test custom output file name."""
        result = plugin_script.makeHklinGemmi(
            ['HKLIN1'],
            output_name='custom_output'
        )

        assert result == temp_dir / "custom_output.mtz"
        assert result.exists()

    def test_merge_strategy_error(self, plugin_script, temp_dir):
        """Test error strategy on column conflict."""
        # Add same file twice - should conflict on F, SIGF columns
        from ccp4i2.core.CCP4Utils import MtzMergeError

        with pytest.raises((ValueError, MtzMergeError), match="conflict"):
            plugin_script.makeHklinGemmi(
                ['HKLIN1', 'HKLIN1'],
                merge_strategy='error'
            )

    def test_merge_strategy_rename(self, plugin_script, temp_dir):
        """Test rename strategy auto-renames conflicts."""
        result = plugin_script.makeHklinGemmi(
            ['HKLIN1', 'HKLIN1', 'HKLIN1'],
            merge_strategy='rename'
        )

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Should have HKLIN1_F, HKLIN1_F_1, HKLIN1_F_2 (auto-renamed on conflict)
        assert 'HKLIN1_F' in column_labels
        assert 'HKLIN1_F_1' in column_labels
        assert 'HKLIN1_F_2' in column_labels

    def test_file_not_found_error(self, plugin_script):
        """Test error when file object not found."""
        with pytest.raises(AttributeError, match="not found"):
            plugin_script.makeHklinGemmi(['NONEXISTENT'])

    def test_no_content_signature_list_error(self, plugin_script, temp_dir):
        """Test error when file has no CONTENT_SIGNATURE_LIST."""
        # Create a mock file without CONTENT_SIGNATURE_LIST
        class BadFile:
            def __init__(self):
                self.contentFlag = CInt(1)
            def getFullPath(self):
                return temp_dir / "bad.mtz"

        plugin_script.container.inputData.BAD = BadFile()

        with pytest.raises(ValueError, match="no CONTENT_SIGNATURE_LIST"):
            plugin_script.makeHklinGemmi(['BAD'])

    def test_invalid_content_flag_error(self, plugin_script, temp_dir, obs_mtz_file):
        """Test error with invalid contentFlag."""
        # Create file with invalid contentFlag
        bad_file = MockCObsDataFile(obs_mtz_file, content_flag=999)
        plugin_script.container.inputData.BAD_FLAG = bad_file

        with pytest.raises(ValueError, match="Invalid contentFlag"):
            plugin_script.makeHklinGemmi(['BAD_FLAG'])

    def test_no_path_error(self, plugin_script, temp_dir):
        """Test error when file has no path set."""
        # Create file that returns empty path
        class NoPathFile(MockCObsDataFile):
            def getFullPath(self):
                return None

        plugin_script.container.inputData.NO_PATH = NoPathFile(temp_dir / "fake.mtz")

        with pytest.raises(ValueError, match="has no path set"):
            plugin_script.makeHklinGemmi(['NO_PATH'])


# ============================================================================
# Tests for makeHklin (backward-compatible API)
# ============================================================================

class TestMakeHklin:
    """Tests for makeHklin method (legacy API)."""

    def test_simple_merge_with_string_names(self, plugin_script, temp_dir):
        """Test basic merge using simple string names (legacy API)."""
        hklin_filename, error = plugin_script.makeHklin(['HKLIN1', 'FREERFLAG'])

        # Should succeed with no errors
        assert error.maxSeverity() == SEVERITY_OK

        # Verify output file exists
        output_path = temp_dir / "hklin.mtz"
        assert output_path.exists()

        # Read and verify
        merged_mtz = gemmi.read_mtz_file(str(output_path))
        column_labels = [col.label for col in merged_mtz.columns]

        # Legacy API uses identity mapping (no prefixing) for backward compatibility
        assert 'F' in column_labels
        assert 'SIGF' in column_labels
        assert 'FREER' in column_labels

    def test_custom_output_name(self, plugin_script, temp_dir):
        """Test custom output name (legacy API)."""
        hklin_filename, error = plugin_script.makeHklin(['HKLIN1'], hklin='custom')

        assert error.maxSeverity() == SEVERITY_OK
        output_path = temp_dir / "custom.mtz"
        assert output_path.exists()

    def test_override_content_flag_triggers_conversion(self, plugin_script, temp_dir, obs_mtz_file):
        """Test that contentFlag override triggers conversion when flags differ."""
        # Create a second obs file with FPAIR columns but contentFlag=4 (FMEAN)
        fpair_mtz_path = temp_dir / "fpair.mtz"

        # Create MTZ with Fplus, SIGFplus, Fminus, SIGFminus
        mtz = gemmi.Mtz(with_base=True)
        mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
        mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))
        mtz.add_dataset("test")
        mtz.add_column("Fplus", "G")
        mtz.add_column("SIGFplus", "L")
        mtz.add_column("Fminus", "G")
        mtz.add_column("SIGFminus", "L")

        data = np.array([
            [1, 0, 0, 105.0, 6.0, 95.0, 6.0],
            [0, 1, 0, 155.0, 8.0, 145.0, 8.0],
            [0, 0, 1, 205.0, 11.0, 195.0, 11.0],
        ], dtype=np.float32)
        mtz.set_data(data)
        mtz.update_reso()
        mtz.write_to_file(str(fpair_mtz_path))

        # Add as HKLIN2 with contentFlag=4 (FMEAN), but file actually has FPAIR data
        # Request contentFlag=2 (FPAIR) - should trigger conversion
        hklin2 = MockCObsDataFile(fpair_mtz_path, content_flag=4)
        plugin_script.container.inputData.HKLIN2 = hklin2

        # Merge with explicit contentFlag override that differs from file's flag
        hklin_filename, error = plugin_script.makeHklin([
            'HKLIN1',
            ['HKLIN2', MockCObsDataFile.CONTENT_FLAG_FPAIR]  # Request FPAIR (flag=2), but file has flag=4
        ])

        # Should get error code 209 (NotImplementedError from conversion)
        assert error.maxSeverity() > SEVERITY_OK
        assert error.count() > 0

        # Check that conversion was attempted
        error_details = error.report()
        assert 'not yet implemented' in error_details.lower() or 'conversion' in error_details.lower()

    def test_file_not_found_error(self, plugin_script):
        """Test error handling when file not found."""
        hklin_filename, error = plugin_script.makeHklin(['NONEXISTENT'])

        assert error.maxSeverity() > SEVERITY_OK
        assert error.count() > 0

    def test_invalid_content_flag_error(self, plugin_script):
        """Test error with invalid contentFlag in [name, flag] syntax."""
        hklin_filename, error = plugin_script.makeHklin([
            ['HKLIN1', 999]  # Invalid flag
        ])

        assert error.maxSeverity() > SEVERITY_OK

    def test_invalid_item_format_error(self, plugin_script):
        """Test error with invalid item format."""
        hklin_filename, error = plugin_script.makeHklin([
            {'invalid': 'dict'}  # Invalid format
        ])

        assert error.maxSeverity() > SEVERITY_OK

    def test_no_conversion_when_flags_match(self, plugin_script, temp_dir):
        """Test that no conversion occurs when contentFlag matches target."""
        # Store original contentFlag (4 = FMEAN)
        original_flag = plugin_script.container.inputData.HKLIN1.contentFlag

        # Request same contentFlag as file has - no conversion needed
        hklin_filename, error = plugin_script.makeHklin([
            ['HKLIN1', MockCObsDataFile.CONTENT_FLAG_FMEAN]  # Same as file's flag
        ])

        assert error.maxSeverity() == SEVERITY_OK

        # Verify original flag unchanged (no conversion occurred)
        assert plugin_script.container.inputData.HKLIN1.contentFlag == original_flag

        # Verify output exists
        output_path = temp_dir / "hklin.mtz"
        assert output_path.exists()


# ============================================================================
# Edge Case Tests
# ============================================================================

class TestMakeHklinEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_file_list(self, plugin_script):
        """Test error with empty file list."""
        from ccp4i2.core.CCP4Utils import MtzMergeError

        with pytest.raises((ValueError, MtzMergeError)):
            plugin_script.makeHklinGemmi([])

    def test_single_file_merge(self, plugin_script, temp_dir):
        """Test merging a single file (should work)."""
        result = plugin_script.makeHklinGemmi(['HKLIN1'])

        assert result.exists()

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Default behavior prefixes columns with name
        assert 'HKLIN1_F' in column_labels
        assert 'HKLIN1_SIGF' in column_labels

    def test_partial_rename(self, plugin_script, temp_dir):
        """Test renaming only some columns."""
        result = plugin_script.makeHklinGemmi([
            {
                'name': 'HKLIN1',
                'rename': {'F': 'F_obs'}  # Only rename F, not SIGF
            }
        ])

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # F gets explicit rename, SIGF gets default prefix
        assert 'F_obs' in column_labels
        assert 'HKLIN1_SIGF' in column_labels  # Default prefix for unrenamed column
        assert 'F' not in column_labels
