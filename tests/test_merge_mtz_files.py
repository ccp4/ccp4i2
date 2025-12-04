"""Tests for merge_mtz_files utility function."""

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

from core.CCP4Utils import merge_mtz_files, MtzMergeError


# Skip all tests if gemmi is not installed
pytestmark = pytest.mark.skipif(not GEMMI_AVAILABLE, reason="gemmi not installed")


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    tmp = tempfile.mkdtemp()
    yield Path(tmp)
    shutil.rmtree(tmp)


@pytest.fixture
def sample_mtz_1(temp_dir):
    """Create a sample MTZ file with F, SIGF columns."""
    if not GEMMI_AVAILABLE:
        return None

    mtz_path = temp_dir / "test1.mtz"

    # Create a minimal MTZ file
    mtz = gemmi.Mtz(with_base=True)

    # Set space group and cell
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))

    # Add a dataset
    mtz.add_dataset("test_dataset")

    # Add F and SIGF columns
    mtz.add_column("F", "F")
    mtz.add_column("SIGF", "Q")

    # Add some dummy data with realistic resolution (100-6 Å range)
    # With cell (50, 60, 70), these Miller indices give reasonable d-spacings
    data = np.array([
        [1, 0, 0, 100.0, 5.0],    # H K L F SIGF - d ≈ 50 Å
        [0, 1, 0, 150.0, 7.0],    # d ≈ 60 Å
        [5, 6, 7, 200.0, 10.0],   # d ≈ 8 Å
        [8, 10, 11, 180.0, 9.0],  # d ≈ 6 Å (high resolution)
    ], dtype=np.float32)
    mtz.set_data(data)

    mtz.update_reso()
    mtz.write_to_file(str(mtz_path))

    return mtz_path


@pytest.fixture
def sample_mtz_2(temp_dir):
    """Create a second sample MTZ file with FreeR_flag."""
    if not GEMMI_AVAILABLE:
        return None

    mtz_path = temp_dir / "test2.mtz"

    # Create a minimal MTZ file with same cell/spacegroup as sample_mtz_1
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name("P 21 21 21")
    mtz.set_cell_for_all(gemmi.UnitCell(50, 60, 70, 90, 90, 90))

    mtz.add_dataset("test_dataset")
    mtz.add_column("FreeR_flag", "I")

    # Same reflections as sample_mtz_1
    data = np.array([
        [1, 0, 0, 0],        # H K L FreeR
        [0, 1, 0, 1],
        [5, 6, 7, 0],
        [8, 10, 11, 1],
    ], dtype=np.float32)
    mtz.set_data(data)

    mtz.update_reso()
    mtz.write_to_file(str(mtz_path))

    return mtz_path


class TestMergeMtzFiles:
    """Tests for merge_mtz_files function."""

    def test_merge_two_files_simple(self, sample_mtz_1, sample_mtz_2, temp_dir):
        """Test basic merge of two MTZ files."""
        output_path = temp_dir / "merged.mtz"

        result = merge_mtz_files(
            input_specs=[
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F', 'SIGF': 'SIGF'}},
                {'path': str(sample_mtz_2), 'column_mapping': {'FreeR_flag': 'FreeR_flag'}}
            ],
            output_path=str(output_path),
            merge_strategy='first'
        )

        # Verify output file exists
        assert result.exists()
        assert result == output_path

        # Read and verify merged file
        merged_mtz = gemmi.read_mtz_file(str(result))

        # Should have H, K, L + F, SIGF, FreeR_flag = 6 columns
        column_labels = [col.label for col in merged_mtz.columns]
        assert 'H' in column_labels
        assert 'K' in column_labels
        assert 'L' in column_labels
        assert 'F' in column_labels
        assert 'SIGF' in column_labels
        assert 'FreeR_flag' in column_labels

        # Check that data is present (merge creates complete reflection set, so > input count)
        assert merged_mtz.nreflections > 0
        # Should have significantly more reflections than the 4 we input (full set for resolution range)
        assert merged_mtz.nreflections >= 4

        # Verify spacegroup and cell are preserved
        assert merged_mtz.spacegroup.hm == 'P 21 21 21'
        assert abs(merged_mtz.cell.a - 50.0) < 0.01

    def test_merge_with_rename(self, sample_mtz_1, temp_dir):
        """Test column renaming."""
        output_path = temp_dir / "renamed.mtz"

        result = merge_mtz_files(
            input_specs=[
                {
                    'path': str(sample_mtz_1),
                    'column_mapping': {'F': 'F_NAT', 'SIGF': 'SIGF_NAT'}
                }
            ],
            output_path=str(output_path)
        )

        # Verify renamed columns
        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        assert 'F_NAT' in column_labels
        assert 'SIGF_NAT' in column_labels
        assert 'F' not in column_labels  # Original name should not be present
        assert 'SIGF' not in column_labels

    def test_merge_strategy_error_on_conflict(self, sample_mtz_1, sample_mtz_2, temp_dir):
        """Test that conflict detection works with strategy='error'."""
        output_path = temp_dir / "conflict.mtz"

        # Add F column to both files - should raise error
        with pytest.raises(ValueError, match="Column conflict"):
            merge_mtz_files(
                input_specs=[
                    {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}},
                    {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}}  # Same file, same column
                ],
                output_path=str(output_path),
                merge_strategy='error'
            )

    def test_merge_strategy_first_keeps_first(self, sample_mtz_1, temp_dir):
        """Test that strategy='first' keeps the first occurrence."""
        output_path = temp_dir / "first.mtz"

        # Try to add F twice - should keep first
        result = merge_mtz_files(
            input_specs=[
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}},
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}}  # Duplicate
            ],
            output_path=str(output_path),
            merge_strategy='first'
        )

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Should only have F once (plus H, K, L)
        assert column_labels.count('F') == 1

    def test_merge_strategy_rename_auto_renames(self, sample_mtz_1, temp_dir):
        """Test that strategy='rename' auto-renames conflicts."""
        output_path = temp_dir / "auto_rename.mtz"

        result = merge_mtz_files(
            input_specs=[
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}},
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}},  # Conflict
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F'}}   # Another conflict
            ],
            output_path=str(output_path),
            merge_strategy='rename'
        )

        merged_mtz = gemmi.read_mtz_file(str(result))
        column_labels = [col.label for col in merged_mtz.columns]

        # Should have F, F_1, F_2
        assert 'F' in column_labels
        assert 'F_1' in column_labels
        assert 'F_2' in column_labels

    def test_file_not_found(self, temp_dir):
        """Test error handling for non-existent file."""
        with pytest.raises(FileNotFoundError):
            merge_mtz_files(
                input_specs=[
                    {'path': '/nonexistent/file.mtz', 'column_mapping': {'F': 'F'}}
                ],
                output_path=str(temp_dir / "output.mtz")
            )

    def test_column_not_found(self, sample_mtz_1, temp_dir):
        """Test error when requested column doesn't exist."""
        with pytest.raises(ValueError, match="Column.*not found"):
            merge_mtz_files(
                input_specs=[
                    {'path': str(sample_mtz_1), 'column_mapping': {'NONEXISTENT': 'NONEXISTENT'}}
                ],
                output_path=str(temp_dir / "output.mtz")
            )

    def test_empty_input_specs(self, temp_dir):
        """Test error with empty input specs."""
        with pytest.raises(ValueError, match="cannot be empty"):
            merge_mtz_files(
                input_specs=[],
                output_path=str(temp_dir / "output.mtz")
            )

    def test_missing_required_keys(self, sample_mtz_1, temp_dir):
        """Test error when spec is missing required keys."""
        with pytest.raises(ValueError, match="missing required key"):
            merge_mtz_files(
                input_specs=[
                    {'path': str(sample_mtz_1)}  # Missing 'columns'
                ],
                output_path=str(temp_dir / "output.mtz")
            )

    def test_preserves_spacegroup_and_cell(self, sample_mtz_1, temp_dir):
        """Test that space group and cell are preserved."""
        output_path = temp_dir / "preserved.mtz"

        result = merge_mtz_files(
            input_specs=[
                {'path': str(sample_mtz_1), 'column_mapping': {'F': 'F', 'SIGF': 'SIGF'}}
            ],
            output_path=str(output_path)
        )

        original_mtz = gemmi.read_mtz_file(str(sample_mtz_1))
        merged_mtz = gemmi.read_mtz_file(str(result))

        # Check space group
        assert merged_mtz.spacegroup.hm == original_mtz.spacegroup.hm

        # Check cell parameters
        assert merged_mtz.cell.a == pytest.approx(original_mtz.cell.a)
        assert merged_mtz.cell.b == pytest.approx(original_mtz.cell.b)
        assert merged_mtz.cell.c == pytest.approx(original_mtz.cell.c)
