"""
Unit tests for MTZ utility functions.

Tests the reusable MTZ column extraction utilities.
"""
import pytest
from pathlib import Path


def test_extract_columns_basic(tmp_path):
    """Test basic column extraction from MTZ file."""
    from core.conversions.mtz_utils import extract_columns
    import gemmi
    import numpy as np

    # Create a simple test MTZ file
    mtz_in = gemmi.Mtz()
    mtz_in.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz_in.cell.set(10, 10, 10, 90, 90, 90)

    mtz_in.add_dataset('crystal')
    mtz_in.add_column('H', 'H')
    mtz_in.add_column('K', 'H')
    mtz_in.add_column('L', 'H')
    mtz_in.add_column('F', 'F')
    mtz_in.add_column('SIGF', 'Q')
    mtz_in.add_column('I', 'J')
    mtz_in.add_column('SIGI', 'Q')

    # Create some test data
    data = np.array([
        [1, 0, 0, 100.0, 5.0, 1000.0, 50.0],
        [0, 1, 0, 200.0, 10.0, 2000.0, 100.0],
        [0, 0, 1, 300.0, 15.0, 3000.0, 150.0],
    ])
    mtz_in.set_data(data)

    # Write test MTZ
    input_mtz = tmp_path / "test_input.mtz"
    mtz_in.write_to_file(str(input_mtz))

    # Extract just F and SIGF columns
    output_mtz = tmp_path / "test_output.mtz"
    extract_columns(
        str(input_mtz),
        str(output_mtz),
        [('F', 'F', 'F'), ('SIGF', 'SIGF', 'Q')]
    )

    # Verify output
    assert output_mtz.exists(), "Output MTZ not created"

    mtz_out = gemmi.read_mtz_file(str(output_mtz))
    column_labels = [col.label for col in mtz_out.columns]

    assert 'H' in column_labels
    assert 'K' in column_labels
    assert 'L' in column_labels
    assert 'F' in column_labels
    assert 'SIGF' in column_labels
    assert 'I' not in column_labels, "I column should not be in output"
    assert 'SIGI' not in column_labels, "SIGI column should not be in output"


def test_extract_columns_rename(tmp_path):
    """Test column extraction with renaming."""
    from core.conversions.mtz_utils import extract_columns
    import gemmi
    import numpy as np

    # Create test MTZ with F(+)/F(-) columns
    mtz_in = gemmi.Mtz()
    mtz_in.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz_in.cell.set(10, 10, 10, 90, 90, 90)

    mtz_in.add_dataset('crystal')
    mtz_in.add_column('H', 'H')
    mtz_in.add_column('K', 'H')
    mtz_in.add_column('L', 'H')
    mtz_in.add_column('F(+)', 'G')
    mtz_in.add_column('SIGF(+)', 'L')

    data = np.array([
        [1, 0, 0, 100.0, 5.0],
        [0, 1, 0, 200.0, 10.0],
    ])
    mtz_in.set_data(data)

    input_mtz = tmp_path / "test_fpair.mtz"
    mtz_in.write_to_file(str(input_mtz))

    # Extract with renaming F(+) → Fplus
    output_mtz = tmp_path / "test_renamed.mtz"
    extract_columns(
        str(input_mtz),
        str(output_mtz),
        [('F(+)', 'Fplus', 'G'), ('SIGF(+)', 'SIGFplus', 'L')]
    )

    # Verify renamed columns
    mtz_out = gemmi.read_mtz_file(str(output_mtz))
    column_labels = [col.label for col in mtz_out.columns]

    assert 'Fplus' in column_labels
    assert 'SIGFplus' in column_labels
    assert 'F(+)' not in column_labels


def test_extract_columns_missing_column(tmp_path):
    """Test that missing column raises error."""
    from core.conversions.mtz_utils import extract_columns
    import gemmi
    import numpy as np

    # Create simple test MTZ
    mtz_in = gemmi.Mtz()
    mtz_in.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz_in.cell.set(10, 10, 10, 90, 90, 90)

    mtz_in.add_dataset('crystal')
    mtz_in.add_column('H', 'H')
    mtz_in.add_column('K', 'H')
    mtz_in.add_column('L', 'H')
    mtz_in.add_column('F', 'F')

    data = np.array([[1, 0, 0, 100.0]])
    mtz_in.set_data(data)

    input_mtz = tmp_path / "test_missing.mtz"
    mtz_in.write_to_file(str(input_mtz))

    output_mtz = tmp_path / "test_output.mtz"

    # Try to extract non-existent column
    with pytest.raises(RuntimeError, match="not found"):
        extract_columns(
            str(input_mtz),
            str(output_mtz),
            [('MISSING_COL', 'F', 'F')]
        )


def test_get_column_labels(tmp_path):
    """Test getting column labels from MTZ."""
    from core.conversions.mtz_utils import get_column_labels
    import gemmi
    import numpy as np

    # Create test MTZ
    mtz = gemmi.Mtz()
    mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz.cell.set(10, 10, 10, 90, 90, 90)

    mtz.add_dataset('crystal')
    mtz.add_column('H', 'H')
    mtz.add_column('K', 'H')
    mtz.add_column('L', 'H')
    mtz.add_column('F', 'F')
    mtz.add_column('SIGF', 'Q')

    data = np.array([[1, 0, 0, 100.0, 5.0]])
    mtz.set_data(data)

    mtz_file = tmp_path / "test.mtz"
    mtz.write_to_file(str(mtz_file))

    # Get labels
    labels = get_column_labels(str(mtz_file))

    assert labels == ['H', 'K', 'L', 'F', 'SIGF']


def test_has_columns(tmp_path):
    """Test checking for required columns."""
    from core.conversions.mtz_utils import has_columns
    import gemmi
    import numpy as np

    # Create test MTZ
    mtz = gemmi.Mtz()
    mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz.cell.set(10, 10, 10, 90, 90, 90)

    mtz.add_dataset('crystal')
    mtz.add_column('H', 'H')
    mtz.add_column('K', 'H')
    mtz.add_column('L', 'H')
    mtz.add_column('F', 'F')
    mtz.add_column('SIGF', 'Q')

    data = np.array([[1, 0, 0, 100.0, 5.0]])
    mtz.set_data(data)

    mtz_file = tmp_path / "test.mtz"
    mtz.write_to_file(str(mtz_file))

    # Test has_columns
    assert has_columns(str(mtz_file), ['F', 'SIGF']) is True
    assert has_columns(str(mtz_file), ['F']) is True
    assert has_columns(str(mtz_file), ['F', 'SIGF', 'I']) is False
    assert has_columns(str(mtz_file), ['MISSING']) is False


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
