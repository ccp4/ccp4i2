"""
MTZ utility functions for column extraction and manipulation.

This module provides reusable utilities for working with MTZ files,
particularly for extracting specific columns to create mini-MTZ files.
Similar to the gemmi-based make_hklin pattern.
"""

from typing import List, Tuple
import gemmi
import numpy as np


def extract_columns(
    input_mtz: str,
    output_mtz: str,
    column_spec: List[Tuple[str, str, str]]
) -> None:
    """
    Extract specific columns from an MTZ file to create a mini-MTZ.

    This is a gemmi-based replacement for splitHklout/sftools operations.
    Can be used from both CPluginScript context and standalone converters.

    Args:
        input_mtz: Path to input MTZ file
        output_mtz: Path to output MTZ file
        column_spec: List of (source_col, target_col, mtz_type) tuples
            - source_col: Column label in input MTZ (e.g., 'F', 'SIGF')
            - target_col: Column label in output MTZ (e.g., 'F', 'SIGF')
            - mtz_type: MTZ column type (e.g., 'F', 'Q', 'J', 'G', 'L')

    Example:
        >>> extract_columns(
        ...     'input.mtz',
        ...     'output.mtz',
        ...     [('F', 'F', 'F'), ('SIGF', 'SIGF', 'Q')]
        ... )

    Raises:
        RuntimeError: If source column not found in input MTZ
        FileNotFoundError: If input MTZ doesn't exist
    """
    # Read input MTZ
    mtz_in = gemmi.read_mtz_file(str(input_mtz))

    # Create output MTZ with same structure
    mtz_out = gemmi.Mtz()
    mtz_out.spacegroup = mtz_in.spacegroup
    mtz_out.cell = mtz_in.cell

    # Add dataset and HKL columns
    mtz_out.add_dataset('crystal')
    mtz_out.add_column('H', 'H')
    mtz_out.add_column('K', 'H')
    mtz_out.add_column('L', 'H')

    # Add data columns
    for target_col, mtz_type in [(spec[1], spec[2]) for spec in column_spec]:
        mtz_out.add_column(target_col, mtz_type)

    # Extract source columns from input
    miller_array = mtz_in.make_miller_array()
    data_arrays = []

    for source_col, target_col, mtz_type in column_spec:
        col = mtz_in.column_with_label(source_col)
        if col is None:
            raise RuntimeError(
                f"Column '{source_col}' not found in {input_mtz}")
        data_arrays.append(np.array(col))

    # Build output data array: [H, K, L, ...data_columns...]
    data = np.column_stack([miller_array] + data_arrays)

    # Gemmi requires float32 arrays
    mtz_out.set_data(data.astype(np.float32))

    # Write output file
    mtz_out.write_to_file(str(output_mtz))


def extract_columns_batch(
    input_mtz: str,
    output_specs: List[Tuple[str, List[Tuple[str, str, str]]]]
) -> None:
    """
    Extract multiple mini-MTZ files from a single input MTZ.

    This is useful when you need to create several output files with
    different column combinations from the same source MTZ.

    Args:
        input_mtz: Path to input MTZ file
        output_specs: List of (output_path, column_spec) tuples
            where column_spec is a list of (source_col, target_col, mtz_type)

    Example:
        >>> extract_columns_batch(
        ...     'servalcat_output.mtz',
        ...     [
        ...         ('fmean.mtz', [('F', 'F', 'F'), ('SIGF', 'SIGF', 'Q')]),
        ...         ('imean.mtz', [('I', 'I', 'J'), ('SIGI', 'SIGI', 'Q')])
        ...     ]
        ... )
    """
    for output_path, column_spec in output_specs:
        extract_columns(input_mtz, output_path, column_spec)


def get_column_labels(mtz_file: str) -> List[str]:
    """
    Get list of column labels from an MTZ file.

    Args:
        mtz_file: Path to MTZ file

    Returns:
        List of column label strings

    Example:
        >>> labels = get_column_labels('data.mtz')
        >>> print(labels)
        ['H', 'K', 'L', 'F', 'SIGF', 'PHIB', 'FOM']
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    return [col.label for col in mtz.columns]


def has_columns(mtz_file: str, required_columns: List[str]) -> bool:
    """
    Check if an MTZ file contains all required columns.

    Args:
        mtz_file: Path to MTZ file
        required_columns: List of column labels to check for

    Returns:
        True if all required columns are present, False otherwise

    Example:
        >>> has_columns('data.mtz', ['F', 'SIGF'])
        True
        >>> has_columns('data.mtz', ['F', 'SIGF', 'MISSING'])
        False
    """
    labels = get_column_labels(mtz_file)
    return all(col in labels for col in required_columns)
