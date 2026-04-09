"""
Plate Layout Utilities

Helper functions and validation for protocol plate layouts.
Handles well coordinate parsing, layout validation, and data extraction.
"""

import re
from typing import Optional


# Plate format dimensions
PLATE_FORMATS = {
    24: {'rows': 4, 'cols': 6, 'row_labels': 'ABCD'},
    96: {'rows': 8, 'cols': 12, 'row_labels': 'ABCDEFGH'},
    384: {'rows': 16, 'cols': 24, 'row_labels': 'ABCDEFGHIJKLMNOP'},
    1536: {'rows': 32, 'cols': 48, 'row_labels': 'ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEF'},
}

# Replicate pattern options
REPLICATE_PATTERNS = [
    'adjacent_rows',      # A1,B1 are reps of compound 1; C1,D1 are reps of compound 2
    'adjacent_columns',   # A1,A2 are reps of compound 1
    'grouped_rows',       # Rows A-H are rep 1, rows I-P are rep 2
    'interleaved_rows',   # Odd rows rep 1, even rows rep 2
]

# Compound source types
COMPOUND_SOURCE_TYPES = [
    'row_order',          # Compounds ordered by row in data file
    'column_header',      # Compound IDs in column headers
    'row_header',         # Compound IDs in row headers
    'adjacent_column',    # Compound names in column immediately after data region
    'plate_map_file',     # Separate plate map file
    'explicit_wells',     # Explicit well-to-compound mapping in layout
]


def parse_well(well: str) -> tuple[str, int]:
    """
    Parse a well coordinate string into row letter and column number.

    Args:
        well: Well coordinate like 'A1', 'B12', 'P24'

    Returns:
        Tuple of (row_letter, column_number)

    Raises:
        ValueError: If well format is invalid
    """
    match = re.match(r'^([A-Za-z]+)(\d+)$', well.strip())
    if not match:
        raise ValueError(f"Invalid well format: {well}")
    return match.group(1).upper(), int(match.group(2))


def well_to_indices(well: str, plate_format: int = 384) -> tuple[int, int]:
    """
    Convert well coordinate to 0-based row/column indices.

    Args:
        well: Well coordinate like 'A1', 'B12'
        plate_format: Plate format (96, 384, 1536)

    Returns:
        Tuple of (row_index, column_index), 0-based
    """
    row_letter, col_num = parse_well(well)
    format_info = PLATE_FORMATS.get(plate_format, PLATE_FORMATS[384])

    row_index = format_info['row_labels'].index(row_letter)
    col_index = col_num - 1  # Convert 1-based to 0-based

    return row_index, col_index


def indices_to_well(row_index: int, col_index: int, plate_format: int = 384) -> str:
    """
    Convert 0-based row/column indices to well coordinate.

    Args:
        row_index: 0-based row index
        col_index: 0-based column index
        plate_format: Plate format (96, 384, 1536)

    Returns:
        Well coordinate string like 'A1', 'B12'
    """
    format_info = PLATE_FORMATS.get(plate_format, PLATE_FORMATS[384])
    row_letter = format_info['row_labels'][row_index]
    col_num = col_index + 1  # Convert 0-based to 1-based
    return f"{row_letter}{col_num}"


def expand_well_range(
    columns: list[int],
    rows: list[str],
    plate_format: int = 384
) -> list[str]:
    """
    Expand column/row specification into list of well coordinates.

    Args:
        columns: List of column numbers (1-based)
        rows: List of row letters
        plate_format: Plate format

    Returns:
        List of well coordinates
    """
    wells = []
    for row in rows:
        for col in columns:
            wells.append(f"{row.upper()}{col}")
    return wells


def get_control_wells(layout: dict) -> dict[str, list[str]]:
    """
    Extract control well positions from plate layout.

    Args:
        layout: Plate layout dictionary

    Returns:
        Dict mapping control type ('max', 'min') to list of well coordinates
    """
    controls = layout.get('controls', {})
    plate_format = layout.get('plate_format', 384)
    result = {}

    for control_type, spec in controls.items():
        if 'wells' in spec:
            # Explicit well list
            result[control_type] = spec['wells']
        elif 'columns' in spec and 'rows' in spec:
            # Column/row expansion
            result[control_type] = expand_well_range(
                spec['columns'],
                spec['rows'],
                plate_format
            )
        else:
            result[control_type] = []

    return result


def get_sample_wells(layout: dict) -> list[tuple[int, int]]:
    """
    Get the sample region as a list of (row_index, col_index) tuples.

    Args:
        layout: Plate layout dictionary

    Returns:
        List of (row, col) index tuples for sample wells
    """
    plate_format = layout.get('plate_format', 384)
    format_info = PLATE_FORMATS.get(plate_format, PLATE_FORMATS[384])

    sample_region = layout.get('sample_region', {})

    start_col = sample_region.get('start_column', 3) - 1  # Convert to 0-based
    end_col = sample_region.get('end_column', format_info['cols'] - 2)
    start_row = format_info['row_labels'].index(
        sample_region.get('start_row', 'A').upper()
    )
    end_row = format_info['row_labels'].index(
        sample_region.get('end_row', format_info['row_labels'][-1]).upper()
    )

    wells = []
    for row in range(start_row, end_row + 1):
        for col in range(start_col, end_col):
            wells.append((row, col))

    return wells


def validate_plate_layout(layout: dict) -> list[str]:
    """
    Validate a plate layout configuration.

    Args:
        layout: Plate layout dictionary

    Returns:
        List of validation error messages (empty if valid)
    """
    errors = []

    if not layout:
        return errors  # Empty layout is valid (not configured)

    # Check plate format
    plate_format = layout.get('plate_format', 384)
    if plate_format not in PLATE_FORMATS:
        errors.append(f"Invalid plate format: {plate_format}. Must be one of {list(PLATE_FORMATS.keys())}")
        return errors  # Can't validate further without valid format

    format_info = PLATE_FORMATS[plate_format]

    # Validate controls
    controls = layout.get('controls', {})
    for control_type, spec in controls.items():
        if 'columns' in spec:
            for col in spec['columns']:
                if col < 1 or col > format_info['cols']:
                    errors.append(
                        f"Control '{control_type}' column {col} out of range "
                        f"(1-{format_info['cols']})"
                    )
        if 'rows' in spec:
            for row in spec['rows']:
                if row.upper() not in format_info['row_labels']:
                    errors.append(
                        f"Control '{control_type}' row {row} not valid for "
                        f"{plate_format}-well plate"
                    )

    # Validate sample region
    sample_region = layout.get('sample_region', {})
    if sample_region:
        start_col = sample_region.get('start_column', 1)
        end_col = sample_region.get('end_column', format_info['cols'])
        if start_col < 1 or start_col > format_info['cols']:
            errors.append(f"Sample region start_column {start_col} out of range")
        if end_col < 1 or end_col > format_info['cols']:
            errors.append(f"Sample region end_column {end_col} out of range")
        if start_col >= end_col:
            errors.append("Sample region start_column must be less than end_column")

        start_row = sample_region.get('start_row', 'A').upper()
        end_row = sample_region.get('end_row', format_info['row_labels'][-1]).upper()
        if start_row not in format_info['row_labels']:
            errors.append(f"Sample region start_row '{start_row}' not valid")
        if end_row not in format_info['row_labels']:
            errors.append(f"Sample region end_row '{end_row}' not valid")

    # Validate replicate pattern
    replicate = layout.get('replicate', {})
    if replicate:
        pattern = replicate.get('pattern')
        if pattern and pattern not in REPLICATE_PATTERNS:
            errors.append(
                f"Invalid replicate pattern: {pattern}. "
                f"Must be one of {REPLICATE_PATTERNS}"
            )
        count = replicate.get('count', 1)
        if count < 1 or count > 4:
            errors.append(f"Replicate count {count} out of range (1-4)")

    # Validate compound source
    compound_source = layout.get('compound_source', {})
    if compound_source:
        source_type = compound_source.get('type')
        if source_type and source_type not in COMPOUND_SOURCE_TYPES:
            errors.append(
                f"Invalid compound source type: {source_type}. "
                f"Must be one of {COMPOUND_SOURCE_TYPES}"
            )

    return errors


def apply_offset(layout: dict, row_offset: int = 0, col_offset: int = 0) -> dict:
    """
    Apply row/column offset to a plate layout.

    Used at import time when the plate pattern is shifted from the default.

    Args:
        layout: Original plate layout
        row_offset: Number of rows to shift (positive = down)
        col_offset: Number of columns to shift (positive = right)

    Returns:
        New layout dict with offset applied
    """
    if row_offset == 0 and col_offset == 0:
        return layout

    plate_format = layout.get('plate_format', 384)
    format_info = PLATE_FORMATS.get(plate_format, PLATE_FORMATS[384])

    result = layout.copy()

    # Offset controls
    if 'controls' in result:
        result['controls'] = {}
        for control_type, spec in layout['controls'].items():
            new_spec = spec.copy()
            if 'columns' in spec:
                new_spec['columns'] = [c + col_offset for c in spec['columns']]
            if 'rows' in spec:
                new_rows = []
                for row in spec['rows']:
                    idx = format_info['row_labels'].index(row.upper())
                    new_idx = idx + row_offset
                    if 0 <= new_idx < len(format_info['row_labels']):
                        new_rows.append(format_info['row_labels'][new_idx])
                new_spec['rows'] = new_rows
            result['controls'][control_type] = new_spec

    # Offset sample region
    if 'sample_region' in result:
        sr = result['sample_region'] = layout['sample_region'].copy()
        if 'start_column' in sr:
            sr['start_column'] = sr['start_column'] + col_offset
        if 'end_column' in sr:
            sr['end_column'] = sr['end_column'] + col_offset
        if 'start_row' in sr:
            idx = format_info['row_labels'].index(sr['start_row'].upper())
            new_idx = idx + row_offset
            if 0 <= new_idx < len(format_info['row_labels']):
                sr['start_row'] = format_info['row_labels'][new_idx]
        if 'end_row' in sr:
            idx = format_info['row_labels'].index(sr['end_row'].upper())
            new_idx = idx + row_offset
            if 0 <= new_idx < len(format_info['row_labels']):
                sr['end_row'] = format_info['row_labels'][new_idx]

    return result


# Default plate layouts for common configurations
DEFAULT_LAYOUTS = {
    '24_standard': {
        'plate_format': 24,
        'controls': {
            'max': {'columns': [1], 'rows': ['A']},
            'min': {'columns': [6], 'rows': ['A']},
        },
        'sample_region': {
            'start_column': 2,
            'end_column': 5,
            'start_row': 'A',
            'end_row': 'D',
        },
        'dilution': {
            'direction': 'horizontal',
            'num_concentrations': 4,
        },
        'replicate': {
            'count': 2,
            'pattern': 'adjacent_rows',
        },
        'compound_source': {
            'type': 'row_order',
        },
    },
    '96_standard': {
        'plate_format': 96,
        'controls': {
            'max': {'columns': [1], 'rows': ['A', 'B', 'C', 'D']},
            'min': {'columns': [12], 'rows': ['A', 'B', 'C', 'D']},
        },
        'sample_region': {
            'start_column': 2,
            'end_column': 11,
            'start_row': 'A',
            'end_row': 'H',
        },
        'dilution': {
            'direction': 'horizontal',
            'num_concentrations': 10,
        },
        'replicate': {
            'count': 2,
            'pattern': 'adjacent_rows',
        },
        'compound_source': {
            'type': 'row_order',
        },
    },
    '384_standard': {
        'plate_format': 384,
        'controls': {
            'max': {'columns': [1, 2], 'rows': ['A', 'B']},
            'min': {'columns': [23, 24], 'rows': ['A', 'B']},
        },
        'sample_region': {
            'start_column': 3,
            'end_column': 22,
            'start_row': 'A',
            'end_row': 'P',
        },
        'dilution': {
            'direction': 'horizontal',
            'num_concentrations': 10,
        },
        'replicate': {
            'count': 2,
            'pattern': 'adjacent_rows',
        },
        'compound_source': {
            'type': 'row_order',
        },
    },
    '1536_standard': {
        'plate_format': 1536,
        'controls': {
            'max': {'columns': [1, 2], 'rows': ['A', 'B', 'C', 'D']},
            'min': {'columns': [47, 48], 'rows': ['A', 'B', 'C', 'D']},
        },
        'sample_region': {
            'start_column': 3,
            'end_column': 46,
            'start_row': 'A',
            'end_row': 'AF',  # Row 32
        },
        'dilution': {
            'direction': 'horizontal',
            'num_concentrations': 22,
        },
        'replicate': {
            'count': 2,
            'pattern': 'adjacent_rows',
        },
        'compound_source': {
            'type': 'row_order',
        },
    },
}
