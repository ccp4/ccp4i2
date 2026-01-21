"""
Base utilities for NCU ADME parsers.

This module provides common functionality shared by all NCU parsers,
including Excel file reading and Summary sheet parsing.
"""

import re
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..base import ADMEParser, ParseResult, ParsedResult, ValidationError


class NCUBaseParser(ADMEParser):
    """
    Base class for NCU ADME parsers.

    Provides common functionality for:
    - File detection by naming convention
    - Reading Excel Summary sheets
    - Parsing multi-table Summary sheets
    - Handling NCU-specific special values
    """

    vendor = "NCU"

    # NCU file naming pattern: ADME-NCU-{AssayType}-{YYYYMMDD}.xlsx
    FILE_PATTERN = r'^ADME-NCU-(.+)-(\d{8})\.xlsx$'

    # Special values that need conversion
    SPECIAL_VALUES = {
        'BLOD': (None, 'below_lod'),      # Below Limit of Detection
        '-': (None, None),                 # Not measured
        '∞': (None, 'stable'),             # Infinite (stable)
        'inf': (None, 'stable'),           # Alternative infinite
        'N/A': (None, None),               # Not applicable
        'NA': (None, None),                # Not applicable
    }

    def detect(self, filepath: Path) -> bool:
        """
        Detect if this parser can handle the file.

        Checks:
        1. File matches NCU naming convention
        2. Assay type in filename matches this parser's assay_code
        """
        if not filepath.exists():
            return False

        if filepath.suffix.lower() != '.xlsx':
            return False

        match = re.match(self.FILE_PATTERN, filepath.name, re.IGNORECASE)
        if not match:
            return False

        # Check if assay type matches
        assay_type_from_file = match.group(1).lower().replace(' ', '-')
        return self._matches_assay_type(assay_type_from_file)

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        """
        Check if the assay type from filename matches this parser.

        Override in subclasses if needed for complex matching.
        """
        return assay_type_from_file.lower() == self.assay_code.lower()

    def read_excel(self, filepath: Path) -> dict[str, pd.DataFrame]:
        """
        Read all sheets from an Excel file.

        Args:
            filepath: Path to the Excel file.

        Returns:
            Dict mapping sheet names to DataFrames.
        """
        return pd.read_excel(filepath, sheet_name=None, header=None)

    def find_summary_sheet(self, sheets: dict[str, pd.DataFrame]) -> Optional[pd.DataFrame]:
        """
        Find the Summary sheet in the workbook.

        NCU files typically have:
        - First sheet: Metadata/cover sheet (often named with project info)
        - Second sheet: "Summary" with actual data tables

        We prefer exact "Summary" match, then fall back to partial match,
        avoiding metadata-looking sheets.

        Args:
            sheets: Dict of sheet names to DataFrames.

        Returns:
            The Summary sheet DataFrame or None if not found.
        """
        import logging
        logger = logging.getLogger(__name__)

        sheet_names = list(sheets.keys())
        logger.debug(f"NCU Parser: Available sheets: {sheet_names}")

        # First, try exact match (case-insensitive)
        for name, df in sheets.items():
            if name.lower().strip() == 'summary':
                logger.debug(f"NCU Parser: Found exact Summary sheet: '{name}'")
                return df

        # Second, try "Data Summary" or similar common patterns
        for name, df in sheets.items():
            lower_name = name.lower().strip()
            if lower_name in ('data summary', 'summary data', 'results summary'):
                logger.debug(f"NCU Parser: Found Summary sheet variant: '{name}'")
                return df

        # Third, try partial match but skip obvious metadata sheets
        metadata_keywords = ['info', 'metadata', 'cover', 'front', 'introduction', 'about']
        for name, df in sheets.items():
            lower_name = name.lower()
            if 'summary' in lower_name:
                # Skip if it looks like a metadata/info sheet
                if any(kw in lower_name for kw in metadata_keywords):
                    logger.debug(f"NCU Parser: Skipping metadata-like sheet: '{name}'")
                    continue
                logger.debug(f"NCU Parser: Found Summary sheet by partial match: '{name}'")
                return df

        # Last resort: if there are exactly 2 sheets, the second is usually Summary
        if len(sheets) == 2:
            second_sheet_name = sheet_names[1]
            logger.info(f"NCU Parser: Falling back to second sheet: '{second_sheet_name}'")
            return sheets[second_sheet_name]

        logger.warning(f"NCU Parser: Could not find Summary sheet in: {sheet_names}")
        return None

    def find_summary_sheet_with_name(self, sheets: dict[str, pd.DataFrame]) -> tuple[Optional[pd.DataFrame], Optional[str], list[str]]:
        """
        Find the Summary sheet and return both the DataFrame and its name.

        Returns:
            Tuple of (DataFrame or None, sheet_name or None, list of all sheet names)
        """
        sheet_names = list(sheets.keys())

        # First, try exact match (case-insensitive)
        for name, df in sheets.items():
            if name.lower().strip() == 'summary':
                return df, name, sheet_names

        # Second, try "Data Summary" or similar common patterns
        for name, df in sheets.items():
            lower_name = name.lower().strip()
            if lower_name in ('data summary', 'summary data', 'results summary'):
                return df, name, sheet_names

        # Third, try partial match but skip obvious metadata sheets
        metadata_keywords = ['info', 'metadata', 'cover', 'front', 'introduction', 'about']
        for name, df in sheets.items():
            lower_name = name.lower()
            if 'summary' in lower_name:
                if any(kw in lower_name for kw in metadata_keywords):
                    continue
                return df, name, sheet_names

        # Last resort: if there are exactly 2 sheets, the second is usually Summary
        if len(sheets) == 2:
            second_sheet_name = sheet_names[1]
            return sheets[second_sheet_name], second_sheet_name, sheet_names

        return None, None, sheet_names

    def parse_summary_tables(
        self,
        df: pd.DataFrame,
        errors: list[ValidationError],
    ) -> list[dict[str, Any]]:
        """
        Parse tables from an NCU Summary sheet.

        NCU Summary sheets have a specific structure:
        - Row 0: "Data Summary"
        - Row 1: "Table 1. {description}"
        - Row 2: Column headers
        - Row 3+: Data rows
        - Tables separated by blank rows or "Note:" rows

        Args:
            df: Summary sheet DataFrame.
            errors: List to append validation errors to.

        Returns:
            List of dicts, each representing a table with 'headers' and 'data' keys.
        """
        tables = []
        current_table = None
        in_data_section = False

        for idx, row in df.iterrows():
            row_values = [str(v).strip() if pd.notna(v) else '' for v in row.values]
            first_cell = row_values[0] if row_values else ''

            # Detect table start
            if first_cell.lower().startswith('table'):
                if current_table and current_table.get('data'):
                    tables.append(current_table)
                current_table = {
                    'title': first_cell,
                    'headers': [],
                    'data': [],
                    'start_row': idx,
                }
                in_data_section = False
                continue

            # Skip header row and capture column names
            if current_table and not current_table['headers']:
                # This should be the header row
                headers = [v for v in row_values if v]
                if headers:
                    current_table['headers'] = row_values
                    in_data_section = True
                continue

            # Check for end of table
            if self._is_table_end(first_cell):
                if current_table and current_table.get('data'):
                    tables.append(current_table)
                current_table = None
                in_data_section = False
                continue

            # Collect data rows
            if current_table and in_data_section:
                # Skip completely empty rows
                if not any(v for v in row_values):
                    continue
                current_table['data'].append({
                    'values': row_values,
                    'row_number': idx,
                })

        # Don't forget the last table
        if current_table and current_table.get('data'):
            tables.append(current_table)

        return tables

    def _is_table_end(self, cell_value: str) -> bool:
        """Check if a cell value indicates end of table."""
        if not cell_value:
            return False
        lower = cell_value.lower()
        return (
            lower.startswith('note:') or
            lower.startswith('table') or
            lower == 'data summary'
        )

    def parse_value(
        self,
        value: str,
        expected_type: type = float,
    ) -> tuple[Any, Optional[str]]:
        """
        Parse a cell value, handling NCU special values.

        Args:
            value: Cell value as string.
            expected_type: Expected type (float, int, str).

        Returns:
            Tuple of (parsed_value, flag) where flag indicates special handling.
        """
        if not value or pd.isna(value):
            return None, None

        value_str = str(value).strip()

        # Check for special values
        if value_str in self.SPECIAL_VALUES:
            return self.SPECIAL_VALUES[value_str]

        # Check for ">" pattern (exceeds measurement window)
        if value_str.startswith('>'):
            try:
                num_part = float(value_str[1:].strip())
                return num_part, 'exceeds_window'
            except ValueError:
                return None, 'exceeds_window'

        # Check for "<" pattern (below quantification)
        if value_str.startswith('<'):
            try:
                num_part = float(value_str[1:].strip())
                return num_part, 'below_loq'
            except ValueError:
                return None, 'below_loq'

        # Parse as expected type
        try:
            if expected_type == float:
                # Handle percentage strings
                if value_str.endswith('%'):
                    return float(value_str[:-1]), None
                return float(value_str), None
            elif expected_type == int:
                return int(float(value_str)), None
            else:
                return value_str, None
        except (ValueError, TypeError):
            return value_str, None

    def forward_fill_compound_ids(
        self,
        rows: list[dict],
        compound_col_idx: int = 0,
    ) -> None:
        """
        Forward-fill compound IDs in multi-row data.

        NCU files often have compound IDs only on the first row
        when multiple species/conditions are tested.

        Args:
            rows: List of row dicts (modified in place).
            compound_col_idx: Index of compound ID column.
        """
        last_compound = None
        for row in rows:
            values = row.get('values', [])
            if compound_col_idx < len(values):
                if values[compound_col_idx]:
                    last_compound = values[compound_col_idx]
                elif last_compound:
                    values[compound_col_idx] = last_compound

    def get_column_index(
        self,
        headers: list[str],
        *possible_names: str,
    ) -> Optional[int]:
        """
        Find column index by name (case-insensitive, partial match).

        Args:
            headers: List of header strings.
            possible_names: Possible column names to match.

        Returns:
            Column index or None if not found.
        """
        for i, header in enumerate(headers):
            if not header:
                continue
            header_lower = header.lower().strip()
            for name in possible_names:
                if name.lower() in header_lower:
                    return i
        return None

    def get_compound_column_index(
        self,
        headers: list[str],
    ) -> int:
        """
        Get the compound ID column index.

        In NCU ADME Summary tables, the compound ID is ALWAYS in the first
        column. This method attempts name-based matching first but falls
        back to column 0 if no match is found.

        This differs from Table-of-Values imports where the user explicitly
        specifies which column contains compound identifiers.

        Args:
            headers: List of header strings.

        Returns:
            Column index (defaults to 0 if not found by name).
        """
        # Try name-based matching first
        idx = self.get_column_index(headers, 'compound id', 'compound', 'sample', 'name')
        if idx is not None:
            return idx

        # Fall back to first column - in NCU ADME files, compound is always first
        return 0

    def extract_time_points(self, headers: list[str]) -> tuple[list[int], list[float]]:
        """
        Extract time point columns from headers.

        Looks for patterns like "0 min", "5 min", "0.5 min", etc.

        Args:
            headers: List of header strings.

        Returns:
            Tuple of (column_indices, time_values_in_minutes).
        """
        indices = []
        times = []
        time_pattern = re.compile(r'^([\d.]+)\s*min', re.IGNORECASE)

        for i, header in enumerate(headers):
            if not header:
                continue
            match = time_pattern.match(header.strip())
            if match:
                try:
                    time_val = float(match.group(1))
                    indices.append(i)
                    times.append(time_val)
                except ValueError:
                    pass

        return indices, times
