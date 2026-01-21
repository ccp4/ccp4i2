"""
NCU GSH Stability Parser.

Parses ADME-NCU-GSH Stability-YYYYMMDD.xlsx files containing
glutathione stability data for reactive metabolite screening.
"""

from pathlib import Path
from typing import Any

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class GSHStabilityParser(NCUBaseParser):
    """
    Parser for NCU GSH Stability data.

    Extracts:
    - Table 1: GSH Stability results (+GSH and -GSH conditions)

    Output format:
    - One ParsedResult per compound
    - Both +GSH and -GSH conditions nested in time_course
    """

    assay_type = "gsh_stability"
    assay_code = "GSH"
    protocol_slug = "ncu-gsh"
    kpi_field = "t1_2_min"  # t1/2 (min)

    # Expected time points (minutes)
    TIME_POINTS = [0, 15, 30, 60, 120]

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        """Match 'gsh' or 'gsh-stability' or 'gsh stability'."""
        return 'gsh' in assay_type_from_file.lower()

    def parse(self, filepath: Path) -> ParseResult:
        """Parse an NCU GSH Stability file."""
        errors: list[ValidationError] = []
        results: list[ParsedResult] = []

        try:
            sheets = self.read_excel(filepath)
        except Exception as e:
            return ParseResult(errors=[ValidationError(
                message=f"Failed to read Excel file: {e}",
                severity="error"
            )])

        # Find Summary sheet with diagnostic info
        summary, summary_sheet_name, all_sheet_names = self.find_summary_sheet_with_name(sheets)

        if summary is None:
            return ParseResult(
                errors=[ValidationError(
                    message=f"Summary sheet not found. Available sheets: {all_sheet_names}",
                    severity="error"
                )],
                metadata={
                    "sheets_found": all_sheet_names,
                }
            )

        # Parse tables from Summary sheet
        tables = self.parse_summary_tables(summary, errors)

        if len(tables) < 1:
            errors.append(ValidationError(
                message=f"No data tables found in sheet '{summary_sheet_name}'. "
                        f"Expected 'Table 1.' header row.",
                severity="error"
            ))
            return ParseResult(
                errors=errors,
                metadata={
                    "sheets_found": all_sheet_names,
                    "summary_sheet_used": summary_sheet_name,
                }
            )

        # Parse Table 1 (GSH Stability Results)
        stability_data = self._parse_gsh_table(tables[0], errors)

        # Build source metadata
        source = self.build_source_metadata(filepath)

        # Build results - combine +GSH and -GSH into single result per compound
        compound_data = self._combine_conditions(stability_data)

        for compound_id, metrics in compound_data.items():
            is_control = self.is_control_compound(compound_id)

            # Build results dict
            result_data = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                **metrics,
                "source": source,
                "flags": [],
            }

            # Determine if GSH reactive
            self._assess_gsh_reactivity(result_data)

            # Add quality flags
            self._add_quality_flags(result_data)

            results.append(ParsedResult(
                compound_id=compound_id,
                results=result_data,
                is_control=is_control,
            ))

        # Validate results
        validation_errors = self.validate(results)
        errors.extend(validation_errors)

        return ParseResult(
            results=results,
            errors=errors,
            metadata={
                "vendor": self.vendor,
                "assay_type": self.assay_type,
                "file": filepath.name,
                "sheets_found": all_sheet_names,
                "summary_sheet_used": summary_sheet_name,
                "tables_found": len(tables),
            }
        )

    def _parse_gsh_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[tuple[str, str], dict[str, Any]]:
        """
        Parse Table 1: GSH Stability Results.

        Returns:
            Dict mapping (compound_id, assay_format) to metrics dict.
        """
        results = {}
        headers = table.get('headers', [])
        rows = table.get('data', [])

        # Find column indices
        # In NCU ADME files, compound ID is always the first column
        compound_idx = self.get_compound_column_index(headers)
        format_idx = self.get_column_index(headers, 'assay format', 'format', 'condition')
        t12_idx = self.get_column_index(headers, 't1/2', 't½', 'half-life')

        # Find time point columns
        time_indices, time_values = self.extract_time_points(headers)

        # Forward-fill compound IDs
        self.forward_fill_compound_ids(rows, compound_idx)

        for row_data in rows:
            values = row_data['values']

            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            assay_format = values[format_idx] if format_idx and format_idx < len(values) else ''

            # Parse t1/2
            metrics = {}
            flags = []

            if t12_idx and t12_idx < len(values):
                val, flag = self.parse_value(values[t12_idx])
                metrics['t1_2_min'] = val
                if flag:
                    flags.append(flag)

            # Parse time course values
            remaining_pct = []
            for idx in time_indices:
                if idx < len(values):
                    val, flag = self.parse_value(values[idx])
                    remaining_pct.append(val)
                else:
                    remaining_pct.append(None)

            metrics['remaining_pct'] = remaining_pct
            metrics['time_points_min'] = time_values if time_values else self.TIME_POINTS

            if flags:
                metrics['_flags'] = flags

            results[(compound_id, assay_format)] = metrics

        return results

    def _combine_conditions(
        self,
        raw_data: dict[tuple[str, str], dict[str, Any]],
    ) -> dict[str, dict[str, Any]]:
        """
        Combine +GSH and -GSH conditions into single result per compound.

        Returns:
            Dict mapping compound_id to combined metrics.
        """
        combined = {}

        for (compound_id, assay_format), metrics in raw_data.items():
            if compound_id not in combined:
                combined[compound_id] = {
                    'time_course': {
                        'time_points_min': metrics.get('time_points_min', self.TIME_POINTS),
                        'with_gsh': {'remaining_pct': []},
                        'without_gsh': {'remaining_pct': []},
                    },
                    't1_2_min': None,
                    't1_2_min_no_gsh': None,
                }

            # Determine if +GSH or -GSH
            is_with_gsh = '+gsh' in str(assay_format).lower() or assay_format.startswith('+')

            if is_with_gsh:
                combined[compound_id]['time_course']['with_gsh']['remaining_pct'] = metrics.get('remaining_pct', [])
                combined[compound_id]['t1_2_min'] = metrics.get('t1_2_min')
            else:
                combined[compound_id]['time_course']['without_gsh']['remaining_pct'] = metrics.get('remaining_pct', [])
                combined[compound_id]['t1_2_min_no_gsh'] = metrics.get('t1_2_min')

            # Collect flags
            if '_flags' in metrics:
                combined[compound_id].setdefault('_flags', []).extend(metrics['_flags'])

        return combined

    def _assess_gsh_reactivity(self, result_data: dict) -> None:
        """
        Assess if compound is GSH-reactive based on t1/2 comparison.

        GSH-reactive if t1/2 with GSH is significantly lower than without GSH.
        """
        t12_gsh = result_data.get('t1_2_min')
        t12_no_gsh = result_data.get('t1_2_min_no_gsh')

        if t12_gsh is not None and t12_no_gsh is not None:
            # If t1/2 with GSH is less than half of t1/2 without GSH
            if t12_no_gsh > 0 and t12_gsh < t12_no_gsh * 0.5:
                result_data['gsh_reactive'] = True
            else:
                result_data['gsh_reactive'] = False
        elif t12_gsh is not None and t12_no_gsh is None:
            # If we only have +GSH data and t1/2 is short, likely reactive
            result_data['gsh_reactive'] = t12_gsh < 30  # Threshold
        else:
            result_data['gsh_reactive'] = None

    def _add_quality_flags(self, result_data: dict) -> None:
        """Add quality flags based on result values."""
        flags = result_data.get('flags', [])

        t12 = result_data.get('t1_2_min')
        if t12 is not None:
            if t12 < 10:
                flags.append('unstable')
            elif t12 > 120:
                flags.append('stable')

        if result_data.get('gsh_reactive'):
            flags.append('gsh_reactive')

        # Add any flags from parsing
        if '_flags' in result_data:
            flags.extend(result_data.pop('_flags'))

        result_data['flags'] = list(set(flags))  # Deduplicate
