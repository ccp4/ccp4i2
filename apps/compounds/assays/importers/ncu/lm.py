"""
NCU Liver Microsome Stability (LM) Parser.

Parses ADME-NCU-LM-YYYYMMDD.xlsx files containing metabolic stability
data from human and mouse liver microsomes.
"""

from pathlib import Path
from typing import Any, Optional

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class LiverMicrosomeParser(NCUBaseParser):
    """
    Parser for NCU Liver Microsome Stability data.

    Extracts:
    - Table 1: Summary results (t1/2, CLint, CLH, ER)
    - Table 2: Time course data (% remaining at each time point)

    Output format:
    - One ParsedResult per compound per species
    - Time course data nested with +/- cofactors conditions
    """

    assay_type = "liver_microsome_stability"
    assay_code = "LM"
    protocol_slug = "ncu-lm"
    kpi_field = "clint_ul_min_mg"  # In vitro CLint (μL/min/mg)

    # Expected time points (minutes)
    TIME_POINTS = [0.5, 5, 15, 30, 60]

    def parse(self, filepath: Path) -> ParseResult:
        """Parse an NCU Liver Microsome Stability file."""
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

        # Parse Table 1 (Summary Results)
        summary_data = self._parse_summary_table(tables[0], errors)

        # Parse Table 2 (Time Course) if present
        time_course_data = {}
        if len(tables) >= 2:
            time_course_data = self._parse_time_course_table(tables[1], errors)

        # Build source metadata
        source = self.build_source_metadata(filepath)

        # Combine data into results
        for (compound_id, species), metrics in summary_data.items():
            is_control = self.is_control_compound(compound_id)

            # Build results dict
            result_data = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "species": species,
                **metrics,
                "source": source,
                "flags": [],
            }

            # Add time course if available
            tc_key = (compound_id, species)
            if tc_key in time_course_data:
                result_data["time_course"] = time_course_data[tc_key]

            # Add quality flags
            self._add_quality_flags(result_data)

            results.append(ParsedResult(
                compound_id=compound_id,
                species=species,
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

    def _parse_summary_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[tuple[str, str], dict[str, Any]]:
        """
        Parse Table 1: Summary Results.

        Returns:
            Dict mapping (compound_id, species) to metrics dict.
        """
        results = {}
        headers = table.get('headers', [])
        rows = table.get('data', [])

        # Find column indices
        # In NCU ADME files, compound ID is always the first column
        compound_idx = self.get_compound_column_index(headers)
        species_idx = self.get_column_index(headers, 'species')
        t12_idx = self.get_column_index(headers, 't1/2', 't½', 'half-life')
        clint_idx = self.get_column_index(headers, 'in vitro clint', 'clint (μl')
        scaled_clint_idx = self.get_column_index(headers, 'scale-up clint', 'scaled clint')
        clh_idx = self.get_column_index(headers, 'clh', 'hepatic cl')
        er_idx = self.get_column_index(headers, 'extraction ratio', 'er')

        # Forward-fill compound IDs
        self.forward_fill_compound_ids(rows, compound_idx)

        for row_data in rows:
            values = row_data['values']
            row_num = row_data.get('row_number')

            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            species = values[species_idx] if species_idx and species_idx < len(values) else 'Unknown'
            if not species:
                species = 'Unknown'

            # Parse metrics
            metrics = {}

            if t12_idx and t12_idx < len(values):
                val, flag = self.parse_value(values[t12_idx])
                metrics['t1_2_min'] = val
                if flag:
                    metrics.setdefault('_flags', []).append(f"t1_2:{flag}")

            if clint_idx and clint_idx < len(values):
                val, flag = self.parse_value(values[clint_idx])
                metrics['clint_ul_min_mg'] = val

            if scaled_clint_idx and scaled_clint_idx < len(values):
                val, flag = self.parse_value(values[scaled_clint_idx])
                metrics['clint_scaled_ml_min_kg'] = val

            if clh_idx and clh_idx < len(values):
                val, flag = self.parse_value(values[clh_idx])
                metrics['clh_predicted_ml_min_kg'] = val

            if er_idx and er_idx < len(values):
                val, flag = self.parse_value(values[er_idx])
                metrics['extraction_ratio'] = val

            results[(compound_id, species)] = metrics

        return results

    def _parse_time_course_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[tuple[str, str], dict[str, Any]]:
        """
        Parse Table 2: Time Course Data.

        Returns:
            Dict mapping (compound_id, species) to time course dict.
        """
        results = {}
        headers = table.get('headers', [])
        rows = table.get('data', [])

        # Find column indices
        # In NCU ADME files, compound ID is always the first column
        compound_idx = self.get_compound_column_index(headers)
        species_idx = self.get_column_index(headers, 'species')
        format_idx = self.get_column_index(headers, 'assay format', 'format')

        # Find time point columns
        time_indices, time_values = self.extract_time_points(headers)

        if not time_indices:
            errors.append(ValidationError(
                message="No time point columns found in Table 2",
                severity="warning"
            ))
            return results

        # Forward-fill compound IDs
        self.forward_fill_compound_ids(rows, compound_idx)

        for row_data in rows:
            values = row_data['values']

            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            species = values[species_idx] if species_idx and species_idx < len(values) else 'Unknown'
            if not species:
                species = 'Unknown'

            assay_format = values[format_idx] if format_idx and format_idx < len(values) else ''

            # Parse time course values
            remaining_pct = []
            for idx in time_indices:
                if idx < len(values):
                    val, flag = self.parse_value(values[idx])
                    remaining_pct.append(val)
                else:
                    remaining_pct.append(None)

            # Determine if +Cofactors or -Cofactors
            is_with_cofactors = '+' in str(assay_format).lower() or 'cofactor' in str(assay_format).lower()
            condition_key = 'with_cofactors' if is_with_cofactors or '+' in str(assay_format) else 'without_cofactors'

            # Get or create time course entry
            key = (compound_id, species)
            if key not in results:
                results[key] = {
                    'time_points_min': time_values,
                    'with_cofactors': {'remaining_pct': []},
                    'without_cofactors': {'remaining_pct': []},
                }

            results[key][condition_key]['remaining_pct'] = remaining_pct

        return results

    def _add_quality_flags(self, result_data: dict) -> None:
        """Add quality flags based on result values."""
        flags = result_data.get('flags', [])

        t12 = result_data.get('t1_2_min')
        if t12 is not None:
            if t12 < 10:
                flags.append('unstable')
            elif t12 > 120:
                flags.append('stable')

        clint = result_data.get('clint_ul_min_mg')
        if clint is not None and clint > 100:
            flags.append('high_clearance')

        # Add any flags from parsing
        if '_flags' in result_data:
            flags.extend(result_data.pop('_flags'))

        result_data['flags'] = list(set(flags))  # Deduplicate
