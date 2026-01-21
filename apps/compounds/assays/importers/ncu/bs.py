"""
NCU Blood/Serum Stability (BS) Parser.

Parses ADME-NCU-BS-YYYYMMDD.xlsx files containing stability data
from blood or serum incubations.
"""

from pathlib import Path
from typing import Any

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class BloodSerumStabilityParser(NCUBaseParser):
    """
    Parser for NCU Blood/Serum Stability data.

    Extracts:
    - Table 1: Stability results (% remaining at time points, t1/2)

    Output format:
    - One ParsedResult per compound per species
    """

    assay_type = "blood_serum_stability"
    assay_code = "BS"
    protocol_slug = "ncu-bs"
    kpi_field = "t1_2_min"  # t1/2 (min)

    # Expected time points (minutes)
    TIME_POINTS = [0, 15, 30, 60, 120]

    def parse(self, filepath: Path) -> ParseResult:
        """Parse an NCU Blood/Serum Stability file."""
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

        # Parse Table 1 (Stability Results)
        stability_data = self._parse_stability_table(tables[0], errors)

        # Build source metadata
        source = self.build_source_metadata(filepath)

        # Build results
        for (compound_id, species), metrics in stability_data.items():
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

    def _parse_stability_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[tuple[str, str], dict[str, Any]]:
        """
        Parse Table 1: Stability Results.

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

        # Find time point columns
        time_indices, time_values = self.extract_time_points(headers)

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
                    if flag:
                        flags.append(f"t{time_values[time_indices.index(idx)]}:{flag}")
                else:
                    remaining_pct.append(None)

            metrics['time_course'] = {
                'time_points_min': time_values if time_values else self.TIME_POINTS,
                'remaining_pct': remaining_pct,
            }

            if flags:
                metrics['_flags'] = flags

            results[(compound_id, species)] = metrics

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

        # Add any flags from parsing
        if '_flags' in result_data:
            flags.extend(result_data.pop('_flags'))

        result_data['flags'] = list(set(flags))  # Deduplicate
