"""
NCU Caco-2 Permeability Parser.

Parses ADME-NCU-Caco-2 Permeability-YYYYMMDD.xlsx files containing
intestinal permeability data using Caco-2 cell monolayers.
"""

from pathlib import Path
from typing import Any, Optional

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class Caco2PermeabilityParser(NCUBaseParser):
    """
    Parser for NCU Caco-2 Permeability data.

    Extracts:
    - Table 1: Permeability results (Papp A-B, Papp B-A, Efflux Ratio, Recovery)
    - Table 2: Monolayer integrity (TEER, LY Leakage)

    Output format:
    - One ParsedResult per compound
    """

    assay_type = "caco2_permeability"
    assay_code = "Caco-2"
    protocol_slug = "ncu-caco2"
    kpi_field = "efflux_ratio"  # Efflux Ratio

    # Permeability classification thresholds (10^-6 cm/s)
    PAPP_HIGH_THRESHOLD = 10.0
    PAPP_MEDIUM_THRESHOLD = 1.0

    # Efflux ratio threshold
    EFFLUX_THRESHOLD = 2.0
    HIGH_EFFLUX_THRESHOLD = 10.0

    # Recovery threshold
    LOW_RECOVERY_THRESHOLD = 70.0

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        """Match 'caco-2' or 'caco2' or 'caco-2-permeability'."""
        return 'caco' in assay_type_from_file.lower()

    def parse(self, filepath: Path) -> ParseResult:
        """Parse an NCU Caco-2 Permeability file."""
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

        # Parse Table 1 (Permeability Results)
        permeability_data = self._parse_permeability_table(tables[0], errors)

        # Parse Table 2 (Monolayer Integrity) if present
        integrity_data = {}
        if len(tables) >= 2:
            integrity_data = self._parse_integrity_table(tables[1], errors)

        # Build source metadata
        source = self.build_source_metadata(filepath)

        # Build results
        for compound_id, metrics in permeability_data.items():
            is_control = self.is_control_compound(compound_id)

            # Build results dict
            result_data = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "papp_unit": "1e-6 cm/s",
                **metrics,
                "source": source,
                "flags": [],
            }

            # Add monolayer integrity if available
            if compound_id in integrity_data:
                result_data["monolayer_integrity"] = integrity_data[compound_id]

            # Calculate derived fields
            self._calculate_derived_fields(result_data)

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

    def _parse_permeability_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[str, dict[str, Any]]:
        """
        Parse Table 1: Permeability Results.

        Returns:
            Dict mapping compound_id to metrics dict.
        """
        results = {}
        headers = table.get('headers', [])
        rows = table.get('data', [])

        # Find column indices
        # In NCU ADME files, compound ID is always the first column
        compound_idx = self.get_compound_column_index(headers)
        papp_ab_idx = self.get_column_index(headers, 'papp (a-b)', 'papp a-b', 'papp ab', 'ap-bl')
        papp_ba_idx = self.get_column_index(headers, 'papp (b-a)', 'papp b-a', 'papp ba', 'bl-ap')
        efflux_idx = self.get_column_index(headers, 'efflux ratio', 'efflux')
        recovery_ab_idx = self.get_column_index(headers, 'recovery (%) ap-bl', 'recovery ap-bl', 'recovery a-b', 'recovery ab')
        recovery_ba_idx = self.get_column_index(headers, 'recovery (%) bl-ap', 'recovery bl-ap', 'recovery b-a', 'recovery ba')

        for row_data in rows:
            values = row_data['values']

            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            metrics = {}

            if papp_ab_idx is not None and papp_ab_idx < len(values):
                val, flag = self.parse_value(values[papp_ab_idx])
                metrics['papp_ab'] = val

            if papp_ba_idx is not None and papp_ba_idx < len(values):
                val, flag = self.parse_value(values[papp_ba_idx])
                metrics['papp_ba'] = val

            if efflux_idx is not None and efflux_idx < len(values):
                val, flag = self.parse_value(values[efflux_idx])
                metrics['efflux_ratio'] = val

            if recovery_ab_idx is not None and recovery_ab_idx < len(values):
                val, flag = self.parse_value(values[recovery_ab_idx])
                metrics['recovery_ab_pct'] = val

            if recovery_ba_idx is not None and recovery_ba_idx < len(values):
                val, flag = self.parse_value(values[recovery_ba_idx])
                metrics['recovery_ba_pct'] = val

            results[compound_id] = metrics

        return results

    def _parse_integrity_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[str, dict[str, Any]]:
        """
        Parse Table 2: Monolayer Integrity.

        Returns:
            Dict mapping compound_id to integrity metrics dict.
        """
        results = {}
        headers = table.get('headers', [])
        rows = table.get('data', [])

        # Find column indices
        # In NCU ADME files, compound ID is always the first column
        compound_idx = self.get_compound_column_index(headers)
        teer_ab_idx = self.get_column_index(headers, 'teer a-b', 'teer ab', 'teer ap-bl')
        teer_ba_idx = self.get_column_index(headers, 'teer b-a', 'teer ba', 'teer bl-ap')
        ly_ab_idx = self.get_column_index(headers, 'ly leakage a-b', 'ly a-b', 'ly ap-bl')
        ly_ba_idx = self.get_column_index(headers, 'ly leakage b-a', 'ly b-a', 'ly bl-ap')

        for row_data in rows:
            values = row_data['values']

            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            integrity = {
                'teer_unit': 'ohm_cm2',
            }

            if teer_ab_idx is not None and teer_ab_idx < len(values):
                val, flag = self.parse_value(values[teer_ab_idx])
                integrity['teer_ab'] = val

            if teer_ba_idx is not None and teer_ba_idx < len(values):
                val, flag = self.parse_value(values[teer_ba_idx])
                integrity['teer_ba'] = val

            if ly_ab_idx is not None and ly_ab_idx < len(values):
                val, flag = self.parse_value(values[ly_ab_idx])
                integrity['ly_leakage_ab_pct'] = val

            if ly_ba_idx is not None and ly_ba_idx < len(values):
                val, flag = self.parse_value(values[ly_ba_idx])
                integrity['ly_leakage_ba_pct'] = val

            results[compound_id] = integrity

        return results

    def _calculate_derived_fields(self, result_data: dict) -> None:
        """Calculate derived classification fields."""
        papp_ab = result_data.get('papp_ab')
        efflux_ratio = result_data.get('efflux_ratio')

        # Permeability classification
        if papp_ab is not None:
            if papp_ab >= self.PAPP_HIGH_THRESHOLD:
                result_data['permeability_class'] = 'high'
            elif papp_ab >= self.PAPP_MEDIUM_THRESHOLD:
                result_data['permeability_class'] = 'medium'
            else:
                result_data['permeability_class'] = 'low'

        # Efflux substrate classification
        if efflux_ratio is not None:
            result_data['efflux_substrate'] = efflux_ratio > self.EFFLUX_THRESHOLD

    def _add_quality_flags(self, result_data: dict) -> None:
        """Add quality flags based on result values."""
        flags = result_data.get('flags', [])

        # Check permeability
        papp_ab = result_data.get('papp_ab')
        if papp_ab is not None and papp_ab < self.PAPP_MEDIUM_THRESHOLD:
            flags.append('poor_permeability')

        # Check efflux
        efflux_ratio = result_data.get('efflux_ratio')
        if efflux_ratio is not None:
            if efflux_ratio > self.HIGH_EFFLUX_THRESHOLD:
                flags.append('high_efflux')
            elif efflux_ratio > self.EFFLUX_THRESHOLD:
                flags.append('efflux_substrate')

        # Check recovery
        recovery_ab = result_data.get('recovery_ab_pct')
        recovery_ba = result_data.get('recovery_ba_pct')
        if recovery_ab is not None and recovery_ab < self.LOW_RECOVERY_THRESHOLD:
            flags.append('low_recovery')
        if recovery_ba is not None and recovery_ba < self.LOW_RECOVERY_THRESHOLD:
            flags.append('low_recovery')

        result_data['flags'] = list(set(flags))  # Deduplicate
