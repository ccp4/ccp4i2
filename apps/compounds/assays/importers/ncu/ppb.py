"""
NCU Plasma Protein Binding (PPB) Parser.

Parses ADME-NCU-PPB-YYYYMMDD.xlsx files. Reports the fraction of the
test compound bound to plasma proteins at equilibrium, along with
recovery and stability checks over the assay duration.

Summary layout:
    Table 1: Compound ID | Species | % Bound | % Unbound |
             % Recovery | % Remaining at 6 hr
"""

from pathlib import Path
from typing import Any

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class PlasmaProteinBindingParser(NCUBaseParser):
    assay_type = "plasma_protein_binding"
    assay_code = "PPB"
    protocol_slug = "ncu-ppb"
    kpi_field = "percent_bound"
    kpi_unit = "%"

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        return assay_type_from_file.lower() == "ppb"

    def parse(self, filepath: Path) -> ParseResult:
        errors: list[ValidationError] = []
        results: list[ParsedResult] = []

        try:
            sheets = self.read_excel(filepath)
        except Exception as e:
            return ParseResult(errors=[ValidationError(
                message=f"Failed to read Excel file: {e}",
                severity="error",
            )])

        summary, summary_name, all_names = self.find_summary_sheet_with_name(sheets)
        if summary is None:
            return ParseResult(
                errors=[ValidationError(
                    message=f"Summary sheet not found. Available sheets: {all_names}",
                    severity="error",
                )],
                metadata={"sheets_found": all_names},
            )

        tables = self.parse_summary_tables(summary, errors)
        if not tables:
            errors.append(ValidationError(
                message=f"No data tables found in sheet '{summary_name}'.",
                severity="error",
            ))
            return ParseResult(errors=errors, metadata={"sheets_found": all_names})

        headers = tables[0].get("headers", [])
        rows = tables[0].get("data", [])
        compound_idx = self.get_compound_column_index(headers)
        species_idx = self.get_column_index(headers, "species")
        bound_idx = self.get_column_index(headers, "% bound", "percent bound", "fraction bound")
        unbound_idx = self.get_column_index(headers, "% unbound", "percent unbound", "% free", "fraction unbound")
        recovery_idx = self.get_column_index(headers, "% recovery", "recovery")
        remaining_idx = self.get_column_index(headers, "% remaining", "remaining at")

        self.forward_fill_compound_ids(rows, compound_idx)

        source = self.build_source_metadata(filepath)
        for row_data in rows:
            values = row_data["values"]
            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue
            species = (values[species_idx] if species_idx is not None and species_idx < len(values) else "Unknown") or "Unknown"

            metrics: dict[str, Any] = {}

            def _take(idx, key):
                if idx is None or idx >= len(values):
                    return
                val, _flag = self.parse_value(values[idx])
                metrics[key] = val

            _take(bound_idx, "percent_bound")
            _take(unbound_idx, "percent_unbound")
            _take(recovery_idx, "percent_recovery")
            _take(remaining_idx, "percent_remaining_6h")

            result_data: dict[str, Any] = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "kpi_unit": self.kpi_unit,
                "species": species,
                **metrics,
                "source": source,
                "flags": [],
            }

            # Quality flag: very high binding (>99%) is often reported as "highly bound"
            if isinstance(metrics.get("percent_bound"), (int, float)) and metrics["percent_bound"] > 99:
                result_data["flags"].append("highly_bound")

            results.append(ParsedResult(
                compound_id=compound_id,
                species=species,
                results=result_data,
                is_control=self.is_control_compound(compound_id),
            ))

        errors.extend(self.validate(results))
        return ParseResult(
            results=results,
            errors=errors,
            metadata={
                "vendor": self.vendor,
                "assay_type": self.assay_type,
                "file": filepath.name,
                "sheets_found": all_names,
                "summary_sheet_used": summary_name,
                "tables_found": len(tables),
            },
        )
