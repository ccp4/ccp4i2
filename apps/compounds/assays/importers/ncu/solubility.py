"""
NCU Solubility Parser.

Parses ADME-NCU-Solubility-YYYYMMDD.xlsx files. Kinetic solubility in
phosphate buffer at pH 7.4, reported in μM. The assay has an upper
cut-off (typically 300 μM) above which values are reported as e.g.
">300" — handled via the base class's parse_value() special-value
logic which flags the row with `exceeds_window`.

Summary layout:
    Table 1: Compound ID | Solubility in phosphate buffer pH 7.4 (μM)
"""

from pathlib import Path
from typing import Any

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class SolubilityParser(NCUBaseParser):
    assay_type = "kinetic_solubility_ph74"
    assay_code = "Solubility"
    protocol_slug = "ncu-solubility"
    kpi_field = "solubility_um"
    kpi_unit = "uM"

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        return "solubility" in assay_type_from_file.lower()

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
        sol_idx = self.get_column_index(headers, "solubility")

        self.forward_fill_compound_ids(rows, compound_idx)

        source = self.build_source_metadata(filepath)
        for row_data in rows:
            values = row_data["values"]
            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            sol_val, sol_flag = (None, None)
            if sol_idx is not None and sol_idx < len(values):
                sol_val, sol_flag = self.parse_value(values[sol_idx])

            flags: list[str] = []
            if sol_flag:
                flags.append(sol_flag)
            if isinstance(sol_val, (int, float)):
                if sol_val < 10:
                    flags.append("low_solubility")
                elif sol_val >= 300:
                    flags.append("highly_soluble")

            result_data: dict[str, Any] = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "kpi_unit": self.kpi_unit,
                "solubility_um": sol_val,
                "source": source,
                "flags": flags,
            }

            results.append(ParsedResult(
                compound_id=compound_id,
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
