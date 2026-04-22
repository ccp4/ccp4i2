"""
NCU LogD Parser.

Parses ADME-NCU-LogD-YYYYMMDD.xlsx files. LogD is the measured
distribution coefficient between 1-octanol and phosphate buffer at
pH 7.4 — a dimensionless log value.

Summary layout:
    Table 1: Compound ID | LogD
"""

from pathlib import Path
from typing import Any

from compounds.assays.kpi_utils import UNITLESS

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class LogDParser(NCUBaseParser):
    assay_type = "logd_ph74"
    assay_code = "LogD"
    protocol_slug = "ncu-logd"
    kpi_field = "log_d"
    kpi_unit = UNITLESS

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        return "logd" in assay_type_from_file.lower().replace("-", "").replace("_", "")

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
        logd_idx = self.get_column_index(headers, "logd", "log d")

        self.forward_fill_compound_ids(rows, compound_idx)

        source = self.build_source_metadata(filepath)
        for row_data in rows:
            values = row_data["values"]
            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue

            logd_val, logd_flag = (None, None)
            if logd_idx is not None and logd_idx < len(values):
                logd_val, logd_flag = self.parse_value(values[logd_idx])

            result_data: dict[str, Any] = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "kpi_unit": self.kpi_unit,
                "log_d": logd_val,
                "source": source,
                "flags": [logd_flag] if logd_flag else [],
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
