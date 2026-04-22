"""
NCU Hepatocyte Stability Parser.

Parses ADME-NCU-Hepatocytes-YYYYMMDD.xlsx files containing metabolic
stability data from human and mouse hepatocytes.

Summary layout (same NCU pattern as LM):
    Table 1: Compound ID | Species | in vitro t1/2 (min) |
             in vitro CLint (μL/min/10⁶ cells) | Scale-up CLint (mL/min/kg) |
             Predicted Hepatic CLH (mL/min/kg) | Hepatic Extraction Ratio (ER)
    Table 2: Remaining percentage time course (optional — not modelled here)

KPI choice: ``t1_2_min`` with unit ``min``. CLint is reported in
μL/min/10⁶ cells — that unit is not in VALID_UNITS, so we use t1/2 as
the KPI and keep CLint (both raw and scaled) in the results for
downstream reporting.
"""

from pathlib import Path
from typing import Any

from ..base import ParseResult, ParsedResult, ValidationError
from ..registry import register_parser
from .base import NCUBaseParser


@register_parser
class HepatocyteStabilityParser(NCUBaseParser):
    assay_type = "hepatocyte_stability"
    assay_code = "Hepatocytes"
    protocol_slug = "ncu-hepatocytes"
    kpi_field = "t1_2_min"
    kpi_unit = "min"

    def _matches_assay_type(self, assay_type_from_file: str) -> bool:
        return "hepatocyte" in assay_type_from_file.lower()

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

        summary_data = self._parse_summary_table(tables[0], errors)

        source = self.build_source_metadata(filepath)
        for (compound_id, species), metrics in summary_data.items():
            result_data = {
                "assay_type": self.assay_type,
                "KPI": self.kpi_field,
                "kpi_unit": self.kpi_unit,
                "species": species,
                **metrics,
                "source": source,
                "flags": [],
            }
            self._add_quality_flags(result_data)
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

    def _parse_summary_table(
        self,
        table: dict,
        errors: list[ValidationError],
    ) -> dict[tuple[str, str], dict[str, Any]]:
        results: dict[tuple[str, str], dict[str, Any]] = {}
        headers = table.get("headers", [])
        rows = table.get("data", [])

        compound_idx = self.get_compound_column_index(headers)
        species_idx = self.get_column_index(headers, "species")
        t12_idx = self.get_column_index(headers, "t1/2", "t½", "half-life")
        clint_idx = self.get_column_index(headers, "in vitro clint", "clint (μl", "clint (ul")
        scaled_clint_idx = self.get_column_index(headers, "scale-up clint", "scaled clint")
        clh_idx = self.get_column_index(headers, "clh", "hepatic cl")
        er_idx = self.get_column_index(headers, "extraction ratio", "er")

        self.forward_fill_compound_ids(rows, compound_idx)

        for row_data in rows:
            values = row_data["values"]
            compound_id = values[compound_idx] if compound_idx < len(values) else None
            if not compound_id:
                continue
            species = (values[species_idx] if species_idx is not None and species_idx < len(values) else "Unknown") or "Unknown"

            metrics: dict[str, Any] = {}

            def _take(idx, key, flag_prefix=None):
                if idx is None or idx >= len(values):
                    return
                val, flag = self.parse_value(values[idx])
                metrics[key] = val
                if flag and flag_prefix:
                    metrics.setdefault("_flags", []).append(f"{flag_prefix}:{flag}")

            _take(t12_idx, "t1_2_min", flag_prefix="t1_2")
            _take(clint_idx, "clint_ul_min_Mcells")
            _take(scaled_clint_idx, "clint_scaled_ml_min_kg")
            _take(clh_idx, "clh_predicted_ml_min_kg")
            _take(er_idx, "extraction_ratio")

            results[(compound_id, species)] = metrics

        return results

    def _add_quality_flags(self, result_data: dict) -> None:
        flags = result_data.get("flags", [])
        t12 = result_data.get("t1_2_min")
        if isinstance(t12, (int, float)):
            if t12 < 10:
                flags.append("unstable")
            elif t12 > 120:
                flags.append("stable")
        clint = result_data.get("clint_ul_min_Mcells")
        if isinstance(clint, (int, float)) and clint > 50:
            flags.append("high_clearance")
        if "_flags" in result_data:
            flags.extend(result_data.pop("_flags"))
        result_data["flags"] = list(set(flags))
