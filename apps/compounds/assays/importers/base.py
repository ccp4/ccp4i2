"""
Base classes for ADME data importers.

This module defines the abstract base class and common data structures
used by all ADME parsers.

ADME Import vs Table-of-Values (TOV) Import
==========================================

The assay system supports two types of pre-analyzed data imports:

1. ADME Import (/assays/import-adme/)
   - For vendor-specific ADME data files (e.g., NCU Excel reports)
   - Parser auto-detects file type from filename pattern
   - KPI is IMPLICIT - determined by assay type (e.g., t1/2 for stability assays)
   - Compound column is ALWAYS the first column in data tables
   - No user configuration needed - just upload and confirm

2. Table-of-Values Import (/assays/import-tov/)
   - For generic spreadsheet data from any source
   - User must select a protocol and specify column mappings
   - KPI is EXPLICIT - user selects a column containing the KPI field name
     (e.g., a "KPI" column where every row says "EC50")
   - User explicitly selects which column contains compound identifiers
   - More flexible but requires more user configuration

Key Differences:
- ADME: kpi_field is a parser attribute, used directly
- TOV: kpi_column contains a field name that must match another column header
- ADME: compound_id always in first column (standardized vendor format)
- TOV: compound_column is user-selected (variable spreadsheet formats)
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Optional


@dataclass
class ParsedResult:
    """
    Represents a single parsed result from an ADME data file.

    Each ParsedResult corresponds to one DataSeries + AnalysisResult pair
    that will be created in the database.
    """

    compound_id: str
    """Original compound identifier from the data file."""

    results: dict[str, Any]
    """JSON-serializable dict to store in AnalysisResult.results."""

    species: Optional[str] = None
    """Species for multi-species assays (e.g., 'Human', 'Mouse')."""

    assay_format: Optional[str] = None
    """Assay format/condition (e.g., '+Cofactors', '+GSH')."""

    source_row: Optional[int] = None
    """Row number in source file for traceability."""

    is_control: bool = False
    """True if this is a control compound (should be skipped for import)."""

    def __post_init__(self):
        # Ensure results is a dict
        if self.results is None:
            self.results = {}


@dataclass
class ValidationError:
    """
    Represents a validation error or warning during parsing.
    """

    message: str
    """Human-readable error message."""

    row: Optional[int] = None
    """Row number in source file (if applicable)."""

    column: Optional[str] = None
    """Column name (if applicable)."""

    severity: str = "error"
    """Severity level: 'error', 'warning', or 'info'."""

    def __str__(self) -> str:
        location = ""
        if self.row is not None:
            location = f"Row {self.row}"
        if self.column:
            location = f"{location}, Column '{self.column}'" if location else f"Column '{self.column}'"
        if location:
            return f"[{self.severity.upper()}] {location}: {self.message}"
        return f"[{self.severity.upper()}] {self.message}"


@dataclass
class ParseResult:
    """
    Container for the complete result of parsing a file.
    """

    results: list[ParsedResult] = field(default_factory=list)
    """Successfully parsed results."""

    errors: list[ValidationError] = field(default_factory=list)
    """Validation errors and warnings."""

    metadata: dict[str, Any] = field(default_factory=dict)
    """File-level metadata (vendor, date, assay type, etc.)."""

    @property
    def success(self) -> bool:
        """Returns True if there were no errors (warnings allowed)."""
        return not any(e.severity == "error" for e in self.errors)

    @property
    def compound_count(self) -> int:
        """Number of unique compounds parsed."""
        return len(set(r.compound_id for r in self.results if not r.is_control))


class ADMEParser(ABC):
    """
    Abstract base class for ADME data parsers.

    Each parser handles a specific vendor/assay type combination.
    Parsers are responsible for:
    1. Detecting if they can handle a given file
    2. Parsing the file into structured results
    3. Validating the parsed data

    The kpi_field attribute defines the key performance indicator for this
    assay type. Unlike Table-of-Values imports where the user specifies
    which column contains the KPI, ADME parsers have an implicit KPI
    determined by the assay type (e.g., t1/2 for stability assays,
    efflux ratio for Caco-2 permeability).
    """

    # Subclasses must define these
    vendor: str = ""
    """Vendor identifier (e.g., 'NCU')."""

    assay_type: str = ""
    """Assay type identifier (e.g., 'liver_microsome_stability')."""

    assay_code: str = ""
    """Short code for assay type (e.g., 'LM')."""

    protocol_slug: str = ""
    """Protocol slug for database matching (e.g., 'ncu-lm')."""

    kpi_field: str = ""
    """
    Primary KPI field name (e.g., 't1_2_min').

    This is the key quantity for this assay type, used as the primary
    metric for compound comparison. Unlike Table-of-Values imports
    where the KPI is user-specified per row, ADME imports have a
    fixed KPI determined by the assay type.
    """

    # Known control compounds to skip
    CONTROL_COMPOUNDS = {
        'verapamil', 'metoprolol', 'digoxin', 'propantheline',
        'afatinib', 'warfarin', 'testosterone', 'midazolam',
        'propranolol', 'atenolol', 'caffeine', 'lucifer yellow',
    }

    @abstractmethod
    def detect(self, filepath: Path) -> bool:
        """
        Determine if this parser can handle the given file.

        Args:
            filepath: Path to the data file.

        Returns:
            True if this parser can handle the file, False otherwise.
        """
        pass

    @abstractmethod
    def parse(self, filepath: Path) -> ParseResult:
        """
        Parse the data file and return structured results.

        Args:
            filepath: Path to the data file.

        Returns:
            ParseResult containing parsed data and any errors.
        """
        pass

    def validate(self, results: list[ParsedResult]) -> list[ValidationError]:
        """
        Validate parsed results.

        Override in subclasses to add assay-specific validation.

        Args:
            results: List of parsed results to validate.

        Returns:
            List of validation errors/warnings.
        """
        errors = []

        for result in results:
            # Check for required fields
            if not result.compound_id:
                errors.append(ValidationError(
                    message="Missing compound ID",
                    row=result.source_row,
                    severity="error"
                ))

            # Check results dict
            if not result.results:
                errors.append(ValidationError(
                    message=f"No results data for {result.compound_id}",
                    row=result.source_row,
                    severity="warning"
                ))

        return errors

    def is_control_compound(self, compound_id: str) -> bool:
        """
        Check if a compound is a known control.

        Args:
            compound_id: Compound identifier to check.

        Returns:
            True if this is a control compound.
        """
        if not compound_id:
            return False
        return compound_id.lower().strip() in self.CONTROL_COMPOUNDS

    def extract_date_from_filename(self, filepath: Path) -> Optional[date]:
        """
        Extract date from filename in YYYYMMDD format.

        Args:
            filepath: Path to the file.

        Returns:
            Extracted date or None if not found.
        """
        import re

        # Look for YYYYMMDD pattern
        match = re.search(r'(\d{4})(\d{2})(\d{2})', filepath.stem)
        if match:
            try:
                year = int(match.group(1))
                month = int(match.group(2))
                day = int(match.group(3))
                return date(year, month, day)
            except ValueError:
                pass
        return None

    def build_source_metadata(self, filepath: Path) -> dict[str, Any]:
        """
        Build source metadata dict for inclusion in results.

        Args:
            filepath: Path to the source file.

        Returns:
            Dict with vendor, file, and date information.
        """
        file_date = self.extract_date_from_filename(filepath)
        return {
            "vendor": self.vendor,
            "file": filepath.name,
            "date": file_date.isoformat() if file_date else None,
        }
