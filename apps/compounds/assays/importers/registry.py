"""
Parser Registry for ADME Data Importers.

This module provides functions for registering, discovering, and
invoking ADME parsers based on file characteristics.
"""

from pathlib import Path
from typing import Optional, Type

from .base import ADMEParser, ParseResult

# Registry of all available parsers
_PARSERS: list[Type[ADMEParser]] = []


def register_parser(parser_class: Type[ADMEParser]) -> Type[ADMEParser]:
    """
    Decorator to register a parser class.

    Usage:
        @register_parser
        class MyParser(ADMEParser):
            ...
    """
    if parser_class not in _PARSERS:
        _PARSERS.append(parser_class)
    return parser_class


def list_parsers() -> list[Type[ADMEParser]]:
    """
    List all registered parser classes.

    Returns:
        List of parser classes.
    """
    # Ensure NCU parsers are loaded
    _ensure_parsers_loaded()
    return list(_PARSERS)


def detect_parser(filepath: Path | str) -> Optional[ADMEParser]:
    """
    Auto-detect the appropriate parser for a file.

    Args:
        filepath: Path to the data file.

    Returns:
        Parser instance if one can handle the file, None otherwise.
    """
    filepath = Path(filepath)

    # Ensure parsers are loaded
    _ensure_parsers_loaded()

    for parser_class in _PARSERS:
        parser = parser_class()
        if parser.detect(filepath):
            return parser

    return None


def get_parser(
    vendor: str,
    assay_type: Optional[str] = None,
    assay_code: Optional[str] = None,
) -> Optional[ADMEParser]:
    """
    Get a specific parser by vendor and assay type.

    Args:
        vendor: Vendor identifier (e.g., 'NCU').
        assay_type: Assay type identifier (e.g., 'liver_microsome_stability').
        assay_code: Short assay code (e.g., 'LM').

    Returns:
        Parser instance if found, None otherwise.
    """
    _ensure_parsers_loaded()

    vendor_lower = vendor.lower()

    for parser_class in _PARSERS:
        if parser_class.vendor.lower() != vendor_lower:
            continue

        if assay_type and parser_class.assay_type == assay_type:
            return parser_class()

        if assay_code and parser_class.assay_code.upper() == assay_code.upper():
            return parser_class()

    return None


def get_parser_by_slug(protocol_slug: str) -> Optional[ADMEParser]:
    """
    Get a specific parser by its protocol_slug.

    Args:
        protocol_slug: Protocol slug (e.g., 'ncu-bs', 'ncu-lm').

    Returns:
        Parser instance if found, None otherwise.
    """
    _ensure_parsers_loaded()

    slug_lower = protocol_slug.lower()

    for parser_class in _PARSERS:
        if parser_class.protocol_slug.lower() == slug_lower:
            return parser_class()

    return None


def parse_adme_file(filepath: Path | str) -> ParseResult:
    """
    Parse an ADME data file using auto-detected parser.

    Args:
        filepath: Path to the data file.

    Returns:
        ParseResult containing parsed data and any errors.

    Raises:
        ValueError: If no parser can handle the file.
    """
    filepath = Path(filepath)

    parser = detect_parser(filepath)
    if parser is None:
        from .base import ParseResult, ValidationError
        return ParseResult(
            errors=[ValidationError(
                message=f"No parser found for file: {filepath.name}",
                severity="error"
            )]
        )

    return parser.parse(filepath)


def _ensure_parsers_loaded():
    """
    Ensure all parser modules are imported.

    This triggers the @register_parser decorators to run.
    """
    # Import NCU parsers to register them
    try:
        from . import ncu  # noqa: F401
    except ImportError:
        pass
