"""
Compound ID formatting utilities.

Provides centralized functions for formatting and parsing compound identifiers
based on deployment configuration. This ensures consistent handling across
the entire codebase.

Configuration is via environment variables:
- COMPOUND_ID_PREFIX: The prefix (default: "NCL")
- COMPOUND_ID_DIGITS: Number of digits to pad to (default: 8)

Example:
    COMPOUND_ID_PREFIX=NCL, COMPOUND_ID_DIGITS=8
    reg_number=26042 -> "NCL-00026042"

    COMPOUND_ID_PREFIX=NCLAK, COMPOUND_ID_DIGITS=6
    reg_number=123 -> "NCLAK-000123"
"""

import re
from functools import lru_cache
from typing import Optional

from django.conf import settings


def get_compound_prefix() -> str:
    """Get the configured compound ID prefix."""
    return getattr(settings, 'COMPOUND_ID_PREFIX', 'NCL')


def get_compound_digits() -> int:
    """Get the configured number of digits for compound IDs."""
    return getattr(settings, 'COMPOUND_ID_DIGITS', 8)


def format_compound_id(reg_number: int) -> str:
    """
    Format a registration number as a compound ID.

    Args:
        reg_number: The integer registration number.

    Returns:
        Formatted compound ID (e.g., "NCL-00026042").
    """
    prefix = get_compound_prefix()
    digits = get_compound_digits()
    return f'{prefix}-{reg_number:0{digits}d}'


@lru_cache(maxsize=1)
def _get_compound_pattern_str() -> str:
    """Get the regex pattern string for compound IDs (cached)."""
    prefix = get_compound_prefix()
    # Escape any special regex characters in prefix
    escaped_prefix = re.escape(prefix)
    return escaped_prefix


def get_compound_pattern(capturing: bool = True) -> re.Pattern:
    """
    Get a regex pattern that matches compound IDs with the configured prefix.

    Handles standard formats:
    - PREFIX-00026042 (standard with zeros)
    - PREFIX-26042 (without leading zeros)
    - PREFIX26042 (no dash)
    - PREFIX 26042 (space instead of dash)
    - PREFIX_26042 (underscore instead of dash)

    Args:
        capturing: If True, captures the numeric portion in group(1).

    Returns:
        Compiled regex pattern (case-insensitive).
    """
    prefix = _get_compound_pattern_str()
    if capturing:
        return re.compile(rf'{prefix}[-_\s]?0*(\d+)', re.IGNORECASE)
    return re.compile(rf'{prefix}[-_\s]?\d+', re.IGNORECASE)


def get_malformed_pattern(capturing: bool = True) -> re.Pattern:
    """
    Get a regex pattern for malformed compound IDs.

    Handles user typos where the dash comes AFTER leading zeros:
    - PREFIX000-27421 (should be PREFIX-00027421)
    - PREFIX00-187 (should be PREFIX-00000187)

    Args:
        capturing: If True, captures the numeric portion in group(1).

    Returns:
        Compiled regex pattern (case-insensitive).
    """
    prefix = _get_compound_pattern_str()
    if capturing:
        return re.compile(rf'{prefix}0+-(\d+)', re.IGNORECASE)
    return re.compile(rf'{prefix}0+-\d+', re.IGNORECASE)


def get_exact_pattern() -> re.Pattern:
    """
    Get a regex pattern for exact compound ID matches (anchored).

    For validating complete strings, not searching within text.
    Handles both standard and malformed formats.

    Returns:
        Compiled regex pattern (case-insensitive).
    """
    prefix = _get_compound_pattern_str()
    # Match either malformed (PREFIX000-123) or standard (PREFIX-00123, PREFIX123, etc.)
    return re.compile(
        rf'^(?:{prefix}0+-(\d+)|{prefix}[-_\s]?0*(\d+))$',
        re.IGNORECASE
    )


def extract_reg_number(identifier: str) -> Optional[int]:
    """
    Extract registration number from a compound identifier string.

    Handles various formats users might input:
    - Full format: PREFIX-00026042
    - Compact: PREFIX26042
    - With spaces: PREFIX 26042
    - Malformed: PREFIX000-27421 (dash after zeros)
    - Plain numeric: 26042

    Args:
        identifier: The compound identifier string.

    Returns:
        Integer registration number, or None if not parseable.
    """
    if not identifier:
        return None

    identifier = str(identifier).strip()
    if not identifier:
        return None

    # Try malformed pattern first (PREFIX000-27421)
    # Must check before standard pattern to avoid incorrect matches
    malformed = get_malformed_pattern(capturing=True)
    match = malformed.match(identifier)
    if match:
        return int(match.group(1))

    # Try standard pattern (PREFIX-00026042, PREFIX26042, etc.)
    standard = get_compound_pattern(capturing=True)
    match = standard.match(identifier)
    if match:
        return int(match.group(1))

    # Fall back to plain numeric
    if identifier.isdigit():
        return int(identifier)

    return None


def search_for_compound_id(text: str) -> Optional[int]:
    """
    Search for a compound ID anywhere in text.

    Unlike extract_reg_number(), this searches within longer strings:
    - "sample NCL-26042 test" -> 26042
    - "NCL_26042" -> 26042
    - "text NCL000-27421 more" -> 27421

    Args:
        text: Text that may contain a compound ID.

    Returns:
        Integer registration number if found, None otherwise.
    """
    if not text:
        return None

    text = str(text).strip()
    if not text:
        return None

    # Try malformed pattern first
    malformed = get_malformed_pattern(capturing=True)
    match = malformed.search(text)
    if match:
        return int(match.group(1))

    # Try standard pattern
    standard = get_compound_pattern(capturing=True)
    match = standard.search(text)
    if match:
        return int(match.group(1))

    return None


def parse_compound_list(text: str) -> list[int]:
    """
    Parse a flexible list of compound identifiers from user input.

    Handles various input formats:
    - Whitespace-separated: "NCL-00035625 NCL-30282 56785"
    - Comma-separated: "NCL-00035625, NCL-30282, 56785"
    - Mixed: "NCL-00035625 ncl-30282,56785"
    - Case-insensitive prefixes: "NCL-123", "ncl-123", "Ncl-123"
    - Various prefix formats: "NCL-00123", "NCL123", "NCL 123"
    - Bare registration numbers: "56785"
    - Malformed: "NCL000-27421" (dash after zeros)

    Args:
        text: User input string containing compound identifiers.

    Returns:
        List of unique registration numbers (integers), preserving input order.
    """
    if not text:
        return []

    text = str(text).strip()
    if not text:
        return []

    # Collect all matches with their positions to preserve input order
    matches: list[tuple[int, int]] = []  # (position, reg_number)
    matched_spans: list[tuple[int, int]] = []  # Track matched regions

    # Find malformed pattern matches (PREFIX000-27421)
    malformed = get_malformed_pattern(capturing=True)
    for match in malformed.finditer(text):
        reg_number = int(match.group(1))
        matches.append((match.start(), reg_number))
        matched_spans.append((match.start(), match.end()))

    # Find standard pattern matches (PREFIX-00123, PREFIX123, etc.)
    standard = get_compound_pattern(capturing=True)
    for match in standard.finditer(text):
        # Skip if this region was already matched by malformed pattern
        if any(start <= match.start() < end for start, end in matched_spans):
            continue
        reg_number = int(match.group(1))
        matches.append((match.start(), reg_number))
        matched_spans.append((match.start(), match.end()))

    # Find bare numbers in unmatched regions
    # Split by whitespace and commas, track position
    bare_number = re.compile(r'\b(\d+)\b')
    for match in bare_number.finditer(text):
        # Skip if this region was already matched by prefixed patterns
        if any(start <= match.start() < end for start, end in matched_spans):
            continue
        reg_number = int(match.group(1))
        matches.append((match.start(), reg_number))

    # Sort by position and deduplicate while preserving order
    matches.sort(key=lambda x: x[0])
    result: list[int] = []
    seen: set[int] = set()
    for _, reg_number in matches:
        if reg_number not in seen:
            result.append(reg_number)
            seen.add(reg_number)

    return result


def clear_pattern_cache():
    """
    Clear the cached pattern string.

    Call this if settings change at runtime (mainly for testing).
    """
    _get_compound_pattern_str.cache_clear()
