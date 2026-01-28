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


def clear_pattern_cache():
    """
    Clear the cached pattern string.

    Call this if settings change at runtime (mainly for testing).
    """
    _get_compound_pattern_str.cache_clear()
