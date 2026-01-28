"""
Shared utilities for the compounds app.
"""

import logging
import os
import re
from typing import Optional, TYPE_CHECKING

from compounds.formatting import (
    get_compound_pattern,
    get_malformed_pattern,
)

if TYPE_CHECKING:
    from compounds.registry.models import Compound

logger = logging.getLogger(__name__)


def resolve_compound(identifier: str) -> Optional['Compound']:
    """
    Resolve a compound identifier to a registered Compound instance.

    Attempts matching using multiple strategies, from most specific to least:

    1. Full NCL format: NCL-00026042, ncl-00026042, NCL-26042
    2. Plain numeric: 26042, 00026042
    3. Partial NCL patterns: NCL26042, ncl26042, NCL 26042
    4. Supplier reference lookup
    5. Barcode format: {initials}-{labbook}-{page}
    6. Alias lookup (case-insensitive)

    Args:
        identifier: The compound identifier string from an import file.
                   Can be NCL format, numeric, supplier ref, or alias.

    Returns:
        Compound instance if matched, None otherwise.

    Examples:
        >>> resolve_compound("NCL-00026042")  # Full NCL format
        <Compound: NCL-00026042>
        >>> resolve_compound("ncl-26042")     # Lowercase, no leading zeros
        <Compound: NCL-00026042>
        >>> resolve_compound("26042")         # Plain numeric
        <Compound: NCL-00026042>
        >>> resolve_compound("AST-123-45")    # Barcode format
        <Compound: NCL-00026042>
        >>> resolve_compound("compound-a")    # Alias match
        <Compound: NCL-00026042>
    """
    # Import here to avoid circular imports
    from compounds.registry.models import Compound, Supplier

    if not identifier:
        return None

    identifier = str(identifier).strip()
    if not identifier:
        return None

    # Strategy 1: Handle malformed patterns FIRST like "PREFIX000-27421" (dash after zeros)
    # Must be checked before standard pattern to avoid incorrect matches
    malformed_pattern = get_malformed_pattern(capturing=True)
    ncl_malformed = malformed_pattern.match(identifier)
    if ncl_malformed:
        try:
            reg_number = int(ncl_malformed.group(1))
            compound = Compound.objects.filter(reg_number=reg_number).first()
            if compound:
                return compound
        except ValueError:
            pass

    # Strategy 2: Standard format patterns (handles PREFIX-XXXXXXXX, PREFIX-XX, etc.)
    # Pattern matches: PREFIX-00026042, PREFIX-26042, PREFIX26042, PREFIX 26042
    standard_pattern = get_compound_pattern(capturing=True)
    ncl_match = standard_pattern.match(identifier)
    if ncl_match:
        try:
            reg_number = int(ncl_match.group(1))
            compound = Compound.objects.filter(reg_number=reg_number).first()
            if compound:
                return compound
        except ValueError:
            pass

    # Strategy 3: Plain numeric (could be reg_number without NCL prefix)
    if identifier.isdigit():
        try:
            reg_number = int(identifier)
            compound = Compound.objects.filter(reg_number=reg_number).first()
            if compound:
                return compound
        except ValueError:
            pass

    # Strategy 4: More permissive extraction (handles embedded compound IDs)
    # E.g., "sample PREFIX-26042 test", "PREFIX_26042"
    # First try malformed pattern for embedded malformed IDs
    ncl_malformed_search = malformed_pattern.search(identifier)
    if ncl_malformed_search:
        try:
            reg_number = int(ncl_malformed_search.group(1))
            compound = Compound.objects.filter(reg_number=reg_number).first()
            if compound:
                return compound
        except ValueError:
            pass

    # Then try standard embedded pattern
    ncl_extract = standard_pattern.search(identifier)
    if ncl_extract:
        try:
            reg_number = int(ncl_extract.group(1))
            compound = Compound.objects.filter(reg_number=reg_number).first()
            if compound:
                return compound
        except ValueError:
            pass

    # Strategy 4: Supplier reference lookup (case-insensitive)
    compound = Compound.objects.filter(supplier_ref__iexact=identifier).first()
    if compound:
        return compound

    # Strategy 5: Barcode format - {supplier_initials}-{labbook}-{page}
    # E.g., "AST-123-45" or "NCL-100-5"
    barcode_match = re.match(r'^([A-Za-z]+)-(\d+)-(\d+)$', identifier)
    if barcode_match:
        initials, labbook, page = barcode_match.groups()
        # Find supplier by initials
        supplier = Supplier.objects.filter(initials__iexact=initials).first()
        if supplier:
            compound = Compound.objects.filter(
                supplier=supplier,
                labbook_number=int(labbook),
                page_number=int(page)
            ).first()
            if compound:
                return compound

    # Strategy 6: Alias lookup (case-insensitive)
    # Search through all compounds' aliases JSONField
    identifier_lower = identifier.lower()
    # Use a database query that searches the JSON array
    # PostgreSQL: aliases @> '["value"]' or use __contains
    # For SQLite compatibility, we'll do a case-insensitive search
    for compound in Compound.objects.exclude(aliases=[]).exclude(aliases__isnull=True):
        if compound.aliases:
            for alias in compound.aliases:
                if str(alias).lower() == identifier_lower:
                    return compound

    return None


def resolve_compound_batch(identifiers: list[str]) -> dict[str, Optional['Compound']]:
    """
    Resolve multiple compound identifiers efficiently.

    For large batches, this is more efficient than calling resolve_compound()
    repeatedly as it pre-fetches data and caches lookups.

    Args:
        identifiers: List of compound identifier strings.

    Returns:
        Dictionary mapping each identifier to its resolved Compound (or None).
    """
    from compounds.registry.models import Compound, Supplier

    results = {}
    remaining = []

    # First pass: Extract compound IDs and plain numerics
    # Note: malformed pattern must be tried FIRST to avoid incorrect matches
    malformed_pattern = get_malformed_pattern(capturing=True)
    standard_pattern = get_compound_pattern(capturing=True)
    reg_numbers_needed = set()
    identifier_to_reg = {}

    for identifier in identifiers:
        if not identifier:
            results[identifier] = None
            continue

        identifier = str(identifier).strip()
        if not identifier:
            results[identifier] = None
            continue

        # Try malformed pattern FIRST (PREFIX000-27421)
        # Must check before standard pattern to avoid 0* consuming zeros incorrectly
        match = malformed_pattern.search(identifier)
        if match:
            reg_num = int(match.group(1))
            reg_numbers_needed.add(reg_num)
            identifier_to_reg[identifier] = reg_num
            continue

        # Try standard pattern (PREFIX-XXXXXXXX, etc.)
        match = standard_pattern.search(identifier)
        if match:
            reg_num = int(match.group(1))
            reg_numbers_needed.add(reg_num)
            identifier_to_reg[identifier] = reg_num
            continue

        # Try plain numeric
        if identifier.isdigit():
            reg_num = int(identifier)
            reg_numbers_needed.add(reg_num)
            identifier_to_reg[identifier] = reg_num
            continue

        # No compound ID pattern found - add to remaining for fallback strategies
        remaining.append(identifier)

    # Batch fetch compounds by reg_number
    if reg_numbers_needed:
        compounds_by_reg = {
            c.reg_number: c
            for c in Compound.objects.filter(reg_number__in=reg_numbers_needed)
        }
        for identifier, reg_num in identifier_to_reg.items():
            results[identifier] = compounds_by_reg.get(reg_num)

    # Handle remaining identifiers with fallback strategies
    if remaining:
        # Batch fetch by supplier_ref
        supplier_ref_map = {
            c.supplier_ref.lower(): c
            for c in Compound.objects.filter(supplier_ref__in=remaining)
            if c.supplier_ref
        }

        # Pre-fetch compounds with aliases for alias matching
        compounds_with_aliases = list(
            Compound.objects.exclude(aliases=[]).exclude(aliases__isnull=True)
        )

        # Pre-fetch suppliers for barcode matching
        suppliers_by_initials = {
            s.initials.upper(): s
            for s in Supplier.objects.exclude(initials__isnull=True).exclude(initials='')
        }

        barcode_pattern = re.compile(r'^([A-Za-z]+)-(\d+)-(\d+)$')

        for identifier in remaining:
            # Try supplier_ref
            if identifier.lower() in supplier_ref_map:
                results[identifier] = supplier_ref_map[identifier.lower()]
                continue

            # Try barcode
            barcode_match = barcode_pattern.match(identifier)
            if barcode_match:
                initials, labbook, page = barcode_match.groups()
                supplier = suppliers_by_initials.get(initials.upper())
                if supplier:
                    compound = Compound.objects.filter(
                        supplier=supplier,
                        labbook_number=int(labbook),
                        page_number=int(page)
                    ).first()
                    if compound:
                        results[identifier] = compound
                        continue

            # Try alias
            identifier_lower = identifier.lower()
            found = False
            for compound in compounds_with_aliases:
                if compound.aliases:
                    for alias in compound.aliases:
                        if str(alias).lower() == identifier_lower:
                            results[identifier] = compound
                            found = True
                            break
                if found:
                    break

            if identifier not in results:
                results[identifier] = None

    return results


def delete_file_field(file_field, save=False):
    """
    Delete a file from a FileField, handling both local and cloud storage.

    This function safely deletes files regardless of the storage backend:
    - Local filesystem: Uses os.remove() after checking file exists
    - Cloud storage (Azure Blob, S3, etc.): Uses storage backend's delete()

    Args:
        file_field: A Django FileField instance (e.g., instance.genbank_file)
        save: Whether to save the model after clearing the field (default False)

    Returns:
        True if file was deleted, False if no file to delete
    """
    if not file_field:
        return False

    try:
        # Try local filesystem deletion first
        file_path = file_field.path
        if os.path.isfile(file_path):
            os.remove(file_path)
            logger.debug(f"Deleted local file: {file_path}")
            return True
    except NotImplementedError:
        # Cloud storage backends don't support .path
        # Use the storage backend's delete method instead
        try:
            file_field.delete(save=save)
            logger.debug(f"Deleted cloud file: {file_field.name}")
            return True
        except Exception as e:
            logger.warning(f"Failed to delete cloud file {file_field.name}: {e}")
            return False
    except Exception as e:
        logger.warning(f"Failed to delete file: {e}")
        return False

    return False
