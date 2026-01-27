"""
Custom filters for compound registry search.

Provides NCL-aware search that properly handles compound identifier formats
like NCL-00018710, NCL-000187, 18710, etc.
"""

import re
from django.db.models import Q
from rest_framework import filters


class CompoundSearchFilter(filters.SearchFilter):
    """
    Custom search filter that properly handles NCL compound identifier formats.

    The default SearchFilter uses icontains which doesn't work well with
    IntegerField reg_numbers. This filter:

    1. Detects NCL-format queries (NCL-XXXXX, NCL-000XXXXX, etc.)
    2. Extracts the numeric reg_number and searches by exact or prefix match
    3. ALSO searches text fields (supplier_ref, comments, aliases) for ALL queries
    4. Combines reg_number and text matches with OR

    Examples:
        - "NCL-000187" → reg_numbers starting with 187 OR text fields containing "NCL-000187"
        - "NCL-00018710" → reg_number 18710 OR text fields containing "NCL-00018710"
        - "18710" → reg_number 18710 OR text fields containing "18710"
        - "aspirin" → text fields containing "aspirin"
    """

    def filter_queryset(self, request, queryset, view):
        search_terms = self.get_search_terms(request)

        if not search_terms:
            return queryset

        # For compound searches, handle NCL format specially
        search_query = ' '.join(search_terms)

        # Always build text search conditions (for smiles, supplier_ref, comments, aliases)
        text_conditions = self._build_text_conditions(search_query)

        # Try to extract a reg_number from the query using multiple patterns
        reg_number_str = self._extract_reg_number(search_query)

        # Combine conditions: reg_number matches OR text matches
        if reg_number_str:
            reg_number_conditions = self._build_reg_number_conditions(reg_number_str)
            combined = reg_number_conditions | text_conditions
        else:
            combined = text_conditions

        return queryset.filter(combined)

    def _extract_reg_number(self, query):
        """
        Extract a registration number from various query formats.

        Handles multiple patterns similar to resolve_compound():
        - Malformed: NCL000-27421, NCL00-187
        - Standard: NCL-00018710, ncl-18710, NCL 18710, NCL_18710
        - Plain numeric: 18710
        - Embedded: "sample NCL-26042 test"

        Returns the extracted reg_number string (digits only), or None.
        """
        # Pattern 1: Malformed NCL - dash after zeros (NCL000-27421)
        # Must check first to avoid incorrect partial matches
        ncl_malformed = re.match(r'^NCL0+-(\d+)$', query, re.IGNORECASE)
        if ncl_malformed:
            return ncl_malformed.group(1)

        # Pattern 2: Standard NCL format - NCL-00018710, ncl-18710, NCL 18710, NCL_18710
        ncl_standard = re.match(r'^NCL[-_\s]?0*(\d+)$', query, re.IGNORECASE)
        if ncl_standard:
            return ncl_standard.group(1)

        # Pattern 3: Plain numeric
        if query.isdigit():
            return query

        # Pattern 4: Embedded malformed NCL pattern in longer string
        ncl_malformed_search = re.search(r'NCL0+-(\d+)', query, re.IGNORECASE)
        if ncl_malformed_search:
            return ncl_malformed_search.group(1)

        # Pattern 5: Embedded standard NCL pattern (handles NCL_26042, "text NCL-26042 text")
        ncl_embedded = re.search(r'NCL[-_\s]?0*(\d+)', query, re.IGNORECASE)
        if ncl_embedded:
            return ncl_embedded.group(1)

        return None

    def _build_reg_number_conditions(self, reg_number_str):
        """
        Build Q conditions for registration number matching.

        Supports both exact matches and prefix matches:
        - "18710" matches compound 18710 exactly
        - "187" matches compounds 187, 1870, 1871, ..., 18710, etc.
        """
        try:
            reg_number = int(reg_number_str)
        except ValueError:
            return Q(pk__isnull=True)  # Return empty Q that matches nothing

        # Start with exact match
        conditions = Q(reg_number=reg_number)

        # Add prefix matching - find all compounds where reg_number starts with these digits
        # Numbers starting with "187" are: 187, 1870-1879, 18700-18799, 187000-187999, etc.
        # Generate range conditions for each possible length
        for extra_digits in range(1, 9):  # Up to 8 additional digits
            multiplier = 10 ** extra_digits
            lower = reg_number * multiplier
            upper = (reg_number + 1) * multiplier
            conditions |= Q(reg_number__gte=lower, reg_number__lt=upper)

        return conditions

    def _build_text_conditions(self, search_query):
        """
        Build Q conditions for text field searching.

        Searches: smiles, supplier_ref, comments, aliases (JSON array)
        """
        conditions = Q()

        # Standard text fields
        conditions |= Q(smiles__icontains=search_query)
        conditions |= Q(supplier_ref__icontains=search_query)
        conditions |= Q(comments__icontains=search_query)

        # Aliases is a JSONField containing a list - use __icontains which works
        # on the JSON string representation
        conditions |= Q(aliases__icontains=search_query)

        return conditions
