"""
Crystal selector - parse and resolve crystal specifiers.

Supports:
- Individual crystals: x0001, x0002
- Ranges: x0001-x0010, x0001:x0010
- By compound ID: NCL-00031069
- Wildcards: * (all), x000* (glob pattern)
- Exclusions: !x0005, !x0003-x0007
"""

import os
import re
import glob as globmodule
from typing import List, Set, Optional, Tuple


# Pattern to match crystal directory names like "CDK4CyclinD1-x0001"
CRYSTAL_DIR_PATTERN = re.compile(r'^(.+)-x(\d+)$')

# Pattern to match crystal specifiers like "x0001" or "0001"
CRYSTAL_SPEC_PATTERN = re.compile(r'^x?(\d+)$')

# Pattern to match range specifiers like "x0001-x0010" or "x0001:x0010"
RANGE_PATTERN = re.compile(r'^x?(\d+)[-:]x?(\d+)$')

# Pattern to match compound IDs like "NCL-00031069"
COMPOUND_PATTERN = re.compile(r'^[A-Z]{2,4}-\d+$')


class CrystalSelector:
    """
    Discovers and selects crystal directories based on various specifiers.
    """

    def __init__(self, base_path, crystal_prefix=None):
        # type: (str, Optional[str]) -> None
        """
        Initialize the selector.

        Args:
            base_path: Path to the model_building directory (or similar)
            crystal_prefix: Expected prefix for crystal dirs (e.g., "CDK4CyclinD1")
                           If None, will auto-detect from first crystal found.
        """
        self.base_path = os.path.abspath(base_path)
        self.crystal_prefix = crystal_prefix
        self._crystal_cache = None  # type: Optional[dict]
        self._compound_cache = None  # type: Optional[dict]

    def _scan_crystals(self):
        # type: () -> dict
        """Scan the base path for crystal directories."""
        if self._crystal_cache is not None:
            return self._crystal_cache

        crystals = {}  # index -> (dir_name, full_path)

        for entry in os.listdir(self.base_path):
            full_path = os.path.join(self.base_path, entry)
            if not os.path.isdir(full_path):
                continue

            match = CRYSTAL_DIR_PATTERN.match(entry)
            if match:
                prefix, index_str = match.groups()

                # Auto-detect prefix from first crystal
                if self.crystal_prefix is None:
                    self.crystal_prefix = prefix
                elif prefix != self.crystal_prefix:
                    # Skip directories with different prefix
                    continue

                index = int(index_str)
                crystals[index] = (entry, full_path)

        self._crystal_cache = crystals
        return crystals

    def _scan_compounds(self):
        # type: () -> dict
        """Build a mapping from compound ID to crystal index."""
        if self._compound_cache is not None:
            return self._compound_cache

        crystals = self._scan_crystals()
        compounds = {}  # compound_id -> [indices]

        for index, (dir_name, full_path) in crystals.items():
            compound_dir = os.path.join(full_path, 'compound')
            if not os.path.isdir(compound_dir):
                continue

            # Look for compound files (NCL-*.cif, NCL-*.pdb, etc.)
            for entry in os.listdir(compound_dir):
                # Extract compound ID from filename
                name_part = os.path.splitext(entry)[0]
                # Handle _TMP suffix directories
                if name_part.endswith('_TMP'):
                    name_part = name_part[:-4]

                if COMPOUND_PATTERN.match(name_part):
                    if name_part not in compounds:
                        compounds[name_part] = []
                    if index not in compounds[name_part]:
                        compounds[name_part].append(index)

        self._compound_cache = compounds
        return compounds

    def get_all_indices(self):
        # type: () -> List[int]
        """Get all available crystal indices, sorted."""
        crystals = self._scan_crystals()
        return sorted(crystals.keys())

    def get_index_range(self):
        # type: () -> Tuple[int, int]
        """Get the min and max crystal indices."""
        indices = self.get_all_indices()
        if not indices:
            return (0, 0)
        return (min(indices), max(indices))

    def get_crystal_path(self, index):
        # type: (int) -> Optional[str]
        """Get the full path for a crystal by index."""
        crystals = self._scan_crystals()
        if index in crystals:
            return crystals[index][1]
        return None

    def get_crystal_dir_name(self, index):
        # type: (int) -> Optional[str]
        """Get the directory name for a crystal by index."""
        crystals = self._scan_crystals()
        if index in crystals:
            return crystals[index][0]
        return None

    def parse_specifier(self, spec):
        # type: (str) -> Set[int]
        """
        Parse a single crystal specifier and return matching indices.

        Args:
            spec: A specifier like "x0001", "x0001-x0010", "NCL-00031069", "*"

        Returns:
            Set of crystal indices matching the specifier
        """
        spec = spec.strip()

        # Check for exclusion prefix
        if spec.startswith('!'):
            # This should be handled at a higher level
            raise ValueError("Exclusion specs should be handled separately")

        # Wildcard - all crystals
        if spec == '*' or spec == 'all':
            return set(self.get_all_indices())

        # Range specifier: x0001-x0010 or x0001:x0010
        range_match = RANGE_PATTERN.match(spec)
        if range_match:
            start, end = int(range_match.group(1)), int(range_match.group(2))
            if start > end:
                start, end = end, start
            crystals = self._scan_crystals()
            return set(i for i in range(start, end + 1) if i in crystals)

        # Single crystal: x0001 or 0001
        crystal_match = CRYSTAL_SPEC_PATTERN.match(spec)
        if crystal_match:
            index = int(crystal_match.group(1))
            crystals = self._scan_crystals()
            if index in crystals:
                return {index}
            return set()

        # Compound ID: NCL-00031069
        if COMPOUND_PATTERN.match(spec):
            compounds = self._scan_compounds()
            if spec in compounds:
                return set(compounds[spec])
            return set()

        # Glob pattern: x000*
        if '*' in spec or '?' in spec:
            # Normalize pattern to match directory names
            if spec.startswith('x'):
                pattern = '{}-{}'.format(self.crystal_prefix or '*', spec)
            else:
                pattern = spec

            full_pattern = os.path.join(self.base_path, pattern)
            matches = globmodule.glob(full_pattern)

            result = set()
            crystals = self._scan_crystals()
            for match_path in matches:
                dir_name = os.path.basename(match_path)
                dir_match = CRYSTAL_DIR_PATTERN.match(dir_name)
                if dir_match:
                    index = int(dir_match.group(2))
                    if index in crystals:
                        result.add(index)
            return result

        # Unknown specifier
        raise ValueError("Cannot parse crystal specifier: {}".format(spec))

    def select(self, specifiers):
        # type: (List[str]) -> List[int]
        """
        Select crystals based on a list of specifiers.

        Specifiers starting with '!' are exclusions.

        Args:
            specifiers: List of specifier strings

        Returns:
            Sorted list of selected crystal indices
        """
        included = set()  # type: Set[int]
        excluded = set()  # type: Set[int]

        for spec in specifiers:
            spec = spec.strip()
            if not spec:
                continue

            if spec.startswith('!'):
                # Exclusion
                excluded.update(self.parse_specifier(spec[1:]))
            else:
                # Inclusion
                included.update(self.parse_specifier(spec))

        # If no explicit inclusions, include all
        if not included:
            included = set(self.get_all_indices())

        result = included - excluded
        return sorted(result)

    def summary(self):
        # type: () -> str
        """Return a summary of available crystals."""
        crystals = self._scan_crystals()
        compounds = self._scan_compounds()

        lines = [
            "Crystal directory: {}".format(self.base_path),
            "Prefix: {}".format(self.crystal_prefix or "(auto-detect)"),
            "Total crystals: {}".format(len(crystals)),
        ]

        if crystals:
            indices = sorted(crystals.keys())
            lines.append("Index range: x{:04d} - x{:04d}".format(
                min(indices), max(indices)))

            # Show gaps if any
            expected = set(range(min(indices), max(indices) + 1))
            missing = expected - set(indices)
            if missing:
                lines.append("Missing indices: {}".format(len(missing)))

        lines.append("Unique compounds: {}".format(len(compounds)))

        return '\n'.join(lines)
