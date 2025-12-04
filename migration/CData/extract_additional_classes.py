#!/usr/bin/env python3
"""
Extract class definitions from CCP4ModelData.py and CCP4CootData.py.

This script scans Python source files and extracts class metadata in the format
needed for cdata.json.
"""

import ast
import json
import re
from pathlib import Path
from typing import Dict, Any, List


class CCP4ClassExtractor:
    """Extract class metadata from CCP4i2 Python source files."""

    def __init__(self, source_file: Path):
        """
        Initialize extractor.

        Args:
            source_file: Path to Python source file to extract from
        """
        self.source_file = source_file
        self.source_text = source_file.read_text()

    def extract_all_classes(self) -> Dict[str, Dict[str, Any]]:
        """
        Extract all class definitions from source file.

        Returns:
            Dictionary mapping class names to their metadata
        """
        classes = {}

        # Parse the Python source
        try:
            tree = ast.parse(self.source_text)
        except SyntaxError as e:
            print(f"Error parsing {self.source_file}: {e}")
            return classes

        # Find all class definitions
        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                class_data = self._extract_class(node)
                if class_data:
                    classes[node.name] = class_data

        return classes

    def _extract_class(self, node: ast.ClassDef) -> Dict[str, Any]:
        """
        Extract metadata from a single class definition.

        Args:
            node: AST ClassDef node

        Returns:
            Dictionary of class metadata
        """
        # Get class source code
        class_source = ast.get_source_segment(self.source_text, node)
        if not class_source:
            return None

        # Extract base classes
        base_classes = []
        immediate_parent = None
        for base in node.bases:
            if isinstance(base, ast.Name):
                base_name = base.id
                base_classes.append(base_name)
                if not immediate_parent:
                    immediate_parent = base_name
            elif isinstance(base, ast.Attribute):
                # Handle CCP4File.CDataFile style
                base_name = ast.unparse(base)
                base_classes.append(base_name)
                if not immediate_parent:
                    immediate_parent = base_name

        # Extract docstring
        docstring = ast.get_docstring(node) or ""

        # Extract class-level attributes (CONTENTS, QUALIFIERS, etc.)
        contents = {}
        qualifiers = {}
        qualifiers_def = {}
        error_codes = {}
        subtypes = {}
        content_flags = {}

        for item in node.body:
            if isinstance(item, ast.Assign):
                for target in item.targets:
                    if isinstance(target, ast.Name):
                        attr_name = target.name

                        # Try to evaluate the value
                        try:
                            if attr_name == 'CONTENTS':
                                contents = self._extract_dict_assignment(class_source, 'CONTENTS')
                            elif attr_name == 'QUALIFIERS':
                                qualifiers = self._extract_dict_assignment(class_source, 'QUALIFIERS')
                            elif attr_name == 'QUALIFIERS_DEFINITION':
                                qualifiers_def = self._extract_dict_assignment(class_source, 'QUALIFIERS_DEFINITION')
                            elif attr_name == 'ERROR_CODES':
                                error_codes = self._extract_dict_assignment(class_source, 'ERROR_CODES')
                            elif attr_name.startswith('SUBTYPE_'):
                                # Extract subtype constants
                                value = self._extract_constant_value(class_source, attr_name)
                                if value is not None:
                                    subtypes[attr_name] = value
                            elif attr_name.startswith('CONTENT_FLAG_'):
                                # Extract content flag constants
                                value = self._extract_constant_value(class_source, attr_name)
                                if value is not None:
                                    content_flags[attr_name] = value
                        except Exception as e:
                            print(f"  Warning: Could not extract {attr_name}: {e}")

        # Build metadata structure
        metadata = {
            'module': self.source_file.stem,  # e.g., 'CCP4ModelData'
            'class': node.name,
            'file_path': str(self.source_file),
            'docstring': docstring,
            'base_classes': base_classes,
            'immediate_parent': immediate_parent,
            'mro': [],  # Would need runtime info
            'CONTENTS': contents,
            'CONTENTS_ORDER': list(contents.keys()),
            'QUALIFIERS': qualifiers,
            'QUALIFIERS_ORDER': list(qualifiers.keys()),
            'QUALIFIERS_DEFINITION': qualifiers_def,
            'ERROR_CODES': error_codes
        }

        # Add SUBTYPES/CONTENT_FLAGS if present
        if subtypes or content_flags:
            # Will add these in a separate step using add_metadata script
            pass

        return metadata

    def _extract_dict_assignment(self, source: str, var_name: str) -> Dict:
        """
        Extract dictionary assignment from source code.

        Args:
            source: Class source code
            var_name: Variable name to extract (e.g., 'CONTENTS')

        Returns:
            Extracted dictionary (simplified)
        """
        # This is a simplified extraction - for full accuracy would need
        # runtime evaluation or more sophisticated parsing

        # Look for pattern like: VARNAME = {...}
        pattern = rf'{var_name}\s*=\s*\{{([^}}]*)\}}'
        match = re.search(pattern, source, re.DOTALL)

        if match:
            # For now, return empty dict - would need more sophisticated parsing
            # to handle nested dicts, class references, etc.
            return {}

        # Look for update pattern: VARNAME.update(...)
        pattern = rf'{var_name}\.update\([^)]+\)'
        if re.search(pattern, source):
            return {}

        return {}

    def _extract_constant_value(self, source: str, const_name: str) -> Any:
        """
        Extract value of a class constant.

        Args:
            source: Class source code
            const_name: Constant name

        Returns:
            Constant value or None
        """
        # Look for pattern like: CONST_NAME = value
        pattern = rf'^\s*{const_name}\s*=\s*(\d+)'
        match = re.search(pattern, source, re.MULTILINE)

        if match:
            return int(match.group(1))

        return None


def main():
    """Main extraction routine."""

    # Paths
    ccp4i2_core = Path('/Users/nmemn/Developer/ccp4i2/core')
    cdata_json_path = Path(__file__).parent / 'cdata.json'

    # Files to extract
    files_to_extract = [
        ccp4i2_core / 'CCP4ModelData.py',
        ccp4i2_core / 'CCP4CootData.py'
    ]

    print("="*70)
    print("CCP4 Class Metadata Extractor")
    print("="*70)

    all_extracted = {}

    for source_file in files_to_extract:
        if not source_file.exists():
            print(f"\n⚠ Warning: {source_file} not found, skipping")
            continue

        print(f"\nExtracting from {source_file.name}...")

        extractor = CCP4ClassExtractor(source_file)
        classes = extractor.extract_all_classes()

        print(f"  Found {len(classes)} classes")

        # Filter to only CDataFile descendants (heuristic)
        datafile_classes = {}
        for name, data in classes.items():
            # Check if inherits from CDataFile or similar
            if any('DataFile' in str(base) or 'CList' in str(base)
                   for base in data.get('base_classes', [])):
                datafile_classes[name] = data
                print(f"    • {name} (parent: {data.get('immediate_parent', '?')})")

        all_extracted.update(datafile_classes)

    print(f"\n{'='*70}")
    print(f"Extraction Summary")
    print(f"{'='*70}")
    print(f"Total classes extracted: {len(all_extracted)}")
    print(f"\nClasses: {', '.join(sorted(all_extracted.keys()))}")

    # Note: This is a simplified extractor
    # For production use, would need to handle:
    # - Full CONTENTS/QUALIFIERS parsing
    # - Proper type resolution
    # - MRO calculation
    # - Runtime evaluation of complex expressions

    print(f"\n⚠ Note: This is a basic extraction showing class structure.")
    print(f"For full metadata, manual review of source files is recommended.")
    print(f"\nFor CPdbDataFile and CCootHistoryDataFile, suggest adding minimal")
    print(f"entries manually with just the CONTENT_FLAGS/SUBTYPES metadata.")


if __name__ == '__main__':
    main()
