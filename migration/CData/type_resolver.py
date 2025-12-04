"""
Type resolution system for CData class generation.

Resolves all type references and determines proper imports.
"""

import re
from pathlib import Path
from typing import Dict, Set, Tuple, Optional


class TypeResolver:
    """Resolves type references and generates proper import statements."""

    def __init__(self, classes_data: Dict[str, dict]):
        """
        Initialize with classes data from cdata.json.

        Args:
            classes_data: Dict mapping class names to their metadata
        """
        self.classes_data = classes_data

        # Fundamental types from core.base_object.fundamental_types
        self.fundamental_types = {
            'CInt', 'CFloat', 'CString', 'CBoolean', 'CList'
        }

        # Type aliases - NOTE: This is now empty because all custom types
        # (COneWord, CUUID, CProjectId, etc.) are proper classes in cdata.json
        # with their own stubs/implementations. They're not aliases anymore.
        self.type_aliases = {}

        # Base classes from base_classes.py (only the hand-written ones)
        self.base_classes = {
            'CData', 'CContainer', 'CDataFile', 'CDataFileContent'
        }

        # Build class-to-file mapping
        self.class_to_file = self._build_class_to_file_map()

    def _build_class_to_file_map(self) -> Dict[str, str]:
        """Build mapping of class name to output file name."""
        mapping = {}
        for class_name, class_info in self.classes_data.items():
            file_path = class_info.get('file_path', 'unknown.py')
            basename = Path(file_path).name
            mapping[class_name] = basename
        return mapping

    def resolve_type_name(self, type_str: str) -> str:
        """
        Extract clean type name from various formats.

        Args:
            type_str: Type string like "<class 'CCP4ModelData.CAsuContent'>" or "CString"

        Returns:
            Clean type name like "CAsuContent" or "CString"
        """
        if not isinstance(type_str, str):
            return 'Any'

        # Handle "<class 'module.ClassName'>" format
        if type_str.startswith("<class '") and type_str.endswith("'>"):
            type_str = type_str[8:-2]

        # Handle fully qualified names like "ccp4x.data_scan.CCP4ModelData.CAsuContent"
        # or "core.CCP4Data.CString"
        if '.' in type_str:
            parts = type_str.split('.')
            type_str = parts[-1]  # Take last part as class name

        return type_str

    def get_type_location(self, type_name: str) -> Tuple[str, str]:
        """
        Determine where a type is defined.

        Args:
            type_name: Clean type name like "CString" or "CAsuContent"

        Returns:
            Tuple of (module_path, category) where category is:
            - 'fundamental': core.base_object.fundamental_types
            - 'base': core.base_object.base_classes
            - 'custom': core.generated.<filename>
            - 'unknown': Type not found
        """
        # Check fundamental types
        if type_name in self.fundamental_types:
            return 'core.base_object.fundamental_types', 'fundamental'

        # Check type aliases
        if type_name in self.type_aliases:
            return 'core.base_object.fundamental_types', 'fundamental'

        # Check base classes
        if type_name in self.base_classes:
            return 'core.base_object.base_classes', 'base'

        # Check custom classes
        if type_name in self.class_to_file:
            filename = self.class_to_file[type_name]
            module = f"core.generated.{filename[:-3]}"  # Remove .py
            return module, 'custom'

        return 'unknown', 'unknown'

    def get_imports_for_file(self, filename: str, classes_in_file: list) -> str:
        """
        Generate all import statements needed for a generated file.

        Args:
            filename: Output filename like "CCP4ModelData.py"
            classes_in_file: List of (class_name, class_data) tuples

        Returns:
            String containing all import statements
        """
        imports = []

        # Standard imports
        imports.append('"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.')
        imports.append('')
        imports.append('To extend these classes, create subclasses in core/extensions/')
        imports.append('"""')
        imports.append('')
        imports.append('from __future__ import annotations')
        imports.append('from typing import TYPE_CHECKING, Optional, Any')
        imports.append('')

        # Track what we need to import
        needed_fundamentals = set()
        needed_base_classes = set()
        needed_custom_classes = {}  # module -> set of classes

        # Analyze all classes in this file
        for class_name, class_data in classes_in_file:
            # Check parent class
            parent = class_data.get('immediate_parent', 'CData')
            if parent:
                module, category = self.get_type_location(parent)
                if category == 'fundamental':
                    needed_fundamentals.add(parent)
                elif category == 'base':
                    needed_base_classes.add(parent)
                elif category == 'custom' and module != f"core.generated.{filename[:-3]}":
                    if module not in needed_custom_classes:
                        needed_custom_classes[module] = set()
                    needed_custom_classes[module].add(parent)

            # Check attributes
            for attr_name, attr_info in class_data.get('CONTENTS', {}).items():
                attr_type_str = attr_info.get('class', '')
                attr_type = self.resolve_type_name(attr_type_str)

                module, category = self.get_type_location(attr_type)
                if category == 'fundamental':
                    needed_fundamentals.add(attr_type)
                elif category == 'base':
                    needed_base_classes.add(attr_type)
                elif category == 'custom' and module != f"core.generated.{filename[:-3]}":
                    if module not in needed_custom_classes:
                        needed_custom_classes[module] = set()
                    needed_custom_classes[module].add(attr_type)

        # Always import base decorator and metadata classes
        imports.append('# Metadata system')
        imports.append('from core.base_object.class_metadata import cdata_class, attribute, AttributeType')
        imports.append('')

        # Import base classes
        if needed_base_classes:
            base_imports = ', '.join(sorted(needed_base_classes))
            imports.append('# Base classes')
            imports.append(f'from core.base_object.base_classes import {base_imports}')
            imports.append('')

        # Import fundamental types
        if needed_fundamentals:
            fund_imports = ', '.join(sorted(needed_fundamentals))
            imports.append('# Fundamental types')
            imports.append(f'from core.base_object.fundamental_types import {fund_imports}')
            imports.append('')

        # Import custom classes from other files
        if needed_custom_classes:
            imports.append('# Cross-file class references')
            for module in sorted(needed_custom_classes.keys()):
                classes = ', '.join(sorted(needed_custom_classes[module]))
                imports.append(f'from {module} import {classes}')
            imports.append('')

        return '\n'.join(imports)

    def resolve_attribute_type_for_code(self, attr_type_str: str) -> str:
        """
        Resolve attribute type for use in generated code.

        Args:
            attr_type_str: Type string from CONTENTS

        Returns:
            Python type annotation string like "CString" or "CAsuContent"
        """
        type_name = self.resolve_type_name(attr_type_str)
        module, category = self.get_type_location(type_name)

        if category == 'unknown':
            return 'Any'

        return type_name

    def is_fundamental_or_alias(self, type_name: str) -> bool:
        """Check if type is a fundamental type or alias."""
        clean_name = self.resolve_type_name(type_name)
        return clean_name in self.fundamental_types or clean_name in self.type_aliases


def extract_type_aliases_from_source(fundamental_types_path: Path) -> Dict[str, str]:
    """
    Extract type aliases from fundamental_types.py source code.

    Args:
        fundamental_types_path: Path to fundamental_types.py

    Returns:
        Dict mapping alias name to base type (e.g., {'CUUID': 'CString'})
    """
    aliases = {}

    if not fundamental_types_path.exists():
        return aliases

    content = fundamental_types_path.read_text()

    # Pattern: ALIAS = CType (where ALIAS starts with C and is all caps or PascalCase)
    pattern = r'^([A-Z][A-Za-z0-9]*)\s*=\s*(C[A-Za-z0-9]+)\s*(?:#.*)?$'

    for line in content.split('\n'):
        line = line.strip()
        match = re.match(pattern, line)
        if match:
            alias_name, base_type = match.groups()
            aliases[alias_name] = base_type

    return aliases
