"""
Stub Code Generator for CData Classes

Generates stub classes (with "Stub" suffix) that will be extended by
implementation classes in the core/ directory.
"""

import json
import subprocess
import argparse
from pathlib import Path
from typing import Dict, Any, List, Tuple

from type_resolver import TypeResolver
from class_graph import ClassDependencyGraph


class StubTypeResolver(TypeResolver):
    """Type resolver that uses cdata_stubs instead of generated."""

    def get_type_location(self, type_name: str):
        """Override to use core.cdata_stubs path."""
        # Check fundamental types
        if type_name in self.fundamental_types:
            return 'core.base_object.fundamental_types', 'fundamental'

        # Check type aliases
        if type_name in self.type_aliases:
            return 'core.base_object.fundamental_types', 'fundamental'

        # Check base classes
        if type_name in self.base_classes:
            return 'core.base_object.base_classes', 'base'

        # Check custom classes (use cdata_stubs instead of generated)
        if type_name in self.class_to_file:
            filename = self.class_to_file[type_name]
            module = f"core.cdata_stubs.{filename[:-3]}"  # Remove .py, use cdata_stubs
            return module, 'custom'

        return 'unknown', 'unknown'

    def get_imports_for_file(self, filename: str, classes_in_file: list) -> str:
        """Generate all import statements needed for a stub file."""
        needed_fundamental = set()
        needed_base_classes = set()
        needed_custom_classes = {}  # {module: {class1, class2, ...}}

        # Analyze each class in this file
        for class_name, class_data in classes_in_file:
            # Check parent class
            parent = class_data.get('immediate_parent', 'CData')
            if parent:
                module, category = self.get_type_location(parent)
                if category == 'fundamental':
                    needed_fundamental.add(parent)
                elif category == 'base':
                    needed_base_classes.add(parent)
                elif category == 'custom' and module != f"core.cdata_stubs.{filename[:-3]}":
                    if module not in needed_custom_classes:
                        needed_custom_classes[module] = set()
                    needed_custom_classes[module].add(parent)

            # Check attribute types
            attributes = class_data.get('CONTENTS', {})
            for attr_name, attr_info in attributes.items():
                attr_type_str = attr_info.get('class', 'str')
                attr_type = self.resolve_type_name(attr_type_str)

                module, category = self.get_type_location(attr_type)
                if category == 'fundamental':
                    needed_fundamental.add(attr_type)
                elif category == 'base':
                    needed_base_classes.add(attr_type)
                elif category == 'custom' and module != f"core.cdata_stubs.{filename[:-3]}":
                    if module not in needed_custom_classes:
                        needed_custom_classes[module] = set()
                    needed_custom_classes[module].add(attr_type)

        # Build import statements
        lines = ['"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.']
        lines.append('')
        lines.append('This is a stub file - extend classes in core/ to add methods.')
        lines.append('"""')
        lines.append('')
        lines.append('from __future__ import annotations')
        lines.append('from typing import TYPE_CHECKING, Optional, Any')
        lines.append('')

        # Metadata system
        lines.append('# Metadata system')
        lines.append('from core.base_object.class_metadata import cdata_class, attribute, AttributeType')
        lines.append('')

        # Base classes
        if needed_base_classes:
            lines.append('# Base classes')
            base_imports = ', '.join(sorted(needed_base_classes))
            lines.append(f'from core.base_object.base_classes import {base_imports}')
            lines.append('')

        # Fundamental types
        if needed_fundamental:
            lines.append('# Fundamental types')
            fund_imports = ', '.join(sorted(needed_fundamental))
            lines.append(f'from core.base_object.fundamental_types import {fund_imports}')
            lines.append('')

        # Custom classes from other stub files
        if needed_custom_classes:
            lines.append('# Cross-file stub class references')
            for module in sorted(needed_custom_classes.keys()):
                classes = sorted(needed_custom_classes[module])
                # Add Stub suffix to imported classes
                stub_classes = [f"{c}Stub" for c in classes]
                class_list = ', '.join(stub_classes)
                lines.append(f'from {module} import {class_list}')
            lines.append('')

        return '\n'.join(lines)


class StubCodeGenerator:
    """Generates stub CData classes with 'Stub' suffix."""

    def __init__(self, cdata_json_path: Path):
        """
        Initialize generator.

        Args:
            cdata_json_path: Path to cdata.json metadata file
        """
        self.json_path = cdata_json_path
        self.data = self._load_json()
        all_classes = self.data.get('classes', {})

        # Base classes that don't get Stub suffix (hand-written in base_object/)
        self.base_classes = {
            'CData', 'CContainer', 'CDataFile', 'CDataFileContent',
            'CInt', 'CFloat', 'CString', 'CBoolean', 'CList',
        }

        # Filter out base classes
        self.classes = {
            name: data for name, data in all_classes.items()
            if name not in self.base_classes
        }

        print(f"Loaded {len(all_classes)} classes, generating {len(self.classes)} stubs " +
              f"(excluding {len(self.base_classes)} hand-written base classes)")

        # Initialize type resolver and dependency graph
        print("Initializing stub type resolver...")
        self.type_resolver = StubTypeResolver(all_classes)

        print("Building dependency graph...")
        self.dep_graph = ClassDependencyGraph(self.classes, self.type_resolver)

    def _load_json(self) -> Dict[str, Any]:
        """Load and parse cdata.json."""
        print(f"Loading metadata from {self.json_path}...")
        with open(self.json_path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def add_stub_suffix(self, class_name: str) -> str:
        """
        Add 'Stub' suffix to class name if not a base class.

        Args:
            class_name: Original class name

        Returns:
            Class name with Stub suffix if applicable
        """
        if class_name in self.base_classes:
            return class_name
        return f"{class_name}Stub"

    def render_decorator(self, class_data: dict) -> List[str]:
        """
        Render @cdata_class decorator with all metadata.

        Args:
            class_data: Class metadata dictionary

        Returns:
            List of lines for the decorator
        """
        lines = ['@cdata_class(']

        # Render attributes if present
        attributes = class_data.get('CONTENTS', {})
        if attributes:
            lines.append('    attributes={')
            for attr_name, attr_info in attributes.items():
                attr_type_str = attr_info.get('class', 'str')
                attr_type = self.type_resolver.resolve_type_name(attr_type_str)

                # Map to AttributeType
                type_map = {
                    'CInt': 'AttributeType.INT',
                    'CFloat': 'AttributeType.FLOAT',
                    'CString': 'AttributeType.STRING',
                    'CBoolean': 'AttributeType.BOOLEAN',
                }

                if attr_type in type_map:
                    attr_type_enum = type_map[attr_type]
                    lines.append(f'        "{attr_name}": attribute({attr_type_enum}),')
                else:
                    # For custom classes, add Stub suffix if not a base class
                    custom_class = self.add_stub_suffix(attr_type)
                    lines.append(f'        "{attr_name}": attribute(AttributeType.CUSTOM, custom_class="{custom_class}"),')

            lines.append('    },')

        # Render error codes
        error_codes = class_data.get('ERROR_CODES', {})
        if error_codes:
            lines.append(f'    error_codes={json.dumps(error_codes, indent=8)},')

        # Render qualifiers (filter out "NotImplemented" values)
        qualifiers = class_data.get('QUALIFIERS', {})
        if qualifiers:
            filtered_qualifiers = {k: v for k, v in qualifiers.items() if v != "NotImplemented"}
            if filtered_qualifiers:
                lines.append('    qualifiers={')
                for k, v in filtered_qualifiers.items():
                    lines.append(f'        "{k}": {repr(v)},')
                lines.append('    },')

        # Render qualifiers_order
        qualifiers_order = class_data.get('QUALIFIERS_ORDER', [])
        if qualifiers_order:
            lines.append(f'    qualifiers_order={repr(qualifiers_order)},')

        # Render qualifiers_definition
        qualifiers_definition = class_data.get('QUALIFIERS_DEFINITION', {})
        if qualifiers_definition:
            lines.append('    qualifiers_definition={')
            for qname, qdef in qualifiers_definition.items():
                if isinstance(qdef, dict):
                    qdef_copy = qdef.copy()
                    if 'type' in qdef_copy:
                        type_val = qdef_copy['type']
                        if isinstance(type_val, str) and type_val.startswith("<class '"):
                            type_val = type_val[8:-2]
                        type_map = {'str': 'str', 'int': 'int', 'float': 'float', 'bool': 'bool', 'dict': 'dict', 'list': 'list'}
                        qdef_copy['type'] = type_map.get(type_val, f'"{type_val}"')

                    lines.append(f'        "{qname}": {qdef_copy},')
                else:
                    lines.append(f'        "{qname}": {repr(qdef)},')
            lines.append('    },')

        # Render contents_order
        contents_order = class_data.get('CONTENTS_ORDER', [])
        if contents_order:
            lines.append(f'    contents_order={repr(contents_order)},')

        # Render content_qualifiers (per-field qualifiers from CONTENTS)
        # This captures qualifiers like allowUndefined=false, mustExist=true for specific fields
        content_qualifiers = {}
        for attr_name, attr_info in attributes.items():
            field_qualifiers = attr_info.get('qualifiers', {})
            if field_qualifiers:
                content_qualifiers[attr_name] = field_qualifiers
        if content_qualifiers:
            lines.append('    content_qualifiers={')
            for field_name, field_quals in content_qualifiers.items():
                lines.append(f'        "{field_name}": {repr(field_quals)},')
            lines.append('    },')

        # Render gui_label
        gui_label = class_data.get('gui_label', '')
        if gui_label:
            lines.append(f'    gui_label="{gui_label}",')

        lines.append(')')
        return lines

    def render_class(self, class_name: str, class_data: dict) -> List[str]:
        """
        Render complete stub class definition.

        Args:
            class_name: Name of the class (without Stub suffix)
            class_data: Class metadata dictionary

        Returns:
            List of lines for the complete class
        """
        lines = []

        # Add decorator
        lines.extend(self.render_decorator(class_data))

        # Class definition with Stub suffix
        stub_name = self.add_stub_suffix(class_name)
        parent = class_data.get('immediate_parent', 'CData')
        parent_stub = self.add_stub_suffix(parent)

        lines.append(f'class {stub_name}({parent_stub}):')

        # Docstring
        docstring = class_data.get('docstring', f'Auto-generated stub for {class_name}.')
        lines.append(f'    """')
        lines.append(f'    {docstring}')
        lines.append(f'    ')
        lines.append(f'    This is a pure data class stub. Extend it in core/{class_name}.py')
        lines.append(f'    to add methods and implementation-specific functionality.')
        lines.append(f'    """')
        lines.append('')

        # Add MTZ metadata class constants (for backward compatibility)
        self._render_mtz_metadata(class_data, lines)

        # Type-annotated attributes
        attributes = class_data.get('CONTENTS', {})
        if attributes:
            for attr_name, attr_info in attributes.items():
                attr_type_str = attr_info.get('class', 'Any')
                attr_type = self.type_resolver.resolve_attribute_type_for_code(attr_type_str)
                # Add Stub suffix to custom types
                if attr_type not in ['int', 'float', 'str', 'bool', 'Any', 'Dict', 'List']:
                    attr_type = self.add_stub_suffix(attr_type)
                lines.append(f'    {attr_name}: Optional[{attr_type}] = None')
            lines.append('')

        # __init__ method
        lines.append('    def __init__(self, parent=None, name=None, **kwargs):')
        lines.append('        """')
        lines.append(f'        Initialize {stub_name}.')
        lines.append('')
        lines.append('        Args:')
        lines.append('            parent: Parent object in hierarchy')
        lines.append('            name: Object name')
        lines.append('            **kwargs: Additional keyword arguments')
        lines.append('        """')
        lines.append('        super().__init__(parent=parent, name=name, **kwargs)')

        # Add blank line after class
        lines.append('')
        lines.append('')

        return lines

    def _render_mtz_metadata(self, class_data: dict, lines: List[str]) -> None:
        """
        Render MTZ-specific class constants for CMiniMtzDataFile subclasses.

        Adds CONTENT_FLAG_*, SUBTYPE_*, CONTENT_ANNOTATION, and CONTENT_SIGNATURE_LIST
        for backward compatibility with old CCP4i2 code.

        Args:
            class_data: Class metadata dictionary
            lines: List to append generated lines to
        """
        content_flags = class_data.get('CONTENT_FLAGS', {})
        subtypes = class_data.get('SUBTYPES', {})

        if not content_flags and not subtypes:
            return  # No MTZ metadata, skip

        # Add subtypes first (if present)
        if subtypes:
            lines.append('    # Subtype constants')
            for name, info in sorted(subtypes.items(), key=lambda x: x[1]['value']):
                value = info['value']
                description = info['description']
                lines.append(f'    {name} = {value}  # {description}')
            lines.append('')

        # Add content flags
        if content_flags:
            lines.append('    # Content flag constants')
            for name, info in sorted(content_flags.items(), key=lambda x: x[1]['value']):
                value = info['value']
                annotation = info['annotation']
                lines.append(f'    {name} = {value}  # {annotation}')
            lines.append('')

        # Build CONTENT_ANNOTATION list (indexed by contentFlag - 1)
        if content_flags:
            # Sort by value to ensure correct order
            sorted_flags = sorted(content_flags.items(), key=lambda x: x[1]['value'])
            annotations = [info['annotation'] for name, info in sorted_flags]
            lines.append('    # Content annotations (indexed by contentFlag - 1)')
            lines.append(f'    CONTENT_ANNOTATION = {repr(annotations)}')
            lines.append('')

        # Build CONTENT_SIGNATURE_LIST (indexed by contentFlag - 1)
        if content_flags:
            sorted_flags = sorted(content_flags.items(), key=lambda x: x[1]['value'])
            signatures = [info['columns'] for name, info in sorted_flags]
            lines.append('    # Column signatures for each content flag (indexed by contentFlag - 1)')
            lines.append(f'    CONTENT_SIGNATURE_LIST = {repr(signatures)}')
            lines.append('')

    def generate_file(self, filename: str, classes_in_file: List[Tuple[str, dict]]) -> str:
        """
        Generate complete file content for stubs.

        Args:
            filename: Output filename
            classes_in_file: List of (class_name, class_data) tuples

        Returns:
            Complete file content as string
        """
        lines = []

        # Generate imports
        imports = self.type_resolver.get_imports_for_file(filename, classes_in_file)
        lines.append(imports)

        # Generate each class
        for class_name, class_data in classes_in_file:
            class_lines = self.render_class(class_name, class_data)
            lines.extend(class_lines)

        return '\n'.join(lines)

    def generate_all(self, output_dir: Path, format_code: bool = True):
        """
        Generate all stub files.

        Args:
            output_dir: Directory to write generated files
            format_code: Whether to run autopep8 formatting
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get sorted classes by file
        print("\nGenerating dependency-sorted stub class definitions...")
        sorted_by_file = self.dep_graph.get_sorted_classes_by_file()

        print(f"\nGenerating {len(sorted_by_file)} stub files...")

        for filename, classes_in_file in sorted(sorted_by_file.items()):
            output_path = output_dir / filename
            print(f"  Writing {filename} ({len(classes_in_file)} classes)...")

            # Generate file content
            content = self.generate_file(filename, classes_in_file)

            # Write to file
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(content)

            # Format with autopep8
            if format_code:
                try:
                    subprocess.run(
                        ['autopep8', '--in-place', '--aggressive', str(output_path)],
                        check=True,
                        capture_output=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"    Warning: autopep8 formatting failed: {e}")
                except FileNotFoundError:
                    print("    Warning: autopep8 not found, skipping formatting")

        # Generate __init__.py
        self._generate_init_file(output_dir, sorted_by_file)

        print(f"\n✓ Generated {len(sorted_by_file)} stub files in {output_dir}")

    def _generate_init_file(self, output_dir: Path, sorted_by_file: Dict[str, List[Tuple[str, dict]]]):
        """Generate __init__.py that exports all stub classes."""
        init_path = output_dir / '__init__.py'

        lines = [
            '"""',
            'Auto-generated CData stub classes.',
            '',
            'These are pure data classes with metadata decorators.',
            'Extend them in core/ to add methods and implementation logic.',
            '',
            'DO NOT EDIT - Generated from CCP4i2 metadata.',
            '"""',
            '',
        ]

        # Import and re-export all stub classes
        for filename in sorted(sorted_by_file.keys()):
            module_name = filename[:-3]  # Remove .py
            classes_in_file = sorted_by_file[filename]

            if classes_in_file:
                # Export with Stub suffix
                stub_names = [self.add_stub_suffix(name) for name, _ in classes_in_file]
                class_list = ', '.join(stub_names)
                lines.append(f'from .{module_name} import {class_list}')

        lines.append('')
        lines.append('__all__ = [')
        for filename in sorted(sorted_by_file.keys()):
            for class_name, _ in sorted_by_file[filename]:
                stub_name = self.add_stub_suffix(class_name)
                lines.append(f'    "{stub_name}",')
        lines.append(']')
        lines.append('')

        with open(init_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))

        print(f"  Generated {init_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Generate stub CData classes from metadata'
    )
    parser.add_argument(
        '--input',
        type=Path,
        default=Path(__file__).parent / 'cdata.json',
        help='Path to cdata.json metadata file'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path(__file__).parent.parent.parent / 'core' / 'cdata_stubs',
        help='Output directory for generated stub files'
    )
    parser.add_argument(
        '--no-format',
        action='store_true',
        help='Skip autopep8 formatting'
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='Print dependency analysis report'
    )

    args = parser.parse_args()

    # Create generator
    generator = StubCodeGenerator(args.input)

    # Print dependency report if requested
    if args.report:
        generator.dep_graph.print_dependency_report()

    # Generate all stub files
    generator.generate_all(args.output, format_code=not args.no_format)

    return 0


if __name__ == '__main__':
    exit(main())
