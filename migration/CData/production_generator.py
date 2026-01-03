"""
Production Code Generator for CData Classes

Generates complete, production-ready Python classes from cdata.json metadata.
No manual patching required - output is directly usable.
"""

import json
import subprocess
import argparse
from pathlib import Path
from typing import Dict, Any, List, Tuple

from type_resolver import TypeResolver
from class_graph import ClassDependencyGraph


class ProductionCodeGenerator:
    """Generates complete CData class code from metadata."""

    def __init__(self, cdata_json_path: Path):
        """
        Initialize generator.

        Args:
            cdata_json_path: Path to cdata.json metadata file
        """
        self.json_path = cdata_json_path
        self.data = self._load_json()
        all_classes = self.data.get('classes', {})

        # Filter out classes that are already hand-written in base_object/
        # These include fundamental types and base classes
        exclude_classes = {
            # Fundamental types (in fundamental_types.py)
            'CInt', 'CFloat', 'CString', 'CBoolean', 'CList',
            # Base classes (in base_classes.py) - only the hand-written ones
            'CData', 'CContainer', 'CDataFile', 'CDataFileContent',
        }

        self.classes = {
            name: data for name, data in all_classes.items()
            if name not in exclude_classes
        }

        print(f"Loaded {len(all_classes)} classes, generating {len(self.classes)} " +
              f"(excluding {len(exclude_classes)} hand-written base classes)")

        # Initialize type resolver and dependency graph
        print("Initializing type resolver...")
        self.type_resolver = TypeResolver(all_classes)  # Use all classes for type resolution

        print("Building dependency graph...")
        self.dep_graph = ClassDependencyGraph(self.classes, self.type_resolver)

    def _load_json(self) -> Dict[str, Any]:
        """Load and parse cdata.json."""
        print(f"Loading metadata from {self.json_path}...")
        with open(self.json_path, 'r', encoding='utf-8') as f:
            return json.load(f)

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

                # Determine if this is a fundamental type or custom class
                # Map to AttributeType (only fundamental types have enum values)
                # All other types (CFilePath, CUUID, CProjectId, CList, etc.) use AttributeType.CUSTOM
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
                    # Custom class - includes CFilePath, CUUID, CProjectId, CList, etc.
                    lines.append(f'        "{attr_name}": attribute(AttributeType.CUSTOM, custom_class="{attr_type}"),')

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
                    # Normalize type field
                    qdef_copy = qdef.copy()
                    if 'type' in qdef_copy:
                        type_val = qdef_copy['type']
                        if isinstance(type_val, str) and type_val.startswith("<class '"):
                            type_val = type_val[8:-2]  # Extract from "<class 'int'>"
                        # Map to Python types
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

        # Render gui_label
        gui_label = class_data.get('gui_label', '')
        if gui_label:
            lines.append(f'    gui_label="{gui_label}",')

        lines.append(')')
        return lines

    def render_class(self, class_name: str, class_data: dict) -> List[str]:
        """
        Render complete class definition.

        Args:
            class_name: Name of the class
            class_data: Class metadata dictionary

        Returns:
            List of lines for the complete class
        """
        lines = []

        # Add decorator
        lines.extend(self.render_decorator(class_data))

        # Class definition
        parent = class_data.get('immediate_parent', 'CData')
        lines.append(f'class {class_name}({parent}):')

        # Docstring
        docstring = class_data.get('docstring', f'Generated {class_name} class.')
        lines.append(f'    """{docstring}"""')
        lines.append('')

        # Type-annotated attributes
        attributes = class_data.get('CONTENTS', {})
        if attributes:
            for attr_name, attr_info in attributes.items():
                attr_type_str = attr_info.get('class', 'Any')
                attr_type = self.type_resolver.resolve_attribute_type_for_code(attr_type_str)
                lines.append(f'    {attr_name}: Optional[{attr_type}] = None')
            lines.append('')

        # __init__ method
        lines.append('    def __init__(self, parent=None, name=None, **kwargs):')
        lines.append('        """')
        lines.append(f'        Initialize {class_name}.')
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

    def generate_file(self, filename: str, classes_in_file: List[Tuple[str, dict]]) -> str:
        """
        Generate complete file content.

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
        Generate all files.

        Args:
            output_dir: Directory to write generated files
            format_code: Whether to run autopep8 formatting
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get sorted classes by file
        print("\nGenerating dependency-sorted class definitions...")
        sorted_by_file = self.dep_graph.get_sorted_classes_by_file()

        print(f"\nGenerating {len(sorted_by_file)} files...")

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

        print(f"\n✓ Generated {len(sorted_by_file)} files in {output_dir}")

    def _generate_init_file(self, output_dir: Path, sorted_by_file: Dict[str, List[Tuple[str, dict]]]):
        """Generate __init__.py that exports all classes."""
        init_path = output_dir / '__init__.py'

        lines = [
            '"""',
            'Auto-generated CData classes.',
            '',
            'DO NOT EDIT - Generated from CCP4i2 metadata.',
            '"""',
            '',
        ]

        # Import and re-export all classes
        for filename in sorted(sorted_by_file.keys()):
            module_name = filename[:-3]  # Remove .py
            classes_in_file = sorted_by_file[filename]

            if classes_in_file:
                class_names = ', '.join([name for name, _ in classes_in_file])
                lines.append(f'from .{module_name} import {class_names}')

        lines.append('')
        lines.append('__all__ = [')
        for filename in sorted(sorted_by_file.keys()):
            for class_name, _ in sorted_by_file[filename]:
                lines.append(f'    "{class_name}",')
        lines.append(']')
        lines.append('')

        with open(init_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))

        print(f"  Generated {init_path}")

    def verify_imports(self, output_dir: Path) -> bool:
        """
        Verify that all generated files can be imported.

        Args:
            output_dir: Directory containing generated files

        Returns:
            True if all imports succeed, False otherwise
        """
        import sys
        import importlib

        print("\nVerifying imports...")

        # Add parent of output_dir to path
        sys.path.insert(0, str(output_dir.parent.parent))

        success = True
        for py_file in sorted(output_dir.glob('*.py')):
            if py_file.name == '__init__.py':
                continue

            module_name = f'core.generated.{py_file.stem}'
            try:
                importlib.import_module(module_name)
                print(f"  ✓ {py_file.name}")
            except Exception as e:
                print(f"  ✗ {py_file.name}: {e}")
                success = False

        return success


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Generate production CData classes from metadata'
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
        default=Path(__file__).parent.parent.parent / 'core' / 'generated',
        help='Output directory for generated files'
    )
    parser.add_argument(
        '--no-format',
        action='store_true',
        help='Skip autopep8 formatting'
    )
    parser.add_argument(
        '--verify',
        action='store_true',
        help='Verify imports after generation'
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='Print dependency analysis report'
    )

    args = parser.parse_args()

    # Create generator
    generator = ProductionCodeGenerator(args.input)

    # Print dependency report if requested
    if args.report:
        generator.dep_graph.print_dependency_report()

    # Generate all files
    generator.generate_all(args.output, format_code=not args.no_format)

    # Verify imports if requested
    if args.verify:
        if generator.verify_imports(args.output):
            print("\n✓ All imports verified successfully")
        else:
            print("\n✗ Some imports failed")
            return 1

    return 0


if __name__ == '__main__':
    exit(main())
