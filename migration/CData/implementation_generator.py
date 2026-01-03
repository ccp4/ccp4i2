"""
Implementation Class Generator

Generates skeleton implementation classes in core/ that extend stub classes.
These are meant to be edited to add methods and business logic.
"""

import json
import argparse
from pathlib import Path
from typing import Dict, Any, List, Tuple


class ImplementationGenerator:
    """Generates implementation classes that extend stubs."""

    def __init__(self, cdata_json_path: Path):
        """
        Initialize generator.

        Args:
            cdata_json_path: Path to cdata.json metadata file
        """
        self.json_path = cdata_json_path
        self.data = self._load_json()
        all_classes = self.data.get('classes', {})

        # Base classes that don't need implementation (hand-written)
        self.base_classes = {
            'CData', 'CContainer', 'CDataFile', 'CDataFileContent',
            'CInt', 'CFloat', 'CString', 'CBoolean', 'CList',
        }

        # Filter out base classes
        self.classes = {
            name: data for name, data in all_classes.items()
            if name not in self.base_classes
        }

        print(f"Loaded {len(all_classes)} classes, generating {len(self.classes)} implementations " +
              f"(excluding {len(self.base_classes)} hand-written base classes)")

        # Build class-to-file mapping
        self.class_to_file = self._build_class_to_file_map()

    def _load_json(self) -> Dict[str, Any]:
        """Load and parse cdata.json."""
        print(f"Loading metadata from {self.json_path}...")
        with open(self.json_path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def _build_class_to_file_map(self) -> Dict[str, str]:
        """Build mapping of class name to output file name."""
        mapping = {}
        for class_name, class_info in self.classes.items():
            file_path = class_info.get('file_path', 'unknown.py')
            basename = Path(file_path).name
            mapping[class_name] = basename
        return mapping

    def render_class(self, class_name: str, class_data: dict) -> List[str]:
        """
        Render implementation class definition.

        Args:
            class_name: Name of the class
            class_data: Class metadata dictionary

        Returns:
            List of lines for the implementation class
        """
        lines = []
        stub_name = f"{class_name}Stub"

        # Class definition extending stub
        lines.append(f'class {class_name}({stub_name}):')

        # Docstring
        docstring = class_data.get('docstring', f'Implementation of {class_name}.')
        lines.append(f'    """')
        lines.append(f'    {docstring}')
        lines.append(f'    ')
        lines.append(f'    Extends {stub_name} with implementation-specific methods.')
        lines.append(f'    Add file I/O, validation, and business logic here.')
        lines.append(f'    """')
        lines.append('')

        # Placeholder pass or minimal __init__
        lines.append('    # Add your methods here')
        lines.append('    pass')

        lines.append('')
        lines.append('')

        return lines

    def generate_file(self, filename: str, classes_in_file: List[Tuple[str, dict]]) -> str:
        """
        Generate complete file content for implementation.

        Args:
            filename: Output filename
            classes_in_file: List of (class_name, class_data) tuples

        Returns:
            Complete file content as string
        """
        lines = []

        # File header
        lines.append('"""')
        lines.append(f'Implementation classes for {filename}')
        lines.append('')
        lines.append('Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.')
        lines.append('This file is safe to edit - add your implementation code here.')
        lines.append('"""')
        lines.append('')
        lines.append('from __future__ import annotations')
        lines.append('from typing import Optional, Any')
        lines.append('')

        # Import all stub classes from the corresponding stub file
        stub_module = f"core.cdata_stubs.{filename[:-3]}"
        stub_names = [f"{name}Stub" for name, _ in classes_in_file]
        stub_list = ', '.join(stub_names)
        lines.append(f'from {stub_module} import {stub_list}')
        lines.append('')
        lines.append('')

        # Generate each implementation class
        for class_name, class_data in classes_in_file:
            class_lines = self.render_class(class_name, class_data)
            lines.extend(class_lines)

        return '\n'.join(lines)

    def generate_all(self, output_dir: Path, overwrite: bool = False):
        """
        Generate all implementation files.

        Args:
            output_dir: Directory to write generated files (core/)
            overwrite: Whether to overwrite existing files
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Group classes by file
        classes_by_file: Dict[str, List[Tuple[str, dict]]] = {}
        for class_name, class_data in self.classes.items():
            filename = self.class_to_file[class_name]
            if filename not in classes_by_file:
                classes_by_file[filename] = []
            classes_by_file[filename].append((class_name, class_data))

        print(f"\nGenerating {len(classes_by_file)} implementation files...")

        generated_count = 0
        skipped_count = 0

        for filename, classes_in_file in sorted(classes_by_file.items()):
            output_path = output_dir / filename

            # Skip if file exists and overwrite is False
            if output_path.exists() and not overwrite:
                print(f"  Skipping {filename} (already exists)")
                skipped_count += 1
                continue

            print(f"  Writing {filename} ({len(classes_in_file)} classes)...")

            # Generate file content
            content = self.generate_file(filename, classes_in_file)

            # Write to file
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(content)

            generated_count += 1

        # Generate __init__.py
        self._generate_init_file(output_dir, classes_by_file, overwrite)

        print(f"\n✓ Generated {generated_count} files, skipped {skipped_count} existing files in {output_dir}")

    def _generate_init_file(self, output_dir: Path, classes_by_file: Dict[str, List[Tuple[str, dict]]], overwrite: bool):
        """Generate __init__.py that exports all implementation classes."""
        init_path = output_dir / '__init__.py'

        # Don't overwrite __init__.py if it exists (might have custom exports)
        if init_path.exists() and not overwrite:
            print(f"  Skipping __init__.py (already exists)")
            return

        lines = [
            '"""',
            'CData implementation classes.',
            '',
            'These classes extend the stubs from ccp4i2.core.cdata_stubs with methods and business logic.',
            '"""',
            '',
        ]

        # Import and re-export all classes
        for filename in sorted(classes_by_file.keys()):
            module_name = filename[:-3]  # Remove .py
            classes_in_file = classes_by_file[filename]

            if classes_in_file:
                class_names = ', '.join([name for name, _ in classes_in_file])
                lines.append(f'from .{module_name} import {class_names}')

        lines.append('')
        lines.append('__all__ = [')
        for filename in sorted(classes_by_file.keys()):
            for class_name, _ in classes_by_file[filename]:
                lines.append(f'    "{class_name}",')
        lines.append(']')
        lines.append('')

        with open(init_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))

        print(f"  Generated {init_path}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Generate implementation classes that extend stubs'
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
        default=Path(__file__).parent.parent.parent / 'core',
        help='Output directory for implementation files'
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Overwrite existing implementation files (WARNING: will lose custom code)'
    )

    args = parser.parse_args()

    if args.overwrite:
        print("WARNING: --overwrite is enabled. Existing files will be overwritten!")
        response = input("Continue? [y/N]: ")
        if response.lower() != 'y':
            print("Aborted.")
            return 1

    # Create generator
    generator = ImplementationGenerator(args.input)

    # Generate all implementation files
    generator.generate_all(args.output, overwrite=args.overwrite)

    return 0


if __name__ == '__main__':
    exit(main())
