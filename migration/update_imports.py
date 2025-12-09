#!/usr/bin/env python3
"""
Import Migration Script for CCP4i2

This script updates PySide2 imports to use the baselayer compatibility module,
enabling code to work in both legacy Qt and modern Django environments.

Usage:
    python migration/update_imports.py [--dry-run] [--verbose] [directory]

The script will:
1. Find all Python files with PySide2 imports
2. Replace import statements with baselayer equivalents
3. Create backup files (.bak) for safety

Import patterns handled:
    from PySide2 import QtCore          -> from ccp4i2.baselayer import QtCore
    from PySide2 import QtGui           -> from ccp4i2.baselayer import QtGui
    from PySide2 import QtWidgets       -> from ccp4i2.baselayer import QtWidgets
    from PySide2 import QtSql           -> from ccp4i2.baselayer import QtSql
    from PySide2.QtCore import X, Y     -> from ccp4i2.baselayer.QtCore import X, Y
    from PySide2.QtGui import X, Y      -> from ccp4i2.baselayer.QtGui import X, Y
    from PySide2.QtWidgets import X, Y  -> from ccp4i2.baselayer.QtWidgets import X, Y
    from PySide2.QtSql import X, Y      -> from ccp4i2.baselayer.QtSql import X, Y
"""

import os
import re
import sys
import argparse
import shutil
from pathlib import Path
from typing import List, Tuple, Set


# Import replacement patterns
IMPORT_PATTERNS = [
    # Module imports: from PySide2 import QtCore, QtGui, QtWidgets
    (r'from\s+PySide2\s+import\s+', 'from ccp4i2.baselayer import '),
    # Specific imports: from PySide2.QtCore import Signal, Slot
    (r'from\s+PySide2\.QtCore\s+import\s+', 'from ccp4i2.baselayer.QtCore import '),
    (r'from\s+PySide2\.QtGui\s+import\s+', 'from ccp4i2.baselayer.QtGui import '),
    (r'from\s+PySide2\.QtWidgets\s+import\s+', 'from ccp4i2.baselayer.QtWidgets import '),
    (r'from\s+PySide2\.QtSql\s+import\s+', 'from ccp4i2.baselayer.QtSql import '),
    # Direct module reference (less common)
    (r'import\s+PySide2\.QtCore\b', 'from ccp4i2.baselayer import QtCore'),
    (r'import\s+PySide2\.QtGui\b', 'from ccp4i2.baselayer import QtGui'),
    (r'import\s+PySide2\.QtWidgets\b', 'from ccp4i2.baselayer import QtWidgets'),
    (r'import\s+PySide2\.QtSql\b', 'from ccp4i2.baselayer import QtSql'),
    (r'import\s+PySide2\b', 'import ccp4i2.baselayer'),
]

# Directories to process (locked directories from the plan)
DEFAULT_DIRECTORIES = [
    'wrappers',
    'wrappers2',
    'pipelines',
    'pimple',
    'utils',
]

# Directories to skip (will be handled separately or already migrated)
SKIP_DIRECTORIES = {
    '__pycache__',
    '.git',
    'node_modules',
    'venv',
    '.venv',
    'build',
    'dist',
}

# Files to skip (baselayer itself, already migrated files)
SKIP_FILES = {
    'baselayer/__init__.py',
    'baselayer/QtCore.py',
    'baselayer/QtGui.py',
    'baselayer/QtWidgets.py',
    'baselayer/QtSql.py',
}


def find_python_files(directory: str) -> List[Path]:
    """Find all Python files in directory recursively."""
    python_files = []
    base_path = Path(directory)

    for path in base_path.rglob('*.py'):
        # Skip directories in SKIP_DIRECTORIES
        if any(skip in path.parts for skip in SKIP_DIRECTORIES):
            continue
        # Skip specific files
        relative_path = str(path.relative_to(base_path.parent))
        if relative_path in SKIP_FILES:
            continue
        python_files.append(path)

    return sorted(python_files)


def has_pyside2_imports(content: str) -> bool:
    """Check if content has any PySide2 imports."""
    return 'PySide2' in content


def update_imports(content: str) -> Tuple[str, List[str]]:
    """
    Update PySide2 imports to baselayer imports.

    Returns:
        Tuple of (updated_content, list_of_changes)
    """
    changes = []
    updated = content

    for pattern, replacement in IMPORT_PATTERNS:
        matches = re.findall(pattern, content, re.MULTILINE)
        if matches:
            updated = re.sub(pattern, replacement, updated, flags=re.MULTILINE)
            for match in set(matches):
                changes.append(f"  {match.strip()} -> {replacement.strip()}")

    return updated, changes


def process_file(filepath: Path, dry_run: bool = False, verbose: bool = False) -> Tuple[bool, List[str]]:
    """
    Process a single file, updating imports if needed.

    Returns:
        Tuple of (was_modified, list_of_changes)
    """
    try:
        content = filepath.read_text(encoding='utf-8')
    except UnicodeDecodeError:
        try:
            content = filepath.read_text(encoding='latin-1')
        except Exception as e:
            return False, [f"Error reading file: {e}"]

    if not has_pyside2_imports(content):
        return False, []

    updated, changes = update_imports(content)

    if updated == content:
        return False, []

    if not dry_run:
        # Create backup
        backup_path = filepath.with_suffix('.py.bak')
        if not backup_path.exists():
            shutil.copy2(filepath, backup_path)
        # Write updated content
        filepath.write_text(updated, encoding='utf-8')

    return True, changes


def main():
    parser = argparse.ArgumentParser(
        description='Update PySide2 imports to baselayer imports'
    )
    parser.add_argument(
        'directories',
        nargs='*',
        default=DEFAULT_DIRECTORIES,
        help='Directories to process (default: wrappers, wrappers2, pipelines, pimple, utils)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be changed without making changes'
    )
    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Show detailed changes'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Process all directories (including core, qtgui, qtcore, etc.)'
    )

    args = parser.parse_args()

    # Determine base directory (ccp4i2 root)
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent

    if args.all:
        directories = [base_dir]
    else:
        directories = [base_dir / d for d in args.directories]

    total_files = 0
    modified_files = 0
    all_changes = []

    print(f"{'[DRY RUN] ' if args.dry_run else ''}Updating PySide2 imports to baselayer...")
    print(f"Processing directories: {', '.join(str(d.name) for d in directories)}")
    print()

    for directory in directories:
        if not directory.exists():
            print(f"Warning: Directory {directory} does not exist, skipping")
            continue

        python_files = find_python_files(directory)

        for filepath in python_files:
            total_files += 1
            was_modified, changes = process_file(filepath, args.dry_run, args.verbose)

            if was_modified:
                modified_files += 1
                relative_path = filepath.relative_to(base_dir)
                print(f"{'Would modify' if args.dry_run else 'Modified'}: {relative_path}")

                if args.verbose and changes:
                    for change in changes:
                        print(f"  {change}")

                all_changes.append((relative_path, changes))

    print()
    print(f"Summary:")
    print(f"  Total files scanned: {total_files}")
    print(f"  Files {'that would be modified' if args.dry_run else 'modified'}: {modified_files}")

    if not args.dry_run and modified_files > 0:
        print()
        print("Backup files created with .bak extension")
        print("To remove backups: find . -name '*.py.bak' -delete")


if __name__ == '__main__':
    main()
