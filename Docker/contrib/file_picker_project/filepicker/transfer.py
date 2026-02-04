"""
Transfer manifest generators for rsync, globus, and azcopy.

Generates file lists and command scripts for external transfer tools.
"""

import os
import stat
from typing import Dict, List, Optional

from .resolver import ResolvedFile


def generate_file_list(resolved_files, output_path, absolute=True):
    # type: (Dict[int, List[ResolvedFile]], str, bool) -> str
    """
    Generate a plain text file list (one file per line).

    This can be used with rsync --files-from or other tools.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output file list
        absolute: If True, write absolute paths; if False, write source paths

    Returns:
        Path to the generated file list
    """
    seen = set()
    lines = []

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]
        for f in files:
            if f.source_path not in seen:
                seen.add(f.source_path)
                lines.append(f.source_path)

    with open(output_path, 'w') as fh:
        fh.write('\n'.join(lines))
        fh.write('\n')

    return output_path


def generate_rsync_script(resolved_files, output_path, destination,
                          base_dir=None, rsync_options=None):
    # type: (Dict[int, List[ResolvedFile]], str, str, Optional[str], Optional[List[str]]) -> str
    """
    Generate an rsync script for transferring files.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output script
        destination: rsync destination (e.g., user@host:/path or /local/path)
        base_dir: Base directory to use for relative paths
        rsync_options: Additional rsync options

    Returns:
        Path to the generated script
    """
    if rsync_options is None:
        rsync_options = ['-avP', '--progress']

    # Generate file list
    file_list_path = output_path + '.files'
    generate_file_list(resolved_files, file_list_path, absolute=True)

    # Build the script
    lines = [
        '#!/bin/bash',
        '# Generated rsync transfer script',
        '# Files: {}'.format(sum(len(f) for f in resolved_files.values())),
        '',
        'set -e',
        '',
        'FILE_LIST="{}"'.format(file_list_path),
        'DESTINATION="{}"'.format(destination),
        '',
        '# Ensure destination directory exists (for local transfers)',
        'if [[ "$DESTINATION" != *:* ]]; then',
        '    mkdir -p "$DESTINATION"',
        'fi',
        '',
        '# Run rsync',
        'rsync {} \\'.format(' '.join(rsync_options)),
        '    --files-from="$FILE_LIST" \\',
        '    / \\',
        '    "$DESTINATION"',
        '',
        'echo "Transfer complete!"',
    ]

    with open(output_path, 'w') as fh:
        fh.write('\n'.join(lines))
        fh.write('\n')

    # Make executable
    os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)

    return output_path


def generate_rsync_with_structure(resolved_files, output_path, destination,
                                   staging_dir=None):
    # type: (Dict[int, List[ResolvedFile]], str, str, Optional[str]) -> str
    """
    Generate an rsync script that preserves the output directory structure.

    This first copies files to a staging directory with the desired structure,
    then rsyncs the staging directory.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output script
        destination: rsync destination
        staging_dir: Directory for staging files (default: creates temp dir)

    Returns:
        Path to the generated script
    """
    if staging_dir is None:
        staging_dir = '/tmp/filepicker_staging_$$'

    # Build copy commands
    copy_commands = []
    seen = set()

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]
        for f in files:
            if f.output_path not in seen:
                seen.add(f.output_path)
                dest_path = os.path.join('$STAGING_DIR', f.output_path)
                dest_dir = os.path.dirname(dest_path)
                copy_commands.append('mkdir -p "{}"'.format(dest_dir))
                copy_commands.append('cp -L "{}" "{}"'.format(
                    f.source_path, dest_path))

    lines = [
        '#!/bin/bash',
        '# Generated rsync transfer script with directory structure',
        '# Files: {}'.format(len(seen)),
        '',
        'set -e',
        '',
        'STAGING_DIR="{}"'.format(staging_dir),
        'DESTINATION="{}"'.format(destination),
        '',
        '# Create staging directory',
        'mkdir -p "$STAGING_DIR"',
        '',
        '# Copy files to staging (dereferencing symlinks)',
        'echo "Staging files..."',
    ] + copy_commands + [
        '',
        '# Transfer staged files',
        'echo "Transferring files..."',
        'rsync -avP --progress "$STAGING_DIR/" "$DESTINATION"',
        '',
        '# Cleanup staging',
        'rm -rf "$STAGING_DIR"',
        '',
        'echo "Transfer complete!"',
    ]

    with open(output_path, 'w') as fh:
        fh.write('\n'.join(lines))
        fh.write('\n')

    os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)

    return output_path


def generate_globus_manifest(resolved_files, output_path,
                             source_endpoint=None, dest_endpoint=None,
                             source_base=None, dest_base=None):
    # type: (Dict[int, List[ResolvedFile]], str, Optional[str], Optional[str], Optional[str], Optional[str]) -> str
    """
    Generate a Globus transfer manifest file.

    The manifest is a CSV file that can be used with the Globus CLI:
        globus transfer --batch manifest.csv SOURCE_EP DEST_EP

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output manifest
        source_endpoint: Globus source endpoint UUID (optional, for comments)
        dest_endpoint: Globus destination endpoint UUID (optional)
        source_base: Base path on source endpoint
        dest_base: Base path on destination endpoint

    Returns:
        Path to the generated manifest
    """
    lines = [
        '# Globus transfer manifest',
        '# Generated by filepicker',
        '# Files: {}'.format(sum(len(f) for f in resolved_files.values())),
    ]

    if source_endpoint:
        lines.append('# Source endpoint: {}'.format(source_endpoint))
    if dest_endpoint:
        lines.append('# Dest endpoint: {}'.format(dest_endpoint))

    lines.append('#')
    lines.append('# Format: source_path,dest_path')
    lines.append('')

    seen = set()
    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]
        for f in files:
            if f.output_path not in seen:
                seen.add(f.output_path)

                # Build source and dest paths
                if source_base:
                    src = f.source_path
                    if src.startswith(source_base):
                        src = src[len(source_base):].lstrip('/')
                else:
                    src = f.source_path

                if dest_base:
                    dst = os.path.join(dest_base, f.output_path)
                else:
                    dst = f.output_path

                lines.append('{},{}'.format(src, dst))

    with open(output_path, 'w') as fh:
        fh.write('\n'.join(lines))
        fh.write('\n')

    return output_path


def generate_azcopy_script(resolved_files, output_path, container_url,
                           dest_prefix=None):
    # type: (Dict[int, List[ResolvedFile]], str, str, Optional[str]) -> str
    """
    Generate an azcopy script for transferring to Azure Blob Storage.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output script
        container_url: Azure blob container URL (with SAS token if needed)
        dest_prefix: Prefix path in the container

    Returns:
        Path to the generated script
    """
    lines = [
        '#!/bin/bash',
        '# Generated azcopy transfer script',
        '# Files: {}'.format(sum(len(f) for f in resolved_files.values())),
        '',
        'set -e',
        '',
        'CONTAINER_URL="{}"'.format(container_url),
    ]

    if dest_prefix:
        lines.append('DEST_PREFIX="{}"'.format(dest_prefix))
    else:
        lines.append('DEST_PREFIX=""')

    lines.extend([
        '',
        '# Transfer files',
    ])

    seen = set()
    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]
        for f in files:
            if f.output_path not in seen:
                seen.add(f.output_path)

                if dest_prefix:
                    dest_path = '{}/{}'.format(dest_prefix.rstrip('/'),
                                               f.output_path)
                else:
                    dest_path = f.output_path

                lines.append('azcopy copy "{}" "$CONTAINER_URL/{}"'.format(
                    f.source_path, dest_path))

    lines.extend([
        '',
        'echo "Transfer complete!"',
    ])

    with open(output_path, 'w') as fh:
        fh.write('\n'.join(lines))
        fh.write('\n')

    os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IXUSR | stat.S_IXGRP)

    return output_path


def generate_manifest_json(resolved_files, output_path, include_metadata=True):
    # type: (Dict[int, List[ResolvedFile]], str, bool) -> str
    """
    Generate a JSON manifest of all files.

    Useful for programmatic consumption or as input to other tools.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output JSON file
        include_metadata: Include file sizes and modification times

    Returns:
        Path to the generated manifest
    """
    import json

    manifest = {
        'version': '1.0',
        'crystals': {},
    }

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]
        crystal_key = 'x{:04d}'.format(index)
        manifest['crystals'][crystal_key] = []

        for f in files:
            entry = {
                'source': f.source_path,
                'output': f.output_path,
                'pattern': f.original_pattern,
                'is_sibling': f.is_sibling,
            }

            if include_metadata:
                try:
                    stat_info = os.stat(f.source_path)
                    entry['size'] = stat_info.st_size
                    entry['mtime'] = stat_info.st_mtime
                except OSError:
                    pass

            manifest['crystals'][crystal_key].append(entry)

    with open(output_path, 'w') as fh:
        json.dump(manifest, fh, indent=2)

    return output_path
