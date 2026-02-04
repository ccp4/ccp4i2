"""
Command-line interface for filepicker.

Supports both local and remote (SFTP) sources.

Usage examples:
    # Remote: Show available crystals on Diamond
    filepicker fedid@nx.diamond.ac.uk:/dls/i04-1/data/2024/sw12345-6/processing/analysis/model_building info

    # Remote: Download selected crystals
    filepicker fedid@nx.diamond.ac.uk:/path/to/model_building pull --crystals x0001-x0010 --dest ./local_data

    # Local: Create archive (if you have local access)
    filepicker /local/path/to/model_building archive --crystals x0001-x0010 -o data.tar.gz
"""

import argparse
import os
import re
import sys
from typing import List, Optional, Tuple

from . import __version__
from .config import get_template, parse_custom_template, TEMPLATES
from .resolver import format_size


# Pattern to match remote source: user@host:/path
REMOTE_PATTERN = re.compile(r'^([^@]+)@([^:]+):(.+)$')


def parse_source(source):
    # type: (str) -> Tuple[bool, str, Optional[str], Optional[str]]
    """
    Parse source string to determine if it's local or remote.

    Returns: (is_remote, path, username, host)
    """
    match = REMOTE_PATTERN.match(source)
    if match:
        username, host, path = match.groups()
        return (True, path, username, host)
    else:
        return (False, source, None, None)


def print_progress(progress):
    """Print progress bar to stderr."""
    bar_width = 40
    filled = int(bar_width * progress.percent_bytes / 100)
    bar = '=' * filled + '-' * (bar_width - filled)

    sys.stderr.write('\r[{}] {:.1f}% ({}/{} files)'.format(
        bar, progress.percent_bytes,
        progress.files_done, progress.total_files))
    sys.stderr.flush()

    if progress.files_done == progress.total_files:
        sys.stderr.write('\n')


def print_status(message):
    """Print a status message."""
    sys.stderr.write('\r' + message + '          \n')
    sys.stderr.flush()


# ============================================================================
# Local commands (original functionality)
# ============================================================================

def cmd_info_local(args, path):
    """Show information about available crystals (local)."""
    from .selector import CrystalSelector

    selector = CrystalSelector(path, crystal_prefix=args.prefix)
    print(selector.summary())

    if args.verbose:
        print("\nAvailable templates:")
        for name in TEMPLATES:
            print("  - {}".format(name))


def cmd_list_local(args, path):
    """List selected crystals (local)."""
    from .selector import CrystalSelector

    selector = CrystalSelector(path, crystal_prefix=args.prefix)

    if args.crystals:
        indices = selector.select(args.crystals)
    else:
        indices = selector.get_all_indices()

    print("Selected {} crystals:".format(len(indices)))
    for idx in indices:
        name = selector.get_crystal_dir_name(idx)
        print("  x{:04d}  {}".format(idx, name))


def cmd_resolve_local(args, path):
    """Resolve and display files that would be selected (local)."""
    from .selector import CrystalSelector
    from .resolver import resolve_crystals, calculate_total_size

    selector = CrystalSelector(path, crystal_prefix=args.prefix)

    if args.crystals:
        indices = selector.select(args.crystals)
    else:
        indices = selector.get_all_indices()

    template = _get_template(args)
    resolved = resolve_crystals(selector, indices, template)

    _print_resolved_summary(resolved, args.verbose, local=True)


def cmd_archive_local(args, path):
    """Create an archive of selected files (local)."""
    from .selector import CrystalSelector
    from .resolver import resolve_crystals
    from .archiver import (
        create_tar_archive, create_zip_archive, estimate_archive_size
    )

    selector = CrystalSelector(path, crystal_prefix=args.prefix)

    if args.crystals:
        indices = selector.select(args.crystals)
    else:
        indices = selector.get_all_indices()

    template = _get_template(args)

    print("Resolving files for {} crystals...".format(len(indices)))
    resolved = resolve_crystals(selector, indices, template)

    estimates = estimate_archive_size(resolved)
    print("Files to archive: {}".format(estimates['file_count']))
    print("Raw size: {}".format(format_size(estimates['raw'])))
    print("Estimated compressed: ~{}".format(format_size(estimates['gzip_estimate'])))

    output = args.output
    progress_cb = print_progress if not args.quiet else None

    if output.endswith('.zip'):
        print("\nCreating zip archive: {}".format(output))
        result = create_zip_archive(resolved, output,
                                    progress_callback=progress_cb,
                                    archive_name=args.archive_name)
    elif output.endswith('.tar.bz2'):
        print("\nCreating tar.bz2 archive: {}".format(output))
        result = create_tar_archive(resolved, output, compression='bz2',
                                    progress_callback=progress_cb,
                                    archive_name=args.archive_name)
    else:
        if not output.endswith('.tar.gz'):
            output = output + '.tar.gz'
        print("\nCreating tar.gz archive: {}".format(output))
        result = create_tar_archive(resolved, output, compression='gz',
                                    progress_callback=progress_cb,
                                    archive_name=args.archive_name)

    final_size = os.path.getsize(result)
    print("Archive created: {} ({})".format(result, format_size(final_size)))


# ============================================================================
# Remote commands (SFTP-based)
# ============================================================================

def cmd_info_remote(args, path, username, host):
    """Show information about available crystals (remote)."""
    from .remote import SFTPScanner, RemoteCrystalScanner

    print("Connecting to {}@{}...".format(username, host))

    with SFTPScanner(host, username) as sftp:
        scanner = RemoteCrystalScanner(sftp, path, crystal_prefix=args.prefix)
        print(scanner.summary())

        if args.verbose:
            print("\nAvailable templates:")
            for name in TEMPLATES:
                print("  - {}".format(name))


def cmd_list_remote(args, path, username, host):
    """List selected crystals (remote)."""
    from .remote import SFTPScanner, RemoteCrystalScanner
    from .selector import CrystalSelector

    print("Connecting to {}@{}...".format(username, host))

    with SFTPScanner(host, username) as sftp:
        scanner = RemoteCrystalScanner(sftp, path, crystal_prefix=args.prefix)

        if args.crystals:
            # Use local selector logic with remote data
            indices = _select_crystals_remote(scanner, args.crystals)
        else:
            indices = scanner.get_all_indices()

        print("Selected {} crystals:".format(len(indices)))
        for idx in indices:
            name = scanner.get_crystal_dir_name(idx)
            print("  x{:04d}  {}".format(idx, name))


def cmd_resolve_remote(args, path, username, host):
    """Resolve and display files that would be selected (remote)."""
    from .remote import SFTPScanner, RemoteCrystalScanner
    from .remote_resolver import resolve_crystals_remote, calculate_total_size_remote

    print("Connecting to {}@{}...".format(username, host))

    with SFTPScanner(host, username) as sftp:
        scanner = RemoteCrystalScanner(sftp, path, crystal_prefix=args.prefix)

        if args.crystals:
            indices = _select_crystals_remote(scanner, args.crystals)
        else:
            indices = scanner.get_all_indices()

        template = _get_template(args)

        print("Resolving files for {} crystals...".format(len(indices)))
        resolved = resolve_crystals_remote(
            sftp, scanner, indices, template,
            progress_callback=print_status if args.verbose else None
        )

        _print_resolved_summary_remote(resolved, args.verbose)


def cmd_pull_remote(args, path, username, host):
    """Download selected files from remote server."""
    from .remote import SFTPScanner, RemoteCrystalScanner
    from .remote_resolver import resolve_crystals_remote, calculate_total_size_remote
    from .sftp_transfer import (
        download_files, print_transfer_progress, write_transfer_manifest,
        parse_manifest, download_from_manifest
    )

    # Check if using pre-generated manifest
    from_manifest = getattr(args, 'from_manifest', None)

    print("Connecting to {}@{}...".format(username, host))

    with SFTPScanner(host, username) as sftp:
        if from_manifest:
            # Use pre-generated manifest - skip SFTP scanning
            print("Reading manifest: {}".format(from_manifest))
            entries = parse_manifest(from_manifest)

            total_files = len(entries)
            total_size = sum(e.size for e in entries)

            print("Files to download: {} ({})".format(total_files, format_size(total_size)))

            # If dry-run, stop here
            if args.dry_run:
                print("\nDry run - no files downloaded.")
                if args.manifest:
                    print("Note: --manifest ignored when using --from-manifest")
                return

            # Confirm before large downloads
            if total_size > 1024 * 1024 * 1024 and not args.yes:  # 1GB
                response = input("\nDownload {} of data? [y/N] ".format(
                    format_size(total_size)))
                if response.lower() != 'y':
                    print("Cancelled.")
                    return

            # Download files
            dest = args.dest
            print("\nDownloading to: {}".format(dest))

            progress_cb = print_transfer_progress if not args.quiet else None
            progress = download_from_manifest(
                sftp, entries, dest,
                progress_callback=progress_cb,
                skip_existing=args.skip_existing
            )

        else:
            # Standard mode: scan and resolve via SFTP
            scanner = RemoteCrystalScanner(sftp, path, crystal_prefix=args.prefix)

            if args.crystals:
                indices = _select_crystals_remote(scanner, args.crystals)
            else:
                indices = scanner.get_all_indices()

            template = _get_template(args)

            print("Resolving files for {} crystals...".format(len(indices)))
            resolved = resolve_crystals_remote(sftp, scanner, indices, template)

            total_size = calculate_total_size_remote(resolved)
            total_files = sum(len(f) for f in resolved.values())

            print("Files to download: {} ({})".format(total_files, format_size(total_size)))

            # Write manifest if requested
            if args.manifest:
                manifest_path = args.manifest
                write_transfer_manifest(resolved, manifest_path)
                print("Manifest written: {}".format(manifest_path))

            # If dry-run, stop here
            if args.dry_run:
                print("\nDry run - no files downloaded.")
                return

            # Confirm before large downloads
            if total_size > 1024 * 1024 * 1024 and not args.yes:  # 1GB
                response = input("\nDownload {} of data? [y/N] ".format(
                    format_size(total_size)))
                if response.lower() != 'y':
                    print("Cancelled.")
                    return

            # Download files
            dest = args.dest
            print("\nDownloading to: {}".format(dest))

            progress_cb = print_transfer_progress if not args.quiet else None
            progress = download_files(
                sftp, resolved, dest,
                progress_callback=progress_cb,
                skip_existing=args.skip_existing
            )

        print("\nDownload complete!")
        print("  Files: {}/{}".format(progress.files_done, progress.total_files))
        print("  Size: {}".format(format_size(progress.bytes_done)))
        print("  Time: {:.1f}s".format(progress.elapsed_seconds))

        if progress.errors:
            print("\nErrors ({}):", len(progress.errors))
            for err in progress.errors[:10]:
                print("  - {}".format(err))
            if len(progress.errors) > 10:
                print("  ... and {} more".format(len(progress.errors) - 10))


def cmd_manifest_remote(args, path, username, host):
    """Generate a transfer manifest without downloading."""
    from .remote import SFTPScanner, RemoteCrystalScanner
    from .remote_resolver import resolve_crystals_remote, calculate_total_size_remote
    from .sftp_transfer import write_transfer_manifest, generate_rsync_pull_command

    print("Connecting to {}@{}...".format(username, host))

    with SFTPScanner(host, username) as sftp:
        scanner = RemoteCrystalScanner(sftp, path, crystal_prefix=args.prefix)

        if args.crystals:
            indices = _select_crystals_remote(scanner, args.crystals)
        else:
            indices = scanner.get_all_indices()

        template = _get_template(args)

        print("Resolving files for {} crystals...".format(len(indices)))
        resolved = resolve_crystals_remote(sftp, scanner, indices, template)

        total_size = calculate_total_size_remote(resolved)
        total_files = sum(len(f) for f in resolved.values())

        print("Files: {} ({})".format(total_files, format_size(total_size)))

        output_base = args.output or 'transfer'

        # Write manifest
        manifest_path = output_base + '_manifest.txt'
        write_transfer_manifest(resolved, manifest_path)
        print("Manifest: {}".format(manifest_path))

        # Write rsync command
        if args.rsync:
            file_list = output_base + '_files.txt'
            rsync_cmd = generate_rsync_pull_command(
                resolved, host, username,
                args.rsync, file_list_path=file_list
            )
            print("File list: {}".format(file_list))
            print("\nRsync command:")
            print("  {}".format(rsync_cmd))


# ============================================================================
# Helper functions
# ============================================================================

def _get_template(args):
    """Get file template from args."""
    if hasattr(args, 'template') and args.template:
        return get_template(args.template)
    elif hasattr(args, 'files') and args.files:
        return parse_custom_template(args.files)
    else:
        return get_template('default')


def _select_crystals_remote(scanner, specifiers):
    # type: (any, List[str]) -> List[int]
    """
    Select crystals using specifiers against a remote scanner.

    Simplified version of local selector logic.
    """
    import re

    crystal_spec = re.compile(r'^x?(\d+)$')
    range_spec = re.compile(r'^x?(\d+)[-:]x?(\d+)$')
    compound_spec = re.compile(r'^[A-Z]{2,4}-\d+$')

    all_indices = set(scanner.get_all_indices())
    included = set()
    excluded = set()

    for spec in specifiers:
        spec = spec.strip()
        if not spec:
            continue

        is_exclude = spec.startswith('!')
        if is_exclude:
            spec = spec[1:]

        matches = set()

        if spec == '*' or spec == 'all':
            matches = all_indices.copy()

        elif range_spec.match(spec):
            m = range_spec.match(spec)
            start, end = int(m.group(1)), int(m.group(2))
            if start > end:
                start, end = end, start
            matches = set(i for i in range(start, end + 1) if i in all_indices)

        elif crystal_spec.match(spec):
            m = crystal_spec.match(spec)
            idx = int(m.group(1))
            if idx in all_indices:
                matches = {idx}

        elif compound_spec.match(spec):
            # Would need compound scanning - skip for now
            print("Warning: Compound selection not yet implemented for remote")

        if is_exclude:
            excluded.update(matches)
        else:
            included.update(matches)

    if not included:
        included = all_indices

    return sorted(included - excluded)


def _print_resolved_summary(resolved, verbose, local=True):
    """Print summary of resolved files (local)."""
    from .resolver import calculate_total_size

    total_files = 0

    for idx in sorted(resolved.keys()):
        files = resolved[idx]
        if verbose:
            print("\nx{:04d} ({} files):".format(idx, len(files)))
            for f in files:
                try:
                    size = os.path.getsize(f.source_path)
                    print("  {} -> {} ({})".format(
                        f.source_path, f.output_path, format_size(size)))
                except OSError:
                    print("  {} -> {}".format(f.source_path, f.output_path))
        total_files += len(files)

    total_size = calculate_total_size(resolved)

    print("\nSummary:")
    print("  Crystals: {}".format(len(resolved)))
    print("  Files: {}".format(total_files))
    print("  Total size: {}".format(format_size(total_size)))


def _print_resolved_summary_remote(resolved, verbose):
    """Print summary of resolved files (remote)."""
    from .remote_resolver import calculate_total_size_remote

    total_files = 0
    crystals_with_files = 0
    crystals_empty = 0

    for idx in sorted(resolved.keys()):
        files = resolved[idx]
        if files:
            crystals_with_files += 1
            if verbose:
                print("\nx{:04d} ({} files):".format(idx, len(files)))
                for f in files:
                    print("  {} -> {} ({})".format(
                        f.remote_path, f.local_path, format_size(f.size)))
        else:
            crystals_empty += 1
            if verbose:
                print("\nx{:04d}: (no files)".format(idx))
        total_files += len(files)

    total_size = calculate_total_size_remote(resolved)

    print("\nSummary:")
    print("  Crystals selected: {}".format(len(resolved)))
    print("  Crystals with files: {}".format(crystals_with_files))
    if crystals_empty:
        print("  Crystals without matching files: {}".format(crystals_empty))
    print("  Total files: {}".format(total_files))
    print("  Total size: {}".format(format_size(total_size)))


# ============================================================================
# Argument parser
# ============================================================================

def create_parser():
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        prog='filepicker',
        description='Select and transfer crystallography data files.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (Remote - Diamond Light Source):
  # Show available crystals
  %(prog)s fedid@nx.diamond.ac.uk:/dls/.../model_building info

  # List crystals in a range
  %(prog)s fedid@nx.diamond.ac.uk:/dls/.../model_building list -c x0001-x0010

  # Preview files that would be selected
  %(prog)s fedid@nx.diamond.ac.uk:/dls/.../model_building resolve -c x0001 -v

  # Download selected crystals
  %(prog)s fedid@nx.diamond.ac.uk:/dls/.../model_building pull -c x0001-x0050 -d ./local_data

  # Generate manifest without downloading
  %(prog)s fedid@nx.diamond.ac.uk:/dls/.../model_building manifest -c '*' -o transfer

Examples (Local):
  # Create archive from local path
  %(prog)s /local/model_building archive -c x0001-x0010 -o data.tar.gz

Crystal specifiers:
  x0001           Single crystal by index
  x0001-x0010     Range of crystals (inclusive)
  *               All crystals
  !x0005          Exclude crystal

Source format:
  /local/path                     Local directory
  user@host:/remote/path          Remote via SFTP (e.g., fedid@nx.diamond.ac.uk:/dls/...)
"""
    )

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    parser.add_argument('source', metavar='SOURCE',
                        help='Path to model_building directory (local or user@host:/path)')

    parser.add_argument('--prefix', '-p', metavar='PREFIX',
                        help='Crystal directory prefix (auto-detected if not specified)')

    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Command to run')

    # Info command
    info_parser = subparsers.add_parser('info', help='Show crystal information')
    info_parser.add_argument('-v', '--verbose', action='store_true',
                             help='Show additional details')

    # List command
    list_parser = subparsers.add_parser('list', help='List selected crystals')
    list_parser.add_argument('--crystals', '-c', nargs='+', metavar='SPEC',
                             help='Crystal specifiers')

    # Resolve command
    resolve_parser = subparsers.add_parser('resolve',
                                           help='Show files that would be selected')
    resolve_parser.add_argument('--crystals', '-c', nargs='+', metavar='SPEC',
                                help='Crystal specifiers')
    resolve_parser.add_argument('--template', '-t', choices=list(TEMPLATES.keys()),
                                help='File template preset')
    resolve_parser.add_argument('--files', '-f', metavar='PATTERNS',
                                help='Custom file patterns (comma-separated)')
    resolve_parser.add_argument('-v', '--verbose', action='store_true',
                                help='Show detailed file list')

    # Pull command (remote only)
    pull_parser = subparsers.add_parser('pull',
                                        help='Download files from remote server')
    pull_parser.add_argument('--crystals', '-c', nargs='+', metavar='SPEC',
                             help='Crystal specifiers')
    pull_parser.add_argument('--template', '-t', choices=list(TEMPLATES.keys()),
                             help='File template preset')
    pull_parser.add_argument('--files', '-f', metavar='PATTERNS',
                             help='Custom file patterns (comma-separated)')
    pull_parser.add_argument('--dest', '-d', required=True, metavar='PATH',
                             help='Local destination directory')
    pull_parser.add_argument('--manifest', '-m', metavar='PATH',
                             help='Write transfer manifest to file')
    pull_parser.add_argument('--from-manifest', metavar='PATH',
                             help='Read files from existing manifest (skip SFTP resolution)')
    pull_parser.add_argument('--skip-existing', action='store_true',
                             help='Skip files that already exist with same size')
    pull_parser.add_argument('--dry-run', '-n', action='store_true',
                             help='Show what would be downloaded without downloading')
    pull_parser.add_argument('--yes', '-y', action='store_true',
                             help='Skip confirmation for large downloads')
    pull_parser.add_argument('--quiet', '-q', action='store_true',
                             help='Suppress progress output')

    # Manifest command (remote)
    manifest_parser = subparsers.add_parser('manifest',
                                            help='Generate transfer manifest')
    manifest_parser.add_argument('--crystals', '-c', nargs='+', metavar='SPEC',
                                 help='Crystal specifiers')
    manifest_parser.add_argument('--template', '-t', choices=list(TEMPLATES.keys()),
                                 help='File template preset')
    manifest_parser.add_argument('--files', '-f', metavar='PATTERNS',
                                 help='Custom file patterns (comma-separated)')
    manifest_parser.add_argument('--output', '-o', metavar='PREFIX',
                                 help='Output file prefix')
    manifest_parser.add_argument('--rsync', metavar='DEST',
                                 help='Also generate rsync command for destination')

    # Archive command (local only)
    archive_parser = subparsers.add_parser('archive',
                                           help='Create archive (local only)')
    archive_parser.add_argument('--crystals', '-c', nargs='+', metavar='SPEC',
                                help='Crystal specifiers')
    archive_parser.add_argument('--template', '-t', choices=list(TEMPLATES.keys()),
                                help='File template preset')
    archive_parser.add_argument('--files', '-f', metavar='PATTERNS',
                                help='Custom file patterns (comma-separated)')
    archive_parser.add_argument('--output', '-o', required=True, metavar='PATH',
                                help='Output archive path')
    archive_parser.add_argument('--archive-name', metavar='NAME',
                                help='Root directory name in archive')
    archive_parser.add_argument('--quiet', '-q', action='store_true',
                                help='Suppress progress output')

    return parser


def main(argv=None):
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        return 1

    # Parse source to determine local vs remote
    is_remote, path, username, host = parse_source(args.source)

    try:
        if is_remote:
            # Remote commands
            if args.command == 'info':
                cmd_info_remote(args, path, username, host)
            elif args.command == 'list':
                cmd_list_remote(args, path, username, host)
            elif args.command == 'resolve':
                cmd_resolve_remote(args, path, username, host)
            elif args.command == 'pull':
                cmd_pull_remote(args, path, username, host)
            elif args.command == 'manifest':
                cmd_manifest_remote(args, path, username, host)
            elif args.command == 'archive':
                print("Error: 'archive' command only works with local sources.",
                      file=sys.stderr)
                print("Use 'pull' to download files first, then archive locally.",
                      file=sys.stderr)
                return 1
            else:
                parser.print_help()
                return 1
        else:
            # Local commands
            if args.command == 'info':
                cmd_info_local(args, path)
            elif args.command == 'list':
                cmd_list_local(args, path)
            elif args.command == 'resolve':
                cmd_resolve_local(args, path)
            elif args.command == 'archive':
                cmd_archive_local(args, path)
            elif args.command in ('pull', 'manifest'):
                print("Error: '{}' command only works with remote sources.".format(
                    args.command), file=sys.stderr)
                print("Use user@host:/path format for remote access.", file=sys.stderr)
                return 1
            else:
                parser.print_help()
                return 1

    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    except ImportError as e:
        if 'paramiko' in str(e):
            print("Error: paramiko is required for remote operations.", file=sys.stderr)
            print("Install with: pip install paramiko", file=sys.stderr)
            return 1
        raise
    except Exception as e:
        print("Error: {}".format(str(e)), file=sys.stderr)
        if os.environ.get('DEBUG'):
            raise
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
