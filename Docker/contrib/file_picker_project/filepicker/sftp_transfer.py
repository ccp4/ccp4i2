"""
SFTP-based file transfer for pulling crystallography data.

Handles:
- Direct SFTP download with progress reporting
- Parallel download option
- Resume capability for interrupted transfers
"""

import os
import sys
import time
from typing import Dict, List, Optional, Callable

from .remote import SFTPScanner, format_size
from .remote_resolver import RemoteResolvedFile, calculate_total_size_remote


class TransferProgress:
    """Track progress of file transfers."""

    def __init__(self, total_files, total_bytes):
        # type: (int, int) -> None
        self.total_files = total_files
        self.total_bytes = total_bytes
        self.files_done = 0
        self.bytes_done = 0
        self.current_file = ""
        self.current_file_bytes = 0
        self.current_file_total = 0
        self.start_time = time.time()
        self.errors = []  # type: List[str]

    def start_file(self, filename, size):
        # type: (str, int) -> None
        """Mark the start of a new file transfer."""
        self.current_file = filename
        self.current_file_bytes = 0
        self.current_file_total = size

    def update_file(self, bytes_transferred, total_bytes):
        # type: (int, int) -> None
        """Update progress for current file."""
        self.current_file_bytes = bytes_transferred
        self.current_file_total = total_bytes

    def complete_file(self, size):
        # type: (int) -> None
        """Mark a file as completed."""
        self.files_done += 1
        self.bytes_done += size
        self.current_file = ""
        self.current_file_bytes = 0
        self.current_file_total = 0

    def add_error(self, message):
        # type: (str) -> None
        """Record an error."""
        self.errors.append(message)

    @property
    def percent_bytes(self):
        # type: () -> float
        if self.total_bytes == 0:
            return 100.0
        return 100.0 * self.bytes_done / self.total_bytes

    @property
    def elapsed_seconds(self):
        # type: () -> float
        return time.time() - self.start_time

    @property
    def bytes_per_second(self):
        # type: () -> float
        elapsed = self.elapsed_seconds
        if elapsed == 0:
            return 0
        return self.bytes_done / elapsed

    @property
    def eta_seconds(self):
        # type: () -> float
        bps = self.bytes_per_second
        if bps == 0:
            return 0
        remaining = self.total_bytes - self.bytes_done
        return remaining / bps


def print_transfer_progress(progress):
    # type: (TransferProgress) -> None
    """Print a progress bar to stderr."""
    bar_width = 30
    filled = int(bar_width * progress.percent_bytes / 100)
    bar = '=' * filled + '-' * (bar_width - filled)

    # Speed and ETA
    speed = format_size(int(progress.bytes_per_second)) + "/s"
    eta = int(progress.eta_seconds)
    if eta > 3600:
        eta_str = "{:.0f}h{:.0f}m".format(eta // 3600, (eta % 3600) // 60)
    elif eta > 60:
        eta_str = "{:.0f}m{:.0f}s".format(eta // 60, eta % 60)
    else:
        eta_str = "{:.0f}s".format(eta)

    # Current file (truncated)
    current = progress.current_file
    if len(current) > 30:
        current = "..." + current[-27:]

    sys.stderr.write('\r[{}] {:.1f}% | {}/{} files | {} | ETA {} | {}   '.format(
        bar,
        progress.percent_bytes,
        progress.files_done,
        progress.total_files,
        speed,
        eta_str,
        current,
    ))
    sys.stderr.flush()


def download_files(sftp, resolved_files, output_dir, progress_callback=None, skip_existing=False):
    # type: (SFTPScanner, Dict[int, List[RemoteResolvedFile]], str, Optional[Callable], bool) -> TransferProgress
    """
    Download resolved files from remote server via SFTP.

    Args:
        sftp: Connected SFTPScanner instance
        resolved_files: Dict mapping crystal index to resolved files
        output_dir: Local directory to store files
        progress_callback: Optional callback(TransferProgress) for updates
        skip_existing: If True, skip files that already exist locally with same size

    Returns:
        TransferProgress with final statistics
    """
    # Calculate totals
    total_bytes = calculate_total_size_remote(resolved_files)
    total_files = sum(len(files) for files in resolved_files.values())

    progress = TransferProgress(total_files, total_bytes)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Track seen files to avoid duplicates
    seen_remotes = set()

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]

        for f in files:
            # Skip duplicates
            if f.remote_path in seen_remotes:
                continue
            seen_remotes.add(f.remote_path)

            # Build local path
            local_path = os.path.join(output_dir, f.local_path)

            # Check if file exists and matches size
            if skip_existing and os.path.exists(local_path):
                local_size = os.path.getsize(local_path)
                if local_size == f.size:
                    progress.complete_file(f.size)
                    if progress_callback:
                        progress_callback(progress)
                    continue

            # Ensure parent directory exists
            local_dir = os.path.dirname(local_path)
            if local_dir and not os.path.exists(local_dir):
                os.makedirs(local_dir)

            # Start file transfer
            progress.start_file(f.local_path, f.size)
            if progress_callback:
                progress_callback(progress)

            try:
                # Create a callback for per-file progress
                def file_progress(transferred, total):
                    progress.update_file(transferred, total)
                    if progress_callback:
                        progress_callback(progress)

                sftp.download_file(f.remote_path, local_path, callback=file_progress)
                progress.complete_file(f.size)

            except Exception as e:
                progress.add_error("Failed to download {}: {}".format(
                    f.remote_path, str(e)))

            if progress_callback:
                progress_callback(progress)

    return progress


def generate_rsync_pull_command(resolved_files, host, username, output_dir, file_list_path=None):
    # type: (Dict[int, List[RemoteResolvedFile]], str, str, str, Optional[str]) -> str
    """
    Generate an rsync command to pull files from remote.

    Note: This uses the resolved (dereferenced) paths, so rsync will
    pull the actual files, not symlinks.

    Args:
        resolved_files: Dict mapping crystal index to resolved files
        host: Remote hostname
        username: SSH username
        output_dir: Local destination directory
        file_list_path: Path to write file list (for --files-from)

    Returns:
        rsync command string
    """
    # Collect unique remote paths
    remote_paths = set()
    for files in resolved_files.values():
        for f in files:
            remote_paths.add(f.remote_path)

    if file_list_path:
        # Write file list
        with open(file_list_path, 'w') as fh:
            for path in sorted(remote_paths):
                fh.write(path + '\n')

        cmd = 'rsync -avP --progress --files-from="{}" {}@{}: "{}"'.format(
            file_list_path, username, host, output_dir)
    else:
        # Include paths directly (may be very long)
        paths_str = ' '.join('"{}"'.format(p) for p in sorted(remote_paths))
        cmd = 'rsync -avP --progress {}@{}:{} "{}"'.format(
            username, host, paths_str, output_dir)

    return cmd


def generate_sftp_batch_script(resolved_files, output_dir, batch_file_path):
    # type: (Dict[int, List[RemoteResolvedFile]], str, str) -> str
    """
    Generate an SFTP batch file for downloading files.

    This can be used with: sftp -b batch_file user@host

    Args:
        resolved_files: Dict mapping crystal index to resolved files
        output_dir: Local destination directory
        batch_file_path: Path to write the batch file

    Returns:
        Path to the batch file
    """
    lines = []

    # Track directories we need to create locally
    local_dirs = set()

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]

        for f in files:
            local_path = os.path.join(output_dir, f.local_path)
            local_dir = os.path.dirname(local_path)

            # Track directory for later creation
            if local_dir:
                local_dirs.add(local_dir)

            # SFTP get command
            lines.append('get "{}" "{}"'.format(f.remote_path, local_path))

    # Write batch file
    with open(batch_file_path, 'w') as fh:
        # Header comment
        fh.write('# SFTP batch file generated by filepicker\n')
        fh.write('# Run with: sftp -b {} user@host\n'.format(batch_file_path))
        fh.write('# First create local directories:\n')
        for d in sorted(local_dirs):
            fh.write('#   mkdir -p "{}"\n'.format(d))
        fh.write('\n')

        # Get commands
        for line in lines:
            fh.write(line + '\n')

    return batch_file_path


def write_transfer_manifest(resolved_files, manifest_path, include_sizes=True):
    # type: (Dict[int, List[RemoteResolvedFile]], str, bool) -> str
    """
    Write a manifest file listing all files to transfer.

    Format: remote_path -> local_path [size]

    Args:
        resolved_files: Dict mapping crystal index to resolved files
        manifest_path: Path to write manifest
        include_sizes: Include file sizes in manifest

    Returns:
        Path to manifest file
    """
    with open(manifest_path, 'w') as fh:
        fh.write("# Transfer manifest generated by filepicker\n")
        fh.write("# Format: remote_path -> local_path [size]\n\n")

        total_size = 0
        file_count = 0

        for index in sorted(resolved_files.keys()):
            fh.write("# Crystal x{:04d}\n".format(index))
            files = resolved_files[index]

            for f in files:
                if include_sizes:
                    fh.write("{} -> {} [{}]\n".format(
                        f.remote_path, f.local_path, format_size(f.size)))
                else:
                    fh.write("{} -> {}\n".format(f.remote_path, f.local_path))

                total_size += f.size
                file_count += 1

            fh.write("\n")

        fh.write("# Total: {} files, {}\n".format(file_count, format_size(total_size)))

    return manifest_path


class ManifestEntry:
    """An entry parsed from a transfer manifest."""

    def __init__(self, remote_path, local_path, size=0):
        # type: (str, str, int) -> None
        self.remote_path = remote_path
        self.local_path = local_path
        self.size = size


def parse_size_string(size_str):
    # type: (str) -> int
    """Parse a human-readable size string back to bytes."""
    size_str = size_str.strip()
    if not size_str:
        return 0

    # Handle "123 KB", "1.5 MB", "2 GB", etc.
    import re
    match = re.match(r'^([\d.]+)\s*([KMGT]?B?)$', size_str, re.IGNORECASE)
    if not match:
        return 0

    value = float(match.group(1))
    unit = match.group(2).upper()

    multipliers = {
        'B': 1,
        '': 1,
        'KB': 1024,
        'MB': 1024 * 1024,
        'GB': 1024 * 1024 * 1024,
        'TB': 1024 * 1024 * 1024 * 1024,
    }

    return int(value * multipliers.get(unit, 1))


def parse_manifest(manifest_path):
    # type: (str) -> List[ManifestEntry]
    """
    Parse a transfer manifest file.

    Format: remote_path -> local_path [size]

    Args:
        manifest_path: Path to manifest file

    Returns:
        List of ManifestEntry objects
    """
    import re

    entries = []
    # Pattern: path -> path [optional size]
    line_pattern = re.compile(r'^(.+?)\s+->\s+(.+?)(?:\s+\[([^\]]+)\])?$')

    with open(manifest_path, 'r') as fh:
        for line in fh:
            line = line.strip()

            # Skip comments and empty lines
            if not line or line.startswith('#'):
                continue

            match = line_pattern.match(line)
            if match:
                remote_path = match.group(1).strip()
                local_path = match.group(2).strip()
                size_str = match.group(3) or ''
                size = parse_size_string(size_str)

                entries.append(ManifestEntry(remote_path, local_path, size))

    return entries


def download_from_manifest(sftp, entries, output_dir, progress_callback=None, skip_existing=False):
    # type: (SFTPScanner, List[ManifestEntry], str, Optional[Callable], bool) -> TransferProgress
    """
    Download files using a parsed manifest.

    Args:
        sftp: Connected SFTPScanner instance
        entries: List of ManifestEntry objects from parse_manifest()
        output_dir: Local directory to store files
        progress_callback: Optional callback(TransferProgress) for updates
        skip_existing: If True, skip files that already exist locally with same size

    Returns:
        TransferProgress with final statistics
    """
    # Calculate totals
    total_bytes = sum(e.size for e in entries)
    total_files = len(entries)

    progress = TransferProgress(total_files, total_bytes)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Track seen files to avoid duplicates
    seen_remotes = set()

    for entry in entries:
        # Skip duplicates
        if entry.remote_path in seen_remotes:
            continue
        seen_remotes.add(entry.remote_path)

        # Build local path
        local_path = os.path.join(output_dir, entry.local_path)

        # Check if file exists and matches size
        if skip_existing and os.path.exists(local_path):
            local_size = os.path.getsize(local_path)
            # If manifest has size info, compare; otherwise skip if file exists
            if entry.size == 0 or local_size == entry.size:
                progress.complete_file(entry.size)
                if progress_callback:
                    progress_callback(progress)
                continue

        # Ensure parent directory exists
        local_dir = os.path.dirname(local_path)
        if local_dir and not os.path.exists(local_dir):
            os.makedirs(local_dir)

        # Start file transfer
        progress.start_file(entry.local_path, entry.size)
        if progress_callback:
            progress_callback(progress)

        try:
            # Create a callback for per-file progress
            def file_progress(transferred, total):
                progress.update_file(transferred, total)
                if progress_callback:
                    progress_callback(progress)

            # Get actual file size from remote if not in manifest
            actual_size = entry.size
            if actual_size == 0:
                try:
                    actual_size = sftp.get_file_size(entry.remote_path)
                except Exception:
                    pass

            sftp.download_file(entry.remote_path, local_path, callback=file_progress)
            progress.complete_file(actual_size)

        except Exception as e:
            progress.add_error("Failed to download {}: {}".format(
                entry.remote_path, str(e)))

        if progress_callback:
            progress_callback(progress)

    return progress
