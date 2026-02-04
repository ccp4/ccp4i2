"""
Archiver - create tar.gz and zip archives from resolved files.

All symlinks are dereferenced (actual file content is archived).
"""

import os
import tarfile
import zipfile
import shutil
from typing import Dict, List, Optional, Callable

from .resolver import ResolvedFile, format_size


class ArchiveProgress:
    """Track progress of archive creation."""

    def __init__(self, total_files, total_bytes, callback=None):
        # type: (int, int, Optional[Callable]) -> None
        self.total_files = total_files
        self.total_bytes = total_bytes
        self.files_done = 0
        self.bytes_done = 0
        self.callback = callback

    def update(self, file_path, file_size):
        # type: (str, int) -> None
        """Update progress after processing a file."""
        self.files_done += 1
        self.bytes_done += file_size

        if self.callback:
            self.callback(self)

    @property
    def percent_files(self):
        # type: () -> float
        if self.total_files == 0:
            return 100.0
        return 100.0 * self.files_done / self.total_files

    @property
    def percent_bytes(self):
        # type: () -> float
        if self.total_bytes == 0:
            return 100.0
        return 100.0 * self.bytes_done / self.total_bytes


def create_tar_archive(resolved_files, output_path, compression='gz',
                       progress_callback=None, archive_name=None):
    # type: (Dict[int, List[ResolvedFile]], str, str, Optional[Callable], Optional[str]) -> str
    """
    Create a tar archive from resolved files.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output archive
        compression: 'gz' for gzip, 'bz2' for bzip2, '' for no compression
        progress_callback: Optional callback for progress updates
        archive_name: Optional name for the root directory in the archive

    Returns:
        Path to the created archive
    """
    # Ensure proper extension
    if compression == 'gz' and not output_path.endswith('.tar.gz'):
        output_path = output_path + '.tar.gz'
    elif compression == 'bz2' and not output_path.endswith('.tar.bz2'):
        output_path = output_path + '.tar.bz2'
    elif compression == '' and not output_path.endswith('.tar'):
        output_path = output_path + '.tar'

    # Calculate totals for progress
    total_files = sum(len(files) for files in resolved_files.values())
    total_bytes = 0
    for files in resolved_files.values():
        for f in files:
            try:
                total_bytes += os.path.getsize(f.source_path)
            except OSError:
                pass

    progress = ArchiveProgress(total_files, total_bytes, progress_callback)

    # Determine tar mode
    if compression == 'gz':
        mode = 'w:gz'
    elif compression == 'bz2':
        mode = 'w:bz2'
    else:
        mode = 'w'

    with tarfile.open(output_path, mode) as tar:
        seen_paths = set()

        for index in sorted(resolved_files.keys()):
            files = resolved_files[index]

            for f in files:
                # Determine archive path
                if archive_name:
                    arcname = os.path.join(archive_name, f.output_path)
                else:
                    arcname = f.output_path

                # Skip duplicates
                if arcname in seen_paths:
                    continue
                seen_paths.add(arcname)

                # Add file to archive
                try:
                    file_size = os.path.getsize(f.source_path)
                    tar.add(f.source_path, arcname=arcname)
                    progress.update(f.source_path, file_size)
                except OSError as e:
                    # Log error but continue
                    print("Warning: Could not add {}: {}".format(
                        f.source_path, str(e)))

    return output_path


def create_zip_archive(resolved_files, output_path, compression=zipfile.ZIP_DEFLATED,
                       progress_callback=None, archive_name=None):
    # type: (Dict[int, List[ResolvedFile]], str, int, Optional[Callable], Optional[str]) -> str
    """
    Create a zip archive from resolved files.

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_path: Path for the output archive
        compression: zipfile compression type (ZIP_DEFLATED, ZIP_STORED, etc.)
        progress_callback: Optional callback for progress updates
        archive_name: Optional name for the root directory in the archive

    Returns:
        Path to the created archive
    """
    # Ensure proper extension
    if not output_path.endswith('.zip'):
        output_path = output_path + '.zip'

    # Calculate totals for progress
    total_files = sum(len(files) for files in resolved_files.values())
    total_bytes = 0
    for files in resolved_files.values():
        for f in files:
            try:
                total_bytes += os.path.getsize(f.source_path)
            except OSError:
                pass

    progress = ArchiveProgress(total_files, total_bytes, progress_callback)

    with zipfile.ZipFile(output_path, 'w', compression=compression) as zf:
        seen_paths = set()

        for index in sorted(resolved_files.keys()):
            files = resolved_files[index]

            for f in files:
                # Determine archive path
                if archive_name:
                    arcname = os.path.join(archive_name, f.output_path)
                else:
                    arcname = f.output_path

                # Skip duplicates
                if arcname in seen_paths:
                    continue
                seen_paths.add(arcname)

                # Add file to archive
                try:
                    file_size = os.path.getsize(f.source_path)
                    zf.write(f.source_path, arcname=arcname)
                    progress.update(f.source_path, file_size)
                except OSError as e:
                    print("Warning: Could not add {}: {}".format(
                        f.source_path, str(e)))

    return output_path


def copy_to_directory(resolved_files, output_dir, progress_callback=None):
    # type: (Dict[int, List[ResolvedFile]], str, Optional[Callable]) -> str
    """
    Copy resolved files to a directory (for subsequent rsync/transfer).

    Args:
        resolved_files: Dict mapping crystal index to list of ResolvedFile
        output_dir: Base directory for output
        progress_callback: Optional callback for progress updates

    Returns:
        Path to the output directory
    """
    # Calculate totals for progress
    total_files = sum(len(files) for files in resolved_files.values())
    total_bytes = 0
    for files in resolved_files.values():
        for f in files:
            try:
                total_bytes += os.path.getsize(f.source_path)
            except OSError:
                pass

    progress = ArchiveProgress(total_files, total_bytes, progress_callback)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    seen_paths = set()

    for index in sorted(resolved_files.keys()):
        files = resolved_files[index]

        for f in files:
            # Skip duplicates
            if f.output_path in seen_paths:
                continue
            seen_paths.add(f.output_path)

            # Create full output path
            dest_path = os.path.join(output_dir, f.output_path)

            # Create parent directories
            dest_dir = os.path.dirname(dest_path)
            os.makedirs(dest_dir, exist_ok=True)

            # Copy file
            try:
                file_size = os.path.getsize(f.source_path)
                shutil.copy2(f.source_path, dest_path)
                progress.update(f.source_path, file_size)
            except OSError as e:
                print("Warning: Could not copy {}: {}".format(
                    f.source_path, str(e)))

    return output_dir


def estimate_archive_size(resolved_files):
    # type: (Dict[int, List[ResolvedFile]]) -> Dict[str, int]
    """
    Estimate the size of archives in various formats.

    Returns dict with keys: 'raw', 'gzip_estimate', 'file_count'
    """
    total_raw = 0
    file_count = 0
    seen = set()

    for files in resolved_files.values():
        for f in files:
            if f.source_path not in seen:
                seen.add(f.source_path)
                try:
                    total_raw += os.path.getsize(f.source_path)
                    file_count += 1
                except OSError:
                    pass

    # Rough estimate: crystallography data (binary MTZ, text logs) typically
    # compresses to about 60-70% of original size
    gzip_estimate = int(total_raw * 0.65)

    return {
        'raw': total_raw,
        'gzip_estimate': gzip_estimate,
        'file_count': file_count,
    }
