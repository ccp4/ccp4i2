"""
Remote file resolver - resolves file templates using SFTP.

Handles:
- Remote glob pattern expansion
- Symlink resolution over SFTP
- Sibling file discovery from autoprocessing directories
"""

import posixpath
import fnmatch
from typing import List, Dict, Any, Optional, Set, Callable

from .remote import SFTPScanner, RemoteCrystalScanner, format_size


class RemoteResolvedFile:
    """Represents a resolved file on the remote server."""

    def __init__(self,
                 remote_path,       # type: str
                 local_path,        # type: str
                 original_pattern,  # type: str
                 size=0,            # type: int
                 is_sibling=False   # type: bool
                 ):
        # type: (...) -> None
        self.remote_path = remote_path    # Path on remote server (resolved)
        self.local_path = local_path      # Path for local storage
        self.original_pattern = original_pattern
        self.size = size
        self.is_sibling = is_sibling

    def __repr__(self):
        return "RemoteResolvedFile({!r} -> {!r}, {})".format(
            self.remote_path, self.local_path, format_size(self.size))


class RemoteFileResolver:
    """
    Resolves file templates to actual files on a remote server via SFTP.
    """

    def __init__(self, sftp, crystal_path, crystal_name):
        # type: (SFTPScanner, str, str) -> None
        """
        Initialize the resolver.

        Args:
            sftp: Connected SFTPScanner instance
            crystal_path: Full remote path to the crystal directory
            crystal_name: Name of the crystal directory
        """
        self.sftp = sftp
        self.crystal_path = crystal_path
        self.crystal_name = crystal_name

    def _remote_glob(self, pattern):
        # type: (str) -> List[str]
        """
        Expand a glob pattern on the remote server.

        Handles patterns like:
        - 'dimple/dimple/final.pdb'  -> exact path
        - 'compound/*.cif'           -> wildcard
        - '*.mtz'                    -> wildcard in crystal root
        """
        full_pattern = posixpath.join(self.crystal_path, pattern)

        # Split into directory parts
        parts = pattern.split('/')
        current_paths = [self.crystal_path]

        for part in parts:
            next_paths = []

            for current in current_paths:
                if '*' in part or '?' in part:
                    # Glob this level
                    try:
                        entries = self.sftp.listdir(current)
                        for entry in entries:
                            if fnmatch.fnmatch(entry, part):
                                full = posixpath.join(current, entry)
                                next_paths.append(full)
                    except IOError:
                        pass
                else:
                    # Exact path
                    full = posixpath.join(current, part)
                    if self.sftp.exists(full):
                        next_paths.append(full)

            current_paths = next_paths

        return current_paths

    def _get_sibling_files(self, resolved_path, sibling_patterns):
        # type: (str, List[str]) -> List[tuple]
        """
        Get sibling files from the same directory as a resolved symlink target.

        Returns list of (remote_path, filename, size) tuples.
        """
        siblings = []
        parent_dir = posixpath.dirname(resolved_path)

        for pattern in sibling_patterns:
            try:
                entries = self.sftp.listdir(parent_dir)
                for entry in entries:
                    if fnmatch.fnmatch(entry, pattern):
                        full_path = posixpath.join(parent_dir, entry)
                        if self.sftp.isfile(full_path):
                            size = self.sftp.get_file_size(full_path)
                            siblings.append((full_path, entry, size))
            except IOError:
                pass

        return siblings

    def resolve_spec(self, spec):
        # type: (Dict[str, Any]) -> List[RemoteResolvedFile]
        """
        Resolve a single file spec to actual remote files.

        Args:
            spec: A file spec dictionary from config

        Returns:
            List of RemoteResolvedFile objects
        """
        pattern = spec['pattern']
        output_name = spec.get('output_name')
        dereference = spec.get('dereference', True)
        grab_siblings = spec.get('grab_siblings')
        required = spec.get('required', False)

        # Expand the glob pattern
        matches = self._remote_glob(pattern)

        if not matches and required:
            raise FileNotFoundError(
                "Required file not found: {} in {}".format(
                    pattern, self.crystal_name))

        results = []
        seen_remotes = set()  # type: Set[str]

        for match_path in matches:
            # Skip directories (unless they're symlinks to files)
            try:
                if self.sftp.isdir(match_path) and not self.sftp.islink(match_path):
                    continue
            except IOError:
                continue

            # Determine the actual remote file path
            if dereference and self.sftp.islink(match_path):
                remote_path = self.sftp.resolve_symlink_fully(match_path)
            else:
                remote_path = match_path

            # Skip if already seen
            if remote_path in seen_remotes:
                continue
            seen_remotes.add(remote_path)

            # Verify the file exists
            if not self.sftp.isfile(remote_path):
                if required:
                    raise FileNotFoundError(
                        "Symlink target not found: {} -> {}".format(
                            match_path, remote_path))
                continue

            # Get file size
            size = self.sftp.get_file_size(remote_path)

            # Determine local output path
            if output_name:
                local_path = posixpath.join(self.crystal_name, output_name)
            else:
                # Preserve relative path structure
                rel_path = match_path[len(self.crystal_path):].lstrip('/')
                local_path = posixpath.join(self.crystal_name, rel_path)

            results.append(RemoteResolvedFile(
                remote_path=remote_path,
                local_path=local_path,
                original_pattern=pattern,
                size=size,
            ))

            # Handle sibling files
            if grab_siblings and self.sftp.islink(match_path):
                siblings = self._get_sibling_files(remote_path, grab_siblings)
                for sib_path, sib_name, sib_size in siblings:
                    if sib_path in seen_remotes:
                        continue
                    seen_remotes.add(sib_path)

                    sib_local = posixpath.join(
                        self.crystal_name, 'autoprocessing', sib_name)

                    results.append(RemoteResolvedFile(
                        remote_path=sib_path,
                        local_path=sib_local,
                        original_pattern=pattern,
                        size=sib_size,
                        is_sibling=True,
                    ))

        return results

    def resolve_template(self, template):
        # type: (List[Dict[str, Any]]) -> List[RemoteResolvedFile]
        """
        Resolve a full file template to actual remote files.

        Args:
            template: List of file spec dictionaries

        Returns:
            List of RemoteResolvedFile objects
        """
        all_results = []
        seen_locals = set()  # type: Set[str]

        for spec in template:
            try:
                files = self.resolve_spec(spec)
                for f in files:
                    if f.local_path not in seen_locals:
                        all_results.append(f)
                        seen_locals.add(f.local_path)
            except FileNotFoundError:
                if spec.get('required', False):
                    raise

        return all_results


class ResolveResult:
    """Result of resolving files across multiple crystals."""

    def __init__(self):
        self.resolved = {}  # type: Dict[int, List[RemoteResolvedFile]]
        self.skipped = []   # type: List[tuple]  # (index, reason)
        self.empty = []     # type: List[int]    # crystals with no files

    @property
    def crystals_with_files(self):
        # type: () -> int
        return len([idx for idx, files in self.resolved.items() if files])


def resolve_crystals_remote(sftp, scanner, indices, template, progress_callback=None):
    # type: (SFTPScanner, RemoteCrystalScanner, List[int], List[Dict[str, Any]], Optional[Callable]) -> Dict[int, List[RemoteResolvedFile]]
    """
    Resolve files for multiple crystals via SFTP.

    Args:
        sftp: Connected SFTPScanner
        scanner: RemoteCrystalScanner instance
        indices: List of crystal indices to process
        template: File template to apply
        progress_callback: Optional callback(message) for progress updates

    Returns:
        Dict mapping crystal index to list of RemoteResolvedFile objects
        (crystals with no matching files are included with empty lists)
    """
    results = {}
    skipped = []
    empty = []
    total = len(indices)

    for i, index in enumerate(indices):
        if progress_callback:
            progress_callback("Resolving x{:04d} ({}/{})".format(index, i + 1, total))

        crystal_path = scanner.get_crystal_path(index)
        crystal_name = scanner.get_crystal_dir_name(index)

        if crystal_path is None:
            skipped.append((index, "crystal directory not found"))
            continue

        resolver = RemoteFileResolver(sftp, crystal_path, crystal_name)
        try:
            files = resolver.resolve_template(template)
            if files:
                results[index] = files
            else:
                empty.append(index)
                # Include empty list so we know it was processed
                results[index] = []
        except FileNotFoundError as e:
            # Skip crystals with missing required files
            skipped.append((index, str(e)))
            continue

    # Report skipped/empty crystals
    if skipped and progress_callback:
        progress_callback("Skipped {} crystals with missing required files".format(len(skipped)))
    if empty and progress_callback:
        progress_callback("{} crystals had no matching files".format(len(empty)))

    return results


def calculate_total_size_remote(resolved_files):
    # type: (Dict[int, List[RemoteResolvedFile]]) -> int
    """Calculate total size of all resolved remote files in bytes."""
    total = 0
    seen = set()  # type: Set[str]

    for index, files in resolved_files.items():
        for f in files:
            if f.remote_path not in seen:
                seen.add(f.remote_path)
                total += f.size

    return total


def count_files_remote(resolved_files):
    # type: (Dict[int, List[RemoteResolvedFile]]) -> int
    """Count unique files in resolved set."""
    seen = set()
    for files in resolved_files.values():
        for f in files:
            seen.add(f.remote_path)
    return len(seen)
