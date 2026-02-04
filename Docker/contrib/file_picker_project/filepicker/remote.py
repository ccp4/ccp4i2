"""
Remote SFTP scanner for discovering and resolving files on Diamond Light Source.

Uses paramiko for SFTP operations including:
- Directory listing
- Symlink resolution
- File transfer
"""

import os
import stat as stat_module
import posixpath
import getpass
from typing import List, Dict, Optional, Tuple, Callable, Any

try:
    import paramiko
    HAS_PARAMIKO = True
except ImportError:
    HAS_PARAMIKO = False
    paramiko = None


class RemoteConnectionError(Exception):
    """Raised when connection to remote host fails."""
    pass


class SFTPScanner:
    """
    SFTP-based scanner for remote crystallography data directories.
    """

    def __init__(self, host, username, password=None, key_filename=None, port=22):
        # type: (str, str, Optional[str], Optional[str], int) -> None
        """
        Initialize SFTP connection.

        Args:
            host: Remote hostname (e.g., 'nx.diamond.ac.uk')
            username: SSH username (FedID for Diamond)
            password: SSH password (will prompt if None and no key)
            key_filename: Path to SSH private key
            port: SSH port (default 22)
        """
        if not HAS_PARAMIKO:
            raise ImportError(
                "paramiko is required for remote operations. "
                "Install with: pip install paramiko"
            )

        self.host = host
        self.username = username
        self.port = port

        self._ssh = None  # type: Optional[paramiko.SSHClient]
        self._sftp = None  # type: Optional[paramiko.SFTPClient]

        self._connect(password, key_filename)

    def _connect(self, password, key_filename):
        # type: (Optional[str], Optional[str]) -> None
        """Establish SSH/SFTP connection."""
        self._ssh = paramiko.SSHClient()
        self._ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        try:
            connect_kwargs = {
                'hostname': self.host,
                'port': self.port,
                'username': self.username,
            }

            if key_filename:
                connect_kwargs['key_filename'] = key_filename
            elif password:
                connect_kwargs['password'] = password
            else:
                # Try agent first, then prompt for password
                connect_kwargs['allow_agent'] = True
                connect_kwargs['look_for_keys'] = True

            try:
                self._ssh.connect(**connect_kwargs)
            except paramiko.AuthenticationException:
                # If agent/keys failed, prompt for password
                if not password:
                    password = getpass.getpass(
                        "Password for {}@{}: ".format(self.username, self.host))
                    connect_kwargs['password'] = password
                    connect_kwargs['allow_agent'] = False
                    connect_kwargs['look_for_keys'] = False
                    self._ssh.connect(**connect_kwargs)
                else:
                    raise

            self._sftp = self._ssh.open_sftp()

        except Exception as e:
            raise RemoteConnectionError(
                "Failed to connect to {}@{}: {}".format(
                    self.username, self.host, str(e)))

    def close(self):
        # type: () -> None
        """Close the SFTP and SSH connections."""
        if self._sftp:
            self._sftp.close()
            self._sftp = None
        if self._ssh:
            self._ssh.close()
            self._ssh = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def listdir(self, path):
        # type: (str) -> List[str]
        """List directory contents."""
        return self._sftp.listdir(path)

    def listdir_attr(self, path):
        # type: (str) -> List[paramiko.SFTPAttributes]
        """List directory contents with attributes."""
        return self._sftp.listdir_attr(path)

    def stat(self, path):
        # type: (str) -> paramiko.SFTPAttributes
        """Stat a file (follows symlinks)."""
        return self._sftp.stat(path)

    def lstat(self, path):
        # type: (str) -> paramiko.SFTPAttributes
        """Stat a file (does not follow symlinks)."""
        return self._sftp.lstat(path)

    def isdir(self, path):
        # type: (str) -> bool
        """Check if path is a directory."""
        try:
            return stat_module.S_ISDIR(self._sftp.stat(path).st_mode)
        except IOError:
            return False

    def islink(self, path):
        # type: (str) -> bool
        """Check if path is a symlink."""
        try:
            return stat_module.S_ISLNK(self._sftp.lstat(path).st_mode)
        except IOError:
            return False

    def isfile(self, path):
        # type: (str) -> bool
        """Check if path is a regular file."""
        try:
            return stat_module.S_ISREG(self._sftp.stat(path).st_mode)
        except IOError:
            return False

    def exists(self, path):
        # type: (str) -> bool
        """Check if path exists."""
        try:
            self._sftp.stat(path)
            return True
        except IOError:
            return False

    def readlink(self, path):
        # type: (str) -> str
        """Read the target of a symlink."""
        return self._sftp.readlink(path)

    def realpath(self, path):
        # type: (str) -> str
        """
        Resolve a path to its real path (resolving all symlinks).

        Note: paramiko's normalize() does this, but we implement
        our own to handle relative symlinks properly.
        """
        # Use normalize for the initial resolution
        try:
            return self._sftp.normalize(path)
        except IOError:
            return path

    def resolve_symlink_fully(self, path):
        # type: (str) -> str
        """
        Fully resolve a symlink, handling relative paths.

        Returns the absolute path to the actual file.
        """
        max_depth = 10  # Prevent infinite loops
        current = path

        for _ in range(max_depth):
            if not self.islink(current):
                return current

            target = self.readlink(current)

            # Handle relative symlinks
            if not target.startswith('/'):
                parent = posixpath.dirname(current)
                target = posixpath.normpath(posixpath.join(parent, target))

            current = target

        # If we get here, too many symlinks
        return current

    def get_file_size(self, path):
        # type: (str) -> int
        """Get size of a file in bytes."""
        try:
            return self._sftp.stat(path).st_size
        except IOError:
            return 0

    def glob(self, path, pattern):
        # type: (str, str) -> List[str]
        """
        Simple glob implementation for remote paths.

        Supports * and ? wildcards in the filename part only.
        """
        import fnmatch

        # Split pattern into directory and filename parts
        if '/' in pattern:
            dir_part, file_pattern = pattern.rsplit('/', 1)
            search_path = posixpath.join(path, dir_part)
        else:
            search_path = path
            file_pattern = pattern

        # Handle ** not supported - just use single level
        if '**' in file_pattern:
            file_pattern = file_pattern.replace('**/', '').replace('**', '*')

        results = []
        try:
            for entry in self.listdir(search_path):
                if fnmatch.fnmatch(entry, file_pattern):
                    results.append(posixpath.join(search_path, entry))
        except IOError:
            pass

        return results

    def download_file(self, remote_path, local_path, callback=None):
        # type: (str, str, Optional[Callable]) -> None
        """
        Download a file from remote to local.

        Args:
            remote_path: Path on remote server
            local_path: Path on local machine
            callback: Optional callback(bytes_transferred, total_bytes)
        """
        # Ensure local directory exists
        local_dir = os.path.dirname(local_path)
        if local_dir and not os.path.exists(local_dir):
            os.makedirs(local_dir)

        self._sftp.get(remote_path, local_path, callback=callback)


class RemoteCrystalScanner:
    """
    Scans remote crystallography directories and builds a manifest.
    """

    def __init__(self, sftp, base_path, crystal_prefix=None):
        # type: (SFTPScanner, str, Optional[str]) -> None
        """
        Initialize the scanner.

        Args:
            sftp: Connected SFTPScanner instance
            base_path: Path to model_building directory on remote
            crystal_prefix: Expected crystal prefix (auto-detected if None)
        """
        self.sftp = sftp
        self.base_path = base_path
        self.crystal_prefix = crystal_prefix
        self._crystal_cache = None  # type: Optional[Dict]
        self._compound_cache = None  # type: Optional[Dict]

    def _scan_crystals(self, progress_callback=None):
        # type: (Optional[Callable]) -> Dict[int, Tuple[str, str]]
        """Scan remote directory for crystal directories."""
        import re

        if self._crystal_cache is not None:
            return self._crystal_cache

        crystal_pattern = re.compile(r'^(.+)-x(\d+)$')
        crystals = {}

        if progress_callback:
            progress_callback("Scanning crystal directories...")

        entries = self.sftp.listdir(self.base_path)

        for entry in entries:
            full_path = posixpath.join(self.base_path, entry)

            # Check if it's a directory
            if not self.sftp.isdir(full_path):
                continue

            match = crystal_pattern.match(entry)
            if match:
                prefix, index_str = match.groups()

                if self.crystal_prefix is None:
                    self.crystal_prefix = prefix
                elif prefix != self.crystal_prefix:
                    continue

                index = int(index_str)
                crystals[index] = (entry, full_path)

        if progress_callback:
            progress_callback("Found {} crystals".format(len(crystals)))

        self._crystal_cache = crystals
        return crystals

    def _scan_compounds(self, progress_callback=None):
        # type: (Optional[Callable]) -> Dict[str, List[int]]
        """Build compound ID to crystal index mapping."""
        import re

        if self._compound_cache is not None:
            return self._compound_cache

        crystals = self._scan_crystals(progress_callback)
        compound_pattern = re.compile(r'^[A-Z]{2,4}-\d+$')
        compounds = {}  # type: Dict[str, List[int]]

        for index, (dir_name, full_path) in crystals.items():
            compound_dir = posixpath.join(full_path, 'compound')

            try:
                entries = self.sftp.listdir(compound_dir)
            except IOError:
                continue

            for entry in entries:
                name_part = posixpath.splitext(entry)[0]
                if name_part.endswith('_TMP'):
                    name_part = name_part[:-4]

                if compound_pattern.match(name_part):
                    if name_part not in compounds:
                        compounds[name_part] = []
                    if index not in compounds[name_part]:
                        compounds[name_part].append(index)

        self._compound_cache = compounds
        return compounds

    def get_all_indices(self):
        # type: () -> List[int]
        """Get all crystal indices, sorted."""
        crystals = self._scan_crystals()
        return sorted(crystals.keys())

    def get_crystal_path(self, index):
        # type: (int) -> Optional[str]
        """Get remote path for a crystal by index."""
        crystals = self._scan_crystals()
        if index in crystals:
            return crystals[index][1]
        return None

    def get_crystal_dir_name(self, index):
        # type: (int) -> Optional[str]
        """Get directory name for a crystal by index."""
        crystals = self._scan_crystals()
        if index in crystals:
            return crystals[index][0]
        return None

    def summary(self):
        # type: () -> str
        """Return a summary of available crystals."""
        crystals = self._scan_crystals()
        compounds = self._scan_compounds()

        lines = [
            "Remote: {}@{}".format(self.sftp.username, self.sftp.host),
            "Path: {}".format(self.base_path),
            "Prefix: {}".format(self.crystal_prefix or "(auto-detect)"),
            "Total crystals: {}".format(len(crystals)),
        ]

        if crystals:
            indices = sorted(crystals.keys())
            lines.append("Index range: x{:04d} - x{:04d}".format(
                min(indices), max(indices)))

            expected = set(range(min(indices), max(indices) + 1))
            missing = expected - set(indices)
            if missing:
                lines.append("Missing indices: {}".format(len(missing)))

        lines.append("Unique compounds: {}".format(len(compounds)))

        return '\n'.join(lines)


def format_size(size_bytes):
    # type: (int) -> str
    """Format byte size as human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return "{:.1f} {}".format(size_bytes, unit)
        size_bytes /= 1024.0
    return "{:.1f} PB".format(size_bytes)
