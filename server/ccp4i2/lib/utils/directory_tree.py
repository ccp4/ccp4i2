"""
Beautiful directory tree visualization utility.

Provides formatted directory structure views for projects and jobs,
similar to the Unix 'tree' command but with CCP4-specific formatting
and filtering options.
"""

from pathlib import Path
from typing import Optional, List, Set


class DirectoryTree:
    """
    Creates beautiful directory tree visualizations.

    Features:
    - Customizable depth limit
    - File size display
    - File count summaries
    - Ignore patterns (e.g., *.pyc, __pycache__)
    - Color-coded output (optional)
    """

    # Unicode box drawing characters
    BRANCH = "├── "
    LAST = "└── "
    PIPE = "│   "
    SPACE = "    "

    # Default patterns to ignore
    DEFAULT_IGNORE = {
        '__pycache__',
        '*.pyc',
        '*.pyo',
        '.git',
        '.DS_Store',
        'Thumbs.db',
    }

    def __init__(
        self,
        max_depth: Optional[int] = None,
        show_sizes: bool = True,
        show_hidden: bool = False,
        ignore_patterns: Optional[Set[str]] = None
    ):
        """
        Initialize directory tree visualizer.

        Args:
            max_depth: Maximum depth to traverse (None = unlimited)
            show_sizes: Show file sizes
            show_hidden: Show hidden files/directories
            ignore_patterns: Patterns to ignore (in addition to defaults)
        """
        self.max_depth = max_depth
        self.show_sizes = show_sizes
        self.show_hidden = show_hidden
        self.ignore_patterns = self.DEFAULT_IGNORE.copy()
        if ignore_patterns:
            self.ignore_patterns.update(ignore_patterns)

        # Statistics
        self.dir_count = 0
        self.file_count = 0
        self.total_size = 0

    def should_ignore(self, path: Path) -> bool:
        """Check if a path should be ignored based on patterns."""
        name = path.name

        # Check hidden files
        if not self.show_hidden and name.startswith('.'):
            return True

        # Check ignore patterns
        for pattern in self.ignore_patterns:
            if pattern.startswith('*'):
                # Wildcard pattern
                suffix = pattern[1:]
                if name.endswith(suffix):
                    return True
            elif name == pattern:
                return True

        return False

    def format_size(self, size_bytes: int) -> str:
        """Format file size in human-readable format."""
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if size_bytes < 1024.0:
                if unit == 'B':
                    return f"{size_bytes:3.0f}{unit}"
                else:
                    return f"{size_bytes:3.1f}{unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f}PB"

    def get_entry_line(self, path: Path, is_last: bool, prefix: str) -> str:
        """
        Format a single directory entry line.

        Args:
            path: Path to format
            is_last: Whether this is the last entry in its directory
            prefix: Prefix string for indentation

        Returns:
            Formatted line string
        """
        # Choose connector
        connector = self.LAST if is_last else self.BRANCH

        # Get name and type
        name = path.name
        if path.is_dir():
            name += "/"
            size_str = ""
        else:
            try:
                size = path.stat().st_size
                self.total_size += size
                size_str = f"  ({self.format_size(size)})" if self.show_sizes else ""
            except (OSError, PermissionError):
                size_str = "  (permission denied)" if self.show_sizes else ""

        return f"{prefix}{connector}{name}{size_str}"

    def generate(
        self,
        root_path: Path,
        current_depth: int = 0,
        prefix: str = ""
    ) -> List[str]:
        """
        Generate directory tree lines recursively.

        Args:
            root_path: Root directory to visualize
            current_depth: Current recursion depth
            prefix: Current indentation prefix

        Returns:
            List of formatted lines
        """
        lines = []

        # Check depth limit
        if self.max_depth is not None and current_depth >= self.max_depth:
            return lines

        try:
            # Get entries, sorted with directories first
            entries = sorted(
                root_path.iterdir(),
                key=lambda p: (not p.is_dir(), p.name.lower())
            )

            # Filter ignored entries
            entries = [e for e in entries if not self.should_ignore(e)]

            # Process each entry
            for idx, entry in enumerate(entries):
                is_last = (idx == len(entries) - 1)

                # Add entry line
                lines.append(self.get_entry_line(entry, is_last, prefix))

                # Update statistics
                if entry.is_dir():
                    self.dir_count += 1
                else:
                    self.file_count += 1

                # Recurse into directories
                if entry.is_dir():
                    # Calculate new prefix
                    if is_last:
                        new_prefix = prefix + self.SPACE
                    else:
                        new_prefix = prefix + self.PIPE

                    # Recurse
                    sublines = self.generate(entry, current_depth + 1, new_prefix)
                    lines.extend(sublines)

        except PermissionError:
            lines.append(f"{prefix}{self.BRANCH}[Permission Denied]")

        return lines

    def visualize(self, root_path: Path, title: Optional[str] = None) -> str:
        """
        Generate complete directory tree visualization.

        Args:
            root_path: Root directory to visualize
            title: Optional title for the tree

        Returns:
            Complete formatted tree as string
        """
        root_path = Path(root_path).resolve()

        # Reset statistics
        self.dir_count = 0
        self.file_count = 0
        self.total_size = 0

        # Build output
        output_lines = []

        # Add title if provided
        if title:
            output_lines.append(f"\n{title}")
            output_lines.append("=" * len(title))

        # Add root directory
        output_lines.append(f"\n{root_path.name}/")

        # Generate tree
        tree_lines = self.generate(root_path)
        output_lines.extend(tree_lines)

        # Add summary
        output_lines.append("")
        summary = f"\n{self.dir_count} directories, {self.file_count} files"
        if self.show_sizes and self.total_size > 0:
            summary += f", {self.format_size(self.total_size)} total"
        output_lines.append(summary)

        return "\n".join(output_lines)


def visualize_directory(
    path: Path,
    title: Optional[str] = None,
    max_depth: Optional[int] = None,
    show_sizes: bool = True,
    show_hidden: bool = False
) -> str:
    """
    Convenience function to visualize a directory tree.

    Args:
        path: Directory path to visualize
        title: Optional title
        max_depth: Maximum depth (None = unlimited)
        show_sizes: Show file sizes
        show_hidden: Show hidden files

    Returns:
        Formatted tree string
    """
    tree = DirectoryTree(
        max_depth=max_depth,
        show_sizes=show_sizes,
        show_hidden=show_hidden
    )
    return tree.visualize(path, title)


def visualize_project_directory(project, max_depth: Optional[int] = 3) -> str:
    """
    Visualize a CCP4 project directory structure.

    Args:
        project: Project model instance
        max_depth: Maximum depth to show

    Returns:
        Formatted tree string
    """
    title = f"Project: {project.name} ({project.uuid})"
    path = Path(project.directory)

    return visualize_directory(
        path,
        title=title,
        max_depth=max_depth,
        show_sizes=True,
        show_hidden=False
    )


def visualize_job_directory(job, max_depth: Optional[int] = 2) -> str:
    """
    Visualize a CCP4 job directory structure.

    Args:
        job: Job model instance
        max_depth: Maximum depth to show

    Returns:
        Formatted tree string
    """
    title = f"Job: {job.number} - {job.title} ({job.task_name})"
    path = Path(job.directory)

    return visualize_directory(
        path,
        title=title,
        max_depth=max_depth,
        show_sizes=True,
        show_hidden=False
    )
