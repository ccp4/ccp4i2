"""
CCP4i2 - Graphical interface and scripting environment for CCP4

This package provides the core functionality for CCP4i2, including:
- Task wrappers for crystallographic programs
- Pipelines that chain multiple programs
- Data management utilities
- Report generation

Usage:
    # Get the root directory for resource files (icons, data, etc.)
    from ccp4i2 import get_resource_path
    icons_dir = get_resource_path('qticons')
"""

from os import environ
from datetime import datetime
from pathlib import Path

MAJOR = 3
MINOR = 0
PATCH = 0

__version__ = f"{MAJOR}.{MINOR}.{PATCH}"
__version_date__ = datetime(2025, 12, 10)
__version_info__ = (MAJOR, MINOR, PATCH)


def get_root_dir() -> Path:
    """
    Get the root directory of the ccp4i2 installation.

    This is used to locate non-Python resources like icons and data files.

    Priority:
    1. CCP4I2_ROOT environment variable (for development/editable installs)
    2. Parent of this package's location (for regular installs)

    Returns:
        Path to the ccp4i2 root directory
    """
    # Check environment variable first
    root = environ.get('CCP4I2_ROOT')
    if root:
        path = Path(root)
        if path.is_dir():
            return path

    # Fall back to package location
    return Path(__file__).parent


def get_resource_path(resource_name: str) -> Path:
    """
    Get the path to a resource directory or file.

    Args:
        resource_name: Name of the resource (e.g., 'qticons', 'svgicons', 'data')

    Returns:
        Path to the resource

    Example:
        >>> from ccp4i2 import get_resource_path
        >>> icons = get_resource_path('qticons')
        >>> icon_file = icons / 'refmac.png'
    """
    root = get_root_dir()
    resource_path = root / resource_name

    if not resource_path.exists():
        # Try looking in share/ccp4i2 (for system installs)
        share_path = root / 'share' / 'ccp4i2' / resource_name
        if share_path.exists():
            return share_path

    return resource_path
