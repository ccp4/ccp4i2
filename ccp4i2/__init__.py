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

    # Check if running in Django mode
    from ccp4i2 import is_django_mode
    if is_django_mode():
        # Django-specific code
        pass
"""

import os
from pathlib import Path

__version__ = "3.0.0-dev"

# Determine the root directory of the ccp4i2 installation
# Priority: CCP4I2_ROOT env var > package location
_ROOT_DIR = None


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
    global _ROOT_DIR

    if _ROOT_DIR is not None:
        return _ROOT_DIR

    # Check environment variable first
    env_root = os.environ.get('CCP4I2_ROOT')
    if env_root and os.path.isdir(env_root):
        _ROOT_DIR = Path(env_root)
        return _ROOT_DIR

    # Fall back to package location
    # This file is at ccp4i2/__init__.py, so parent.parent is the root
    _ROOT_DIR = Path(__file__).parent.parent
    return _ROOT_DIR


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


def is_django_mode() -> bool:
    """
    Check if running in Django backend mode (vs Qt mode).

    Returns:
        True if running in Django mode, False otherwise
    """
    # Check explicit environment variable
    backend = os.environ.get('CCP4I2_BACKEND', '').lower()
    if backend == 'django':
        return True
    if backend == 'qt':
        return False

    # Check for Django settings module
    if os.environ.get('DJANGO_SETTINGS_MODULE'):
        return True

    # Try import detection
    try:
        from PySide2 import QtCore
        return False  # Qt is available
    except ImportError:
        return True  # Assume Django mode if Qt not available


def is_qt_mode() -> bool:
    """
    Check if running in Qt mode (vs Django mode).

    Returns:
        True if running in Qt mode, False otherwise
    """
    return not is_django_mode()
