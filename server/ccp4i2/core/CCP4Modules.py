"""
CCP4Modules.py - Modern module accessor for CCP4 components

This module provides access to CCP4 system components in a modern,
Django-based architecture. It replaces the legacy CCP4Modules class
with simple function accessors.

This module maintains backward compatibility with legacy code that
imports from ccp4i2.core.CCP4Modules.
"""

from .CCP4TaskManager import TASKMANAGER
from .CCP4ProjectsManager import PROJECTSMANAGER
from .CCP4ProcessManager import PROCESSMANAGER


__all__ = ['TASKMANAGER', 'PROJECTSMANAGER', 'PROCESSMANAGER', 'PREFERENCES']


class _PreferencesStub:
    """Stub for legacy PREFERENCES object.

    Legacy code checks for attributes like SHELXDIR on the preferences object.
    This stub allows attribute access without raising errors.
    """
    def __init__(self):
        pass

    def __getattr__(self, name):
        # Return None for any requested attribute (preferences not set)
        return None


def PREFERENCES():
    """
    Get the user preferences object (stub for legacy compatibility).

    This replaces the legacy CCP4Modules.PREFERENCES() function from classic ccp4i2.
    Returns a stub object that returns None for any attribute access.

    Returns:
        _PreferencesStub: A stub preferences object

    Example:
        >>> from ccp4i2.core import CCP4Modules
        >>> prefs = CCP4Modules.PREFERENCES()
        >>> shelx_dir = prefs.SHELXDIR  # Returns None
    """
    return _PreferencesStub()
