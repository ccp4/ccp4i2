"""
CCP4Modules.py - Modern module accessor for CCP4 components

This module provides access to CCP4 system components in a modern,
Django-based architecture. It replaces the legacy CCP4Modules class
with simple function accessors.

This module maintains backward compatibility with legacy code that
imports from ccp4i2.core.CCP4Modules.
"""

from .CCP4ProjectsManager import PROJECTSMANAGER
from .CCP4ProcessManager import PROCESSMANAGER


__all__ = ['PROJECTSMANAGER', 'PROCESSMANAGER', 'PREFERENCES']


class _Preferences:
    """Accessor for functional user preferences.

    Legacy/wrapper code reads attributes like ``SHELXDIR``, ``BUSTERDIR``,
    ``COOT_EXECUTABLE``, ``PDB_REDO_TOKEN_ID`` and ``RETAIN_DIAGNOSTIC_FILES`` off
    the preferences object. Each attribute is resolved from the shared store with
    precedence ``env var > ~/.ccp4i2/preferences.json (userPreferences) > None``
    (see ccp4i2.config.preferences). Unknown preferences return None, so attribute
    access never raises — preserving the previous stub's safety.
    """

    def __getattr__(self, name):
        # __getattr__ only fires for attributes not found normally, so this never
        # shadows real methods. Resolve via the shared preferences store.
        from ccp4i2.config import preferences
        return preferences.user_preference(name, default=None)


def PREFERENCES():
    """
    Get the user preferences accessor.

    Replaces the legacy CCP4Modules.PREFERENCES() from classic ccp4i2. Attribute
    access (e.g. ``PREFERENCES().SHELXDIR``) resolves from the shared store —
    ``env var > ~/.ccp4i2/preferences.json (userPreferences) > None`` — so the
    GUI, the i2/i2run CLI, and job execution all see the same values. Unknown
    preferences return None.

    Returns:
        _Preferences: the user preferences accessor

    Example:
        >>> from ccp4i2.core import CCP4Modules
        >>> CCP4Modules.PREFERENCES().SHELXDIR   # value from env/file, else None
    """
    return _Preferences()
