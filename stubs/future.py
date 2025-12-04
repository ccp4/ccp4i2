"""
Minimal stub for future module.

The 'future' module provides Python 2/3 compatibility utilities.
This stub provides just enough to allow imports to succeed.

Most functionality will raise NotImplementedError if actually called.
This stub allows plugins to import future for registration purposes.
"""

__version__ = "stub-0.0.1"


class standard_library:
    """Stub for future.standard_library."""

    @staticmethod
    def install_aliases():
        """Install Python 3 style imports (stub - does nothing)."""
        pass


class builtins:
    """Stub for future.builtins."""
    pass


class utils:
    """Stub for future.utils."""

    @staticmethod
    def with_metaclass(meta, *bases):
        """Create a base class with a metaclass (stub)."""
        class metaclass(type):
            def __new__(cls, name, this_bases, d):
                return meta(name, bases, d)
        return type.__new__(metaclass, 'temporary_class', (), {})
