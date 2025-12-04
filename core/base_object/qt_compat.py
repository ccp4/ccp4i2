"""Qt API Compatibility Layer

This module provides compatibility wrappers for legacy Qt API patterns
to support the migration from Qt-based CCP4i2 to Qt-free Python.

The main use case is handling legacy code that expects Qt string objects
with a .get() method to extract the underlying Python string.
"""


class QtStringCompat(str):
    """String subclass providing Qt QString API compatibility.

    Legacy Qt code often used QString objects which had a .get() method
    to extract the underlying string value. This class provides that
    compatibility while being a normal Python string in all other respects.

    Example:
        # Legacy code pattern:
        path = file_obj.fullPath.get()  # Qt QString had .get() method

        # With QtStringCompat:
        path = file_obj.fullPath.get()  # Still works
        path = str(file_obj.fullPath)   # Also works
        path = file_obj.fullPath        # Also works (it's a string!)
    """

    def get(self):
        """Return the string value (Qt QString compatibility).

        Returns:
            The string value (self)
        """
        return str(self)

    def isSet(self):
        """Check if the string has a value (CData compatibility).

        Returns:
            bool: True if the string is non-empty, False otherwise
        """
        return bool(self)

    def __repr__(self):
        """Return repr showing this is a compatibility wrapper."""
        return f"QtStringCompat({str.__repr__(self)})"
