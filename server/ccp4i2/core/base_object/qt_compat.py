"""Qt API Compatibility Layer

This module provides compatibility wrappers for legacy Qt API patterns
to support the migration from Qt-based CCP4i2 to Qt-free Python.

The main use case is handling legacy code that expects Qt string objects
with a .get() method to extract the underlying Python string.
"""


class QtStringCompat(str):
    """String subclass providing Qt QString/CData API compatibility.

    Legacy Qt code often used QString objects which had a .get() method
    to extract the underlying string value. This class provides that
    compatibility while being a normal Python string in all other respects.

    When created with a parent reference (typically a CDataFile), also
    supports the .set() method to update the parent's path.

    Example:
        # Legacy code pattern:
        path = file_obj.fullPath.get()  # Qt QString had .get() method

        # With QtStringCompat:
        path = file_obj.fullPath.get()  # Still works
        path = str(file_obj.fullPath)   # Also works
        path = file_obj.fullPath        # Also works (it's a string!)

        # Setting path (when parent is provided):
        file_obj.fullPath.set(new_path)  # Delegates to parent.setFullPath()
    """

    def __new__(cls, value, parent=None):
        """Create new QtStringCompat instance.

        Args:
            value: String value
            parent: Parent CDataFile object (optional, required for .set())
        """
        instance = super().__new__(cls, value)
        return instance

    def __init__(self, value, parent=None):
        """Initialize with optional parent reference.

        Args:
            value: String value (handled by __new__)
            parent: Parent CDataFile object for .set() delegation
        """
        # Note: str.__init__ doesn't take arguments, value is set in __new__
        self._parent = parent

    def get(self):
        """Return the string value (Qt QString compatibility).

        Returns:
            The string value (self)
        """
        return str(self)

    def set(self, value):
        """Set the value on the parent CDataFile.

        This delegates to the parent's setFullPath() method, which handles:
        - Parsing the path into baseName, relPath, project
        - Checking if the file exists
        - Updating all related attributes

        Args:
            value: New path value (string or Path-like object)

        Raises:
            RuntimeError: If no parent was provided when this object was created
        """
        if self._parent is None:
            raise RuntimeError(
                "Cannot call .set() on fullPath without parent reference. "
                "Use file_obj.setFullPath(path) or file_obj.fullPath = path instead."
            )
        self._parent.setFullPath(str(value))

    def isSet(self):
        """Check if the string has a value (CData compatibility).

        Returns:
            bool: True if the string is non-empty, False otherwise
        """
        return bool(self)

    def __repr__(self):
        """Return repr showing this is a compatibility wrapper."""
        return f"QtStringCompat({str.__repr__(self)})"
