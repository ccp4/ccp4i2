"""
Implementation classes for CCP4Data.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any

from ccp4i2.core.cdata_stubs.CCP4Data import CBaseDataStub, CCollectionStub, CDictStub, CFloatRangeStub, CFollowFromJobStub, CI2DataTypeStub, CIntRangeStub, CJobStatusStub, CJobTitleStub, COneWordStub, COutputFileListStub, CPatchSelectionStub, CRangeStub, CRangeSelectionStub, CUUIDStub

# Re-export fundamental types for legacy code compatibility
# Many legacy files use "CCP4Data.CList", "CCP4Data.CString", etc.
# which are actually in base_object.fundamental_types
from ccp4i2.core.base_object.fundamental_types import CList, CString, CInt, CFloat, CBoolean


class CBaseData(CBaseDataStub):
    """
    Base class for simple classes
    
    Extends CBaseDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CCollection(CCollectionStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CCollectionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDict(CDictStub, CCollection):
    """
    Dictionary-like CData type for key-value storage.

    Inherits from:
    - CDictStub: Metadata and structure
    - CCollection: Shared full-fat methods

    Used for selection dictionaries in ASU files and similar key-value data.
    Behaves like a Python dict with __getitem__, __setitem__, etc.

    For unset keys, returns True by default (all items selected).
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize CDict with an empty internal dictionary."""
        super().__init__(parent=parent, name=name, **kwargs)
        self._dict_data: dict = {}

    def __getitem__(self, key: str) -> Any:
        """Get value for key. Returns True if key not found (default selected)."""
        return self._dict_data.get(key, True)

    def __setitem__(self, key: str, value: Any) -> None:
        """Set value for key."""
        self._dict_data[key] = value

    def __delitem__(self, key: str) -> None:
        """Delete key from dictionary."""
        del self._dict_data[key]

    def __contains__(self, key: str) -> bool:
        """Check if key exists in dictionary."""
        return key in self._dict_data

    def __len__(self) -> int:
        """Return number of keys in dictionary."""
        return len(self._dict_data)

    def __iter__(self):
        """Iterate over keys."""
        return iter(self._dict_data)

    def keys(self):
        """Return dictionary keys."""
        return self._dict_data.keys()

    def values(self):
        """Return dictionary values."""
        return self._dict_data.values()

    def items(self):
        """Return dictionary items."""
        return self._dict_data.items()

    def get(self, key: str = None, default: Any = True) -> Any:
        """Get value for key with default, or return all data if no key specified.

        When called with no arguments (for CData serialization), returns the
        internal dictionary data. When called with a key argument (dict-style
        access), returns the value for that key with the given default.

        Args:
            key: Optional key to look up. If None, returns all data.
            default: Default value if key not found (default: True for "selected")

        Returns:
            Value for key, or dict of all data if no key specified.
        """
        if key is None:
            # CData serialization - return all data
            return self._dict_data.copy()
        return self._dict_data.get(key, default)

    def setdefault(self, key: str, default: Any = True) -> Any:
        """Set default value for key if not present."""
        return self._dict_data.setdefault(key, default)

    def update(self, other: dict = None, **kwargs) -> None:
        """Update dictionary with another dict or keyword args."""
        if other:
            self._dict_data.update(other)
        self._dict_data.update(kwargs)

    def clear(self) -> None:
        """Clear all items from dictionary."""
        self._dict_data.clear()


class CFloatRange(CFloatRangeStub):
    """
    Two floats defining start and end of range

    Extends CFloatRangeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize CFloatRange with .start and .end not set by default.

        This ensures that range fields are only marked as EXPLICITLY_SET when
        actually assigned by user code, preventing them from appearing in
        serialized XML when using excludeUnset=True.

        Note: Smart assignment is now handled by the base CData class, which
        only copies CData fields that are explicitly set (isSet(allowDefault=False)).
        """
        super().__init__(parent=parent, name=name, **kwargs)

        # Mark .start and .end as NOT_SET so they won't be serialized unless explicitly set
        from ccp4i2.core.base_object.cdata import ValueState
        if hasattr(self, 'start') and hasattr(self.start, '_value_states'):
            self.start._value_states['value'] = ValueState.NOT_SET
        if hasattr(self, 'end') and hasattr(self.end, '_value_states'):
            self.end._value_states['value'] = ValueState.NOT_SET


class CFollowFromJob(CFollowFromJobStub):
    """
    A string
    
    Extends CFollowFromJobStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CI2DataType(CI2DataTypeStub):
    """
    A string
    
    Extends CI2DataTypeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CIntRange(CIntRangeStub):
    """
    Two integers defining start and end of range

    Extends CIntRangeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize CIntRange with .start and .end not set by default.

        This ensures that range fields are only marked as EXPLICITLY_SET when
        actually assigned by user code, preventing them from appearing in
        serialized XML when using excludeUnset=True.

        Note: Smart assignment is handled by the base CData class, which
        only copies CData fields that are explicitly set (isSet(allowDefault=False)).
        """
        super().__init__(parent=parent, name=name, **kwargs)

        # Mark .start and .end as NOT_SET so they won't be serialized unless explicitly set
        from ccp4i2.core.base_object.cdata import ValueState
        if hasattr(self, 'start') and hasattr(self.start, '_value_states'):
            self.start._value_states['value'] = ValueState.NOT_SET
        if hasattr(self, 'end') and hasattr(self.end, '_value_states'):
            self.end._value_states['value'] = ValueState.NOT_SET


class CJobStatus(CJobStatusStub):
    """
    An integer
    
    Extends CJobStatusStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CJobTitle(CJobTitleStub):
    """
    A string
    
    Extends CJobTitleStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class COneWord(COneWordStub):
    """
    A single word string - no white space
    
    Extends COneWordStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class COutputFileList(COutputFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends COutputFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPatchSelection(CPatchSelectionStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CPatchSelectionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRange(CRangeStub):
    """
    Base class for CIntRange and CFloatRange
    
    Extends CRangeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRangeSelection(CRangeSelectionStub):
    """
    Range selection string (e.g., "1-10,15,20-25").

    Extends CRangeSelectionStub with validation for range selection syntax.
    Validates format like "1-10,15,20-25" where ranges are specified as
    comma-separated numbers or hyphen-separated ranges.
    """

    ERROR_CODES = {
        201: {'description': 'Range selection contains invalid character'},
        202: {'description': 'Range selection contains bad syntax'}
    }

    def validity(self, arg):
        """
        Validate range selection string syntax.

        Legacy API compatibility method for validating range selection format.
        Checks for:
        - Only digits, commas, and hyphens
        - Valid range syntax (start-end where start < end)
        - No empty ranges

        Args:
            arg: String to validate (e.g., "1-10,15,20-25")

        Returns:
            CErrorReport: Error report with validation issues
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from ccp4i2.core.base_object.base_classes import CData

        err = CErrorReport()

        # Check for undefined value
        if arg is None:
            if not self.get_qualifier('allowUndefined'):
                err.append(
                    "CData", 2,
                    details="Value is not set",
                    name=self.object_path()
                )
            return err

        # Remove whitespace
        arg = self.removeWhiteSpace(arg)

        # Check for invalid characters (only allow digits, comma, hyphen)
        import re
        s = re.search(r'[^0-9,\-]', arg)
        if s is not None:
            err.append(
                "CRangeSelection", 201,
                details=self.ERROR_CODES[201]['description'],
                name=self.object_path()
            )
        else:
            # Check each range component
            rList = arg.split(',')
            for r in rList:
                if len(r) < 1:
                    # Empty range (e.g., "1,,2")
                    err.append(
                        "CRangeSelection", 202,
                        details=self.ERROR_CODES[202]['description'],
                        name=self.object_path()
                    )
                elif r.count('-') > 1:
                    # Too many hyphens (e.g., "1-2-3")
                    err.append(
                        "CRangeSelection", 202,
                        details=self.ERROR_CODES[202]['description'],
                        name=self.object_path()
                    )
                elif r.count('-') == 1:
                    # Range format "start-end"
                    rr = r.split('-')
                    try:
                        if int(rr[0]) > int(rr[1]):
                            # Start is greater than end (e.g., "10-5")
                            err.append(
                                "CRangeSelection", 202,
                                details=self.ERROR_CODES[202]['description'],
                                name=self.object_path()
                            )
                    except:
                        # Invalid integers
                        err.append(
                            "CRangeSelection", 202,
                            details=self.ERROR_CODES[202]['description'],
                            name=self.object_path()
                        )

        return err


class CUUID(CUUIDStub):
    """
    A string
    
    Extends CUUIDStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

