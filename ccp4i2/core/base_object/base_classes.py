"""
Base classes for the CData hierarchy - Re-export module for backward compatibility.

This module re-exports all base classes from their individual files.
Existing code can continue to import from base_classes.py without changes.

New code should import directly from the specific modules:
    from ccp4i2.core.base_object.cdata import CData, ValueState
    from ccp4i2.core.base_object.cdata_file_content import CDataFileContent
    from ccp4i2.core.base_object.cdata_file import CDataFile
    from ccp4i2.core.base_object.ccontainer import CContainer
"""

# Re-export ValueState and CData
from .cdata import ValueState, CData

# Re-export CDataFileContent
from .cdata_file_content import CDataFileContent

# Re-export CDataFile
from .cdata_file import CDataFile

# Re-export CContainer
from .ccontainer import CContainer

# Expose everything at module level for backward compatibility
__all__ = [
    'ValueState',
    'CData',
    'CDataFileContent',
    'CDataFile',
    'CContainer',
]
