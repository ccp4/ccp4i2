"""
Base object system for CData hierarchy.

Provides hierarchical objects, signals, events, and the core CData classes.
"""

# Re-export core classes for easy importing
from .cdata import CData, ValueState
from .ccontainer import CContainer
from .cdata_file import CDataFile
from .cdata_file_content import CDataFileContent

__all__ = [
    'CData',
    'ValueState',
    'CContainer',
    'CDataFile',
    'CDataFileContent',
]
