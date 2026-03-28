# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
