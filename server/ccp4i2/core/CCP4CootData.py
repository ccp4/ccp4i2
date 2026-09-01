"""
CCP4CootData.py --- CData classes.

These were once two classes each: a generated stub carrying the data
model, and an implementation carrying the methods. The generator is no
longer run and the split cost more than it saved --- the two halves
interleaved in the MRO, so an implementation could drop out of its own
subclass's ancestry and `isinstance` would say no to a file that plainly
was one. They are one class now.
"""

from typing import Optional, Any
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType, content
from ccp4i2.core.base_object.base_classes import CDataFile
from ccp4i2.core.base_object.fundamental_types import CInt, CString
from ccp4i2.core.CCP4Data import CUUID
from ccp4i2.core.CCP4File import CFilePath, CProjectId



class CCootHistoryDataFile(CDataFile):

    # Subtype constants

    """
    Extends CCootHistoryDataFile with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "mimeTypeName": 'application/coot-script',
            "mimeTypeDescription": 'Coot history/script file',
            "fileExtensions": ['scm', 'py'],
            "fileContentClassName": None,
            "guiLabel": 'Coot history',
            "fileLabel": 'coot_history',
            "toolTip": 'history.scm or 0-state.scm file from Coot',
        }
        content_qualifiers = {
            "subType": {'default': 2, 'enumerators': [1, 2], 'onlyEnumerators': True, 'menuText': ['Coot 0-state.scm', 'Coot history.scm']},
            "contentFlag": {'min': 0, 'default': None},
        }
    SUBTYPE_INITIAL = 1  # Coot 0-state.scm
    SUBTYPE_HISTORY = 2  # Coot history.scm


    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCootHistoryDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass

