"""
Implementation classes for CCP4RefmacData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType
from ccp4i2.core.base_object.base_classes import CData, CDataFile
from ccp4i2.core.base_object.fundamental_types import CFloat, CInt, CList, CString
from ccp4i2.core.CCP4Data import CUUID
from ccp4i2.core.CCP4File import CFilePath, CProjectId



@cdata_class(
    error_codes={
    },
    qualifiers={
        "fileLabel": 'restraints',
        "mimeTypeName": 'application/refmac-external-restraints',
        "mimeTypeDescription": 'Refmac external restraints',
        "guiLabel": 'Additional restraints',
        "fileExtensions": ['txt'],
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'fileLabel',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag'],
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CRefmacRestraintsDataFile(CDataFile):


    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRestraintsDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CRefmacRestraintsDataFile with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
        "101": {
            "description": "No sequence identity or structure RMS to target set"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
    contents_order=['chain_id', 'residue_1', 'residue_2'],
    content_qualifiers={
        "chain_id": {'charWidth': 1, 'allowUndefined': False, 'mustExist': True},
        "residue_1": {'mustExist': True, 'allowUndefined': False},
        "residue_2": {'mustExist': True, 'allowUndefined': False},
    },
)
class CRefmacRigidGroupSegment(CData):

    chain_id: Optional["CString"] = None
    residue_1: Optional["CInt"] = None
    residue_2: Optional["CInt"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRigidGroupSegment.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CRefmacRigidGroupSegment with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
    content_qualifiers={
        "atomType": {'charWidth': 5, 'toolTip': 'Element name as in PDB file'},
        "Fp": {'toolTip': "Form factor f' for element at given wavelength"},
        "Fpp": {'toolTip': "Form factor f'' for element at given wavelength"},
    },
)
class CRefmacAnomalousAtom(CData):

    atomType: Optional["CString"] = None
    Fp: Optional["CFloat"] = None
    Fpp: Optional["CFloat"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacAnomalousAtom.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CRefmacAnomalousAtom with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "listMinLength": 0,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
)
class CRefmacRigidGroupList(CList):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRigidGroupList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A list with all items of one CData sub-class
    
    Extends CRefmacRigidGroupList with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
    content_qualifiers={
        "segmentList": {'listMinLength': 1},
    },
)
class CRefmacRigidGroupItem(CData):

    rigid_group_id: Optional["CString"] = None
    segmentList: Optional["CList"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRigidGroupItem.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CRefmacRigidGroupItem with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

