"""
CCP4RefmacData.py --- CData classes.

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
from ccp4i2.core.base_object.base_classes import CData, CDataFile
from ccp4i2.core.base_object.fundamental_types import CFloat, CInt, CList, CString
from ccp4i2.core.CCP4Data import CUUID
from ccp4i2.core.CCP4File import CFilePath, CProjectId



class CRefmacRestraintsDataFile(CDataFile):



    class Meta:
        qualifiers = {
            "fileLabel": 'restraints',
            "mimeTypeName": 'application/refmac-external-restraints',
            "mimeTypeDescription": 'Refmac external restraints',
            "guiLabel": 'Additional restraints',
            "fileExtensions": ['txt'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
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


class CRefmacRigidGroupSegment(CData):


    class Meta:
        error_codes = {
            "101": {
                "description": "No sequence identity or structure RMS to target set"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    chain_id = content("CString", charWidth=1, allowUndefined=False, mustExist=True)
    residue_1 = content("CInt", mustExist=True, allowUndefined=False)
    residue_2 = content("CInt", mustExist=True, allowUndefined=False)

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


class CRefmacAnomalousAtom(CData):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    atomType = content("CString", charWidth=5, toolTip='Element name as in PDB file')
    Fp = content("CFloat", toolTip="Form factor f' for element at given wavelength")
    Fpp = content("CFloat", toolTip="Form factor f'' for element at given wavelength")

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


class CRefmacRigidGroupList(CList):


    class Meta:
        qualifiers = {
            "listMinLength": 0,
        }
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


class CRefmacRigidGroupItem(CData):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    rigid_group_id = content("CString")
    segmentList = content("CList", listMinLength=1)

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

