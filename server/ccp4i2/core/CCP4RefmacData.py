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



    """External restraints supplied to Refmac."""

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

    # Add your methods here
    pass


class CRefmacRigidGroupSegment(CData):


    """One chain segment of a Refmac rigid body group."""

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
    chain_id = content(
        "CString",
        charWidth=1,
        allowUndefined=False,
        mustExist=True,
        guiLabel='Chain',
        toolTip='Chain identifier of the segment')
    residue_1 = content(
        "CInt",
        mustExist=True,
        allowUndefined=False,
        guiLabel='First residue',
        toolTip='Number of the first residue in the segment')
    residue_2 = content(
        "CInt",
        mustExist=True,
        allowUndefined=False,
        guiLabel='Last residue',
        toolTip='Number of the last residue in the segment')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRigidGroupSegment.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CRefmacAnomalousAtom(CData):


    """
    An anomalous scatterer and its scattering factors at the wavelength
    in use.
    """

    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    atomType = content(
        "CString",
        charWidth=5,
        toolTip='Element name as in PDB file',
        guiLabel='Element')
    Fp = content(
        "CFloat",
        toolTip="Form factor f' for element at given wavelength",
        guiLabel="f'")
    Fpp = content(
        "CFloat",
        toolTip="Form factor f'' for element at given wavelength",
        guiLabel="f''")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacAnomalousAtom.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CRefmacRigidGroupList(CList):


    """A list with all items of one CData sub-class"""
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

    # Add your methods here
    pass


class CRefmacRigidGroupItem(CData):


    """One Refmac rigid body group: the chain segments that move together."""

    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    rigid_group_id = content(
        "CString",
        guiLabel='Group',
        toolTip='Identifier for this rigid body group')
    segmentList = content(
        "CList",
        listMinLength=1,
        guiLabel='Segments',
        toolTip='The chain segments that move together as one rigid body')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefmacRigidGroupItem.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass

