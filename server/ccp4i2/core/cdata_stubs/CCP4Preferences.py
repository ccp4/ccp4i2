"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.

This is a stub file - extend classes in core/ to add methods.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any

# Metadata system
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType

# Base classes
from ccp4i2.core.base_object.base_classes import CContainer


@cdata_class(
    error_codes={
        "101": {
            "description": "Error parsing XML"
        },
        "102": {
            "description": "Missing information"
        },
        "103": {
            "description": "Unknown data class"
        },
        "104": {
            "description": "Error creating data object"
        },
        "105": {
            "description": "Error setting data object qualifiers"
        },
        "106": {
            "description": "Error loading container definition"
        },
        "107": {
            "description": "XML file does not have correct function defined in the header"
        },
        "108": {
            "description": "XML undefined error interpreting sub-container"
        },
        "109": {
            "description": "Error attempting to access unknown attribute",
            "severity": 2
        },
        "110": {
            "description": "Error creating sub-container"
        },
        "111": {
            "description": "XML file does not have expected pluginName defined in the header"
        },
        "113": {
            "description": "Attempting to add object that is not a CData"
        },
        "114": {
            "description": "Attempting to add object without valid name"
        },
        "115": {
            "description": "Attempting to add object with name that is already in container"
        },
        "116": {
            "description": "Error while attempting to add object"
        },
        "117": {
            "description": "Attempting to delete object with unrecognised name"
        },
        "118": {
            "description": "Error while attempting to delete object"
        },
        "119": {
            "description": "Error while attempting to set this container as object parent"
        },
        "120": {
            "description": "Attempting to add object of unrecognised class to container contents"
        },
        "121": {
            "description": "Error while attempting to add to container contents"
        },
        "122": {
            "description": "Error while attempting to make object from new content in container"
        },
        "123": {
            "description": "Unknown error while reading container header"
        },
        "124": {
            "description": "Definition of sub-content for data of class that does not require sub-content"
        },
        "125": {
            "description": "Unknown error while reading container content"
        },
        "126": {
            "description": "No id for sub-container in XML file"
        },
        "127": {
            "description": "Attempting to load container data from file that does not exist"
        },
        "128": {
            "description": "Unknown error creating XML for sub-container"
        },
        "129": {
            "description": "Error retieving data object for XML"
        },
        "130": {
            "description": "Error saving data object to XML"
        },
        "131": {
            "description": "Unknown error writing container contents to XML file"
        },
        "132": {
            "description": "Error changing object name - no name given"
        },
        "133": {
            "description": "Error changing object name - object with new name already exists"
        },
        "134": {
            "description": "Error changing object name - no object with old name"
        },
        "135": {
            "description": "Unknown error changing object name"
        },
        "136": {
            "description": "Error inserting object in container data order"
        },
        "137": {
            "description": "Unknown error restoring data from database"
        },
        "138": {
            "description": "Attempting to copy from otherContainer which is not a CContainer"
        },
        "139": {
            "severity": 2,
            "description": "Attempting to copy data which is not in this container"
        },
        "140": {
            "severity": 2,
            "description": "Attempting to copy data which is not in the other container"
        },
        "141": {
            "severity": 2,
            "description": "Unknown error copying data"
        },
        "142": {
            "description": "Unrecognised class name in file"
        },
        "143": {
            "description": "Item in file does not have an id"
        },
        "144": {
            "description": "Item id in file is not unique"
        },
        "145": {
            "description": "Failed setting command line argument"
        },
        "146": {
            "description": "Insufficient arguments at end of command line"
        },
        "147": {
            "description": "Error handling XmlDataFile for file element in def xml"
        },
        "148": {
            "description": "XmlDataFile for file element in def xml: file not found"
        },
        "149": {
            "description": "XmlDataFile for file element in def xml: can not read xml"
        },
        "150": {
            "description": "loadDataFromXml could not find plugin def file"
        },
        "160": {
            "description": "Error in adding guiAdmin to CContainer"
        },
        "161": {
            "description": "Error adding object to guiAdmin"
        },
        "162": {
            "description": "Error adding guiAdmin to CContainer"
        },
        "310": {
            "description": "Different number of file objects to compare"
        },
        "311": {
            "description": "Different number of XData objects to compare"
        },
        "312": {
            "description": "Different number of key-value pairs to compare"
        },
        "313": {
            "description": "Different values of key-value pair"
        },
        "314": {
            "description": "Error running comparison of object"
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool'},
        "default": {'type': 'dict'},
        "toolTip": {'type': 'str'},
        "guiLabel": {'type': 'str'},
        "guiDefinition": {'type': 'dict'},
        "helpFile": {'type': 'str'},
        "saveToDb": {'type': 'bool', 'description': 'Save this data in the database'},
    },
)
class CPreferencesStub(CContainer):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CPreferences.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPreferencesStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
