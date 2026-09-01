"""
CCP4CustomTaskManager.py --- CData classes.

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
from ccp4i2.core.base_object.base_classes import CContainer, CData
from ccp4i2.core.base_object.fundamental_types import CBoolean, CList, CString
from ccp4i2.core.CCP4Data import CI2DataType, COneWord



@cdata_class(
    error_codes={
    },
    qualifiers={
        "enumerators": ['unknown', 'input', 'output', 'control parameter', 'log'],
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
)
class CCustomTaskFileFunction(CString):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomTaskFileFunction.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A string
    
    Extends CCustomTaskFileFunction with implementation-specific methods.
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
class CCustomTaskParamList(CList):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomTaskParamList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A list with all items of one CData sub-class
    
    Extends CCustomTaskParamList with implementation-specific methods.
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
)
class CCustomComFile(CData):

    text = content("CString")
    name = content("CString", default='./com.txt')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomComFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CCustomComFile with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


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
)
class CCustomTaskDefinition(CContainer):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomTaskDefinition.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CCustomTaskDefinition with implementation-specific methods.
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
class CCustomComFileList(CList):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomComFileList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A list with all items of one CData sub-class
    
    Extends CCustomComFileList with implementation-specific methods.
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
)
class CCustomTaskParam(CData):

    name = content("COneWord")
    dataType = content("CI2DataType", default='CPdbDataFile')
    label = content("CString")
    obligatory = content("CBoolean", default=True)
    saveDataToDb = content("CBoolean", default=False)
    function = content("CCustomTaskFileFunction", default='input')
    mergeTo = content("CString")
    splitColumns = content("CString")
    requiredContentType = content("CList")
    outputFilePath = content("CString")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCustomTaskParam.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CCustomTaskParam with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

