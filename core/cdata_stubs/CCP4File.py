"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.

This is a stub file - extend classes in core/ to add methods.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any

# Metadata system
from core.base_object.class_metadata import cdata_class, attribute, AttributeType

# Base classes
from core.base_object.base_classes import CData, CDataFile, CDataFileContent

# Fundamental types
from core.base_object.fundamental_types import CInt, CList, CString

# Cross-file stub class references
from core.cdata_stubs.CCP4Annotation import CHostNameStub, CTimeStub, CUserIdStub
from core.cdata_stubs.CCP4Data import CUUIDStub


@cdata_class(
    error_codes={
        "201": {
            "description": "Unrecognised projectId"
        },
        "202": {
            "description": "Project does not have directory set"
        },
        "203": {
            "description": "Project directory does not exist"
        },
        "205": {
            "severity": 2,
            "description": "Warning - Project does not have directory set"
        },
        "206": {
            "severity": 2,
            "description": "Warning - Project directory does not exist"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "allowUnfound": True,
        "default": None,
    },
    qualifiers_order=['allowUndefined', 'allowUnfound', 'default'],
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if allow undefined value at run time'},
        "allowUnfound": {'type': 'bool', 'description': 'Flag if allow unfound project at run time'},
        "default": {'type': 'str'},
    },
)
class CProjectIdStub(CUUIDStub):
    """
    The CCP4i2 database project id - a global unique id

    This is a pure data class stub. Extend it in core/CProjectId.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProjectIdStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Version is not of form n.m or n.m.i"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "default": None,
        "charWidth": 10,
    },
    qualifiers_order=['allowUndefined', 'default', 'charWidth'],
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if allow an unset value at run time'},
        "default": {'description': 'A default value'},
        "charWidth": {'type': 'int', 'description': 'Number of characters allowed for widget in GUI'},
    },
)
class CVersionStub(CString):
    """
    A (string) version number of the form n.m.i

    This is a pure data class stub. Extend it in core/CVersion.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CVersionStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "0": {
            "severity": 0,
            "description": "OK"
        },
        "1": {
            "severity": 1,
            "description": "Data has undefined value"
        },
        "2": {
            "severity": 3,
            "description": "Data has undefined value"
        },
        "3": {
            "severity": 2,
            "description": "Missing data"
        },
        "4": {
            "description": "Missing data"
        },
        "5": {
            "description": "Attempting to set data of wrong type"
        },
        "6": {
            "description": "Default value does not satisfy validity check"
        },
        "7": {
            "severity": 2,
            "description": "Unrecognised qualifier in data input"
        },
        "8": {
            "severity": 2,
            "description": "Attempting to get inaccessible attribute:"
        },
        "9": {
            "description": "Failed to get property"
        },
        "10": {
            "severity": 2,
            "description": "Attempting to set inaccessible attribute:"
        },
        "11": {
            "description": "Failed to set property:"
        },
        "12": {
            "description": "Undetermined error setting value from XML"
        },
        "13": {
            "description": "Unrecognised class name in qualifier"
        },
        "14": {
            "severity": 2,
            "description": "No object name when saving qualifiers to XML"
        },
        "15": {
            "description": "Error saving qualifier to XML"
        },
        "16": {
            "severity": 2,
            "description": "Unrecognised item in XML data file"
        },
        "17": {
            "description": "Attempting to set unrecognised qualifier"
        },
        "18": {
            "description": "Attempting to set qualifier with wrong type"
        },
        "19": {
            "description": "Attempting to set qualifier with wrong list item type"
        },
        "20": {
            "description": "Error creating a list/dict item object"
        },
        "21": {
            "description": "Unknown error setting qualifiers from Xml file"
        },
        "22": {
            "description": "Unknown error testing validity"
        },
        "23": {
            "description": "Error saving data object to XML"
        },
        "24": {
            "description": "Unable to test validity of default",
            "severity": 2
        },
        "300": {
            "description": "Compared objects are the same",
            "severity": 0
        },
        "315": {
            "description": "Both compared objects are null",
            "severity": 0
        },
        "301": {
            "description": "Unable to compare this class of data",
            "severity": 2
        },
        "302": {
            "description": "Other data has null value"
        },
        "303": {
            "description": "My data has null value"
        },
        "304": {
            "description": "Data has different values"
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
class CMmcifDataStub(CDataFileContent):
    """
    Generic mmCIF data.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.

    This is a pure data class stub. Extend it in core/CMmcifData.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMmcifDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "exeName": attribute(AttributeType.STRING),
        "exePath": attribute(AttributeType.CUSTOM, custom_class="CDataFile"),
    },
    error_codes={
        "0": {
            "severity": 0,
            "description": "OK"
        },
        "1": {
            "severity": 1,
            "description": "Data has undefined value"
        },
        "2": {
            "severity": 3,
            "description": "Data has undefined value"
        },
        "3": {
            "severity": 2,
            "description": "Missing data"
        },
        "4": {
            "description": "Missing data"
        },
        "5": {
            "description": "Attempting to set data of wrong type"
        },
        "6": {
            "description": "Default value does not satisfy validity check"
        },
        "7": {
            "severity": 2,
            "description": "Unrecognised qualifier in data input"
        },
        "8": {
            "severity": 2,
            "description": "Attempting to get inaccessible attribute:"
        },
        "9": {
            "description": "Failed to get property"
        },
        "10": {
            "severity": 2,
            "description": "Attempting to set inaccessible attribute:"
        },
        "11": {
            "description": "Failed to set property:"
        },
        "12": {
            "description": "Undetermined error setting value from XML"
        },
        "13": {
            "description": "Unrecognised class name in qualifier"
        },
        "14": {
            "severity": 2,
            "description": "No object name when saving qualifiers to XML"
        },
        "15": {
            "description": "Error saving qualifier to XML"
        },
        "16": {
            "severity": 2,
            "description": "Unrecognised item in XML data file"
        },
        "17": {
            "description": "Attempting to set unrecognised qualifier"
        },
        "18": {
            "description": "Attempting to set qualifier with wrong type"
        },
        "19": {
            "description": "Attempting to set qualifier with wrong list item type"
        },
        "20": {
            "description": "Error creating a list/dict item object"
        },
        "21": {
            "description": "Unknown error setting qualifiers from Xml file"
        },
        "22": {
            "description": "Unknown error testing validity"
        },
        "23": {
            "description": "Error saving data object to XML"
        },
        "24": {
            "description": "Unable to test validity of default",
            "severity": 2
        },
        "300": {
            "description": "Compared objects are the same",
            "severity": 0
        },
        "315": {
            "description": "Both compared objects are null",
            "severity": 0
        },
        "301": {
            "description": "Unable to compare this class of data",
            "severity": 2
        },
        "302": {
            "description": "Other data has null value"
        },
        "303": {
            "description": "My data has null value"
        },
        "304": {
            "description": "Data has different values"
        }
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
    contents_order=['exeName', 'exePath'],
    content_qualifiers={
        "exePath": {'mustExist': True, 'allowUndefined': False},
    },
)
class CExePathStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CExePath.py
    to add methods and implementation-specific functionality.
    """

    exeName: Optional[CString] = None
    exePath: Optional[CDataFile] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExePathStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Invalid project name"
        },
        "102": {
            "description": "Project does not have directory set"
        },
        "103": {
            "description": "Project directory does not exist"
        },
        "104": {
            "severity": 2,
            "description": "Warning - Project name is a directory alias"
        },
        "105": {
            "severity": 2,
            "description": "Warning - Project does not have directory set"
        },
        "106": {
            "severity": 2,
            "description": "Warning - Project directory does not exist"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "allowAlias": True,
        "allowUnfound": True,
        "default": None,
    },
    qualifiers_order=[
        'allowUndefined',
        'allowAlias',
        'allowUnfound',
        'default'],
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if allow undefined value at run time'},
        "allowAlias": {'type': 'bool', 'description': 'Flag if allow project to be directory alias at run time'},
        "allowUnfound": {'type': 'bool', 'description': 'Flag if allow unfound project at run time'},
        "default": {'type': 'str'},
    },
)
class CProjectNameStub(CString):
    """
    The name of a CCP4i project or directory alias

    This is a pure data class stub. Extend it in core/CProjectName.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProjectNameStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "String too short"
        },
        "102": {
            "description": "String too long"
        },
        "103": {
            "description": "not one of limited allowed values"
        },
        "104": {
            "description": "Contains disallowed characters"
        }
    },
    qualifiers={
        "enumerators": ['DEF', 'PARAMS', 'LOG', 'PROJECTDIRECTORIES', 'COM', 'REFMAC', 'OUTPUT', 'STATUS', 'PROJECTDATABASE', 'MGSCENE', 'JOBSERVERSTATUS', 'WORKFLOW', 'COMFILEPATCH', 'CUSTOMTASK', 'IMPORTEDJOB', 'I1SUPPLEMENT', 'ASUCONTENT', 'UNKNOWN'],
        "onlyEnumerators": True,
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "enumerators": {'type': 'list'},
        "onlyEnumerators": {'type': 'bool', 'editable': False},
    },
)
class CFileFunctionStub(CString):
    """
    List of recognised XML file functions

    This is a pure data class stub. Extend it in core/CFileFunction.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFileFunctionStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "exportId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
    },
    error_codes={
        "0": {
            "severity": 0,
            "description": "OK"
        },
        "1": {
            "severity": 1,
            "description": "Data has undefined value"
        },
        "2": {
            "severity": 3,
            "description": "Data has undefined value"
        },
        "3": {
            "severity": 2,
            "description": "Missing data"
        },
        "4": {
            "description": "Missing data"
        },
        "5": {
            "description": "Attempting to set data of wrong type"
        },
        "6": {
            "description": "Default value does not satisfy validity check"
        },
        "7": {
            "severity": 2,
            "description": "Unrecognised qualifier in data input"
        },
        "8": {
            "severity": 2,
            "description": "Attempting to get inaccessible attribute:"
        },
        "9": {
            "description": "Failed to get property"
        },
        "10": {
            "severity": 2,
            "description": "Attempting to set inaccessible attribute:"
        },
        "11": {
            "description": "Failed to set property:"
        },
        "12": {
            "description": "Undetermined error setting value from XML"
        },
        "13": {
            "description": "Unrecognised class name in qualifier"
        },
        "14": {
            "severity": 2,
            "description": "No object name when saving qualifiers to XML"
        },
        "15": {
            "description": "Error saving qualifier to XML"
        },
        "16": {
            "severity": 2,
            "description": "Unrecognised item in XML data file"
        },
        "17": {
            "description": "Attempting to set unrecognised qualifier"
        },
        "18": {
            "description": "Attempting to set qualifier with wrong type"
        },
        "19": {
            "description": "Attempting to set qualifier with wrong list item type"
        },
        "20": {
            "description": "Error creating a list/dict item object"
        },
        "21": {
            "description": "Unknown error setting qualifiers from Xml file"
        },
        "22": {
            "description": "Unknown error testing validity"
        },
        "23": {
            "description": "Error saving data object to XML"
        },
        "24": {
            "description": "Unable to test validity of default",
            "severity": 2
        },
        "300": {
            "description": "Compared objects are the same",
            "severity": 0
        },
        "315": {
            "description": "Both compared objects are null",
            "severity": 0
        },
        "301": {
            "description": "Unable to compare this class of data",
            "severity": 2
        },
        "302": {
            "description": "Other data has null value"
        },
        "303": {
            "description": "My data has null value"
        },
        "304": {
            "description": "Data has different values"
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
class CExportedFileStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CExportedFile.py
    to add methods and implementation-specific functionality.
    """

    exportId: Optional[CUUIDStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExportedFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "List shorter than required minimum length"
        },
        "102": {
            "description": "List longer than required maximum length"
        },
        "103": {
            "description": "Consecutive values in list fail comparison test"
        },
        "104": {
            "description": "Attempting to add object of wrong type"
        },
        "105": {
            "description": "Attempting to add object of correct type but wrong qualifiers"
        },
        "106": {
            "description": "Attempting to add data which does not satisfy the qualifiers for a list item"
        },
        "107": {
            "description": "Deleting item will reduce list below minimum length"
        },
        "108": {
            "description": "Adding item will extend list beyond maximum length"
        },
        "109": {
            "description": "Invalid item class"
        },
        "110": {
            "description": "etree (XML) list item of wrong type"
        },
        "112": {
            "description": "No list item object set for list"
        }
    },
    qualifiers={
        "listMinLength": 0,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
    qualifiers_definition={
        "default": {'type': 'list'},
        "listMaxLength": {'type': 'int', 'description': 'Inclusive maximum length of list'},
        "listMinLength": {'type': 'int', 'description': 'Inclusive minimum length of list'},
        "listCompare": {'type': 'int', 'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
    },
)
class CExportedFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CExportedFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExportedFileListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "List shorter than required minimum length"
        },
        "102": {
            "description": "List longer than required maximum length"
        },
        "103": {
            "description": "Consecutive values in list fail comparison test"
        },
        "104": {
            "description": "Attempting to add object of wrong type"
        },
        "105": {
            "description": "Attempting to add object of correct type but wrong qualifiers"
        },
        "106": {
            "description": "Attempting to add data which does not satisfy the qualifiers for a list item"
        },
        "107": {
            "description": "Deleting item will reduce list below minimum length"
        },
        "108": {
            "description": "Adding item will extend list beyond maximum length"
        },
        "109": {
            "description": "Invalid item class"
        },
        "110": {
            "description": "etree (XML) list item of wrong type"
        },
        "112": {
            "description": "No list item object set for list"
        }
    },
    qualifiers={
        "listMinLength": 1,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
    qualifiers_definition={
        "default": {'type': 'list'},
        "listMaxLength": {'type': 'int', 'description': 'Inclusive maximum length of list'},
        "listMinLength": {'type': 'int', 'description': 'Inclusive minimum length of list'},
        "listCompare": {'type': 'int', 'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
    },
)
class CExePathListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CExePathList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExePathListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "name": attribute(AttributeType.STRING),
        "path": attribute(AttributeType.CUSTOM, custom_class="CDataFile"),
    },
    error_codes={
        "0": {
            "severity": 0,
            "description": "OK"
        },
        "1": {
            "severity": 1,
            "description": "Data has undefined value"
        },
        "2": {
            "severity": 3,
            "description": "Data has undefined value"
        },
        "3": {
            "severity": 2,
            "description": "Missing data"
        },
        "4": {
            "description": "Missing data"
        },
        "5": {
            "description": "Attempting to set data of wrong type"
        },
        "6": {
            "description": "Default value does not satisfy validity check"
        },
        "7": {
            "severity": 2,
            "description": "Unrecognised qualifier in data input"
        },
        "8": {
            "severity": 2,
            "description": "Attempting to get inaccessible attribute:"
        },
        "9": {
            "description": "Failed to get property"
        },
        "10": {
            "severity": 2,
            "description": "Attempting to set inaccessible attribute:"
        },
        "11": {
            "description": "Failed to set property:"
        },
        "12": {
            "description": "Undetermined error setting value from XML"
        },
        "13": {
            "description": "Unrecognised class name in qualifier"
        },
        "14": {
            "severity": 2,
            "description": "No object name when saving qualifiers to XML"
        },
        "15": {
            "description": "Error saving qualifier to XML"
        },
        "16": {
            "severity": 2,
            "description": "Unrecognised item in XML data file"
        },
        "17": {
            "description": "Attempting to set unrecognised qualifier"
        },
        "18": {
            "description": "Attempting to set qualifier with wrong type"
        },
        "19": {
            "description": "Attempting to set qualifier with wrong list item type"
        },
        "20": {
            "description": "Error creating a list/dict item object"
        },
        "21": {
            "description": "Unknown error setting qualifiers from Xml file"
        },
        "22": {
            "description": "Unknown error testing validity"
        },
        "23": {
            "description": "Error saving data object to XML"
        },
        "24": {
            "description": "Unable to test validity of default",
            "severity": 2
        },
        "300": {
            "description": "Compared objects are the same",
            "severity": 0
        },
        "315": {
            "description": "Both compared objects are null",
            "severity": 0
        },
        "301": {
            "description": "Unable to compare this class of data",
            "severity": 2
        },
        "302": {
            "description": "Other data has null value"
        },
        "303": {
            "description": "My data has null value"
        },
        "304": {
            "description": "Data has different values"
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
class CSearchPathStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CSearchPath.py
    to add methods and implementation-specific functionality.
    """

    name: Optional[CString] = None
    path: Optional[CDataFile] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSearchPathStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "List shorter than required minimum length"
        },
        "102": {
            "description": "List longer than required maximum length"
        },
        "103": {
            "description": "Consecutive values in list fail comparison test"
        },
        "104": {
            "description": "Attempting to add object of wrong type"
        },
        "105": {
            "description": "Attempting to add object of correct type but wrong qualifiers"
        },
        "106": {
            "description": "Attempting to add data which does not satisfy the qualifiers for a list item"
        },
        "107": {
            "description": "Deleting item will reduce list below minimum length"
        },
        "108": {
            "description": "Adding item will extend list beyond maximum length"
        },
        "109": {
            "description": "Invalid item class"
        },
        "110": {
            "description": "etree (XML) list item of wrong type"
        },
        "112": {
            "description": "No list item object set for list"
        }
    },
    qualifiers={
        "listMinLength": 0,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
    qualifiers_definition={
        "default": {'type': 'list'},
        "listMaxLength": {'type': 'int', 'description': 'Inclusive maximum length of list'},
        "listMinLength": {'type': 'int', 'description': 'Inclusive minimum length of list'},
        "listCompare": {'type': 'int', 'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
    },
)
class CSearchPathListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CSearchPathList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSearchPathListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Invalid characters in file name"
        },
        "102": {
            "severity": 2,
            "description": "Invalid characters in file name"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "allowedCharacters": '',
        "allowedCharactersMode": 1,
        "default": None,
    },
    qualifiers_order=[
        'allowUndefined',
        'allowedCharacters',
        'allowedCharactersMode',
        'default'],
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if allow undefined value at run time'},
        "allowedCharacters": {'type': 'str', 'description': 'Set of characters allowed in file name'},
        "allowedCharactersMode": {'type': 'int', 'description': 'Handling of violation of allowed characters'},
        "default": {'type': 'str', 'description': 'Default file path'},
    },
)
class CFilePathStub(CString):
    """
    A file path

    This is a pure data class stub. Extend it in core/CFilePath.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFilePathStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "function": attribute(AttributeType.CUSTOM, custom_class="CFileFunctionStub"),
        "userId": attribute(AttributeType.CUSTOM, custom_class="CUserIdStub"),
        "hostName": attribute(AttributeType.CUSTOM, custom_class="CHostNameStub"),
        "creationTime": attribute(AttributeType.CUSTOM, custom_class="CTimeStub"),
        "ccp4iVersion": attribute(AttributeType.CUSTOM, custom_class="CVersionStub"),
        "pluginName": attribute(AttributeType.STRING),
        "pluginVersion": attribute(AttributeType.CUSTOM, custom_class="CVersionStub"),
        "pluginTitle": attribute(AttributeType.STRING),
        "projectName": attribute(AttributeType.CUSTOM, custom_class="CProjectNameStub"),
        "projectId": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "jobId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "jobNumber": attribute(AttributeType.STRING),
        "comment": attribute(AttributeType.STRING),
        "OS": attribute(AttributeType.STRING),
    },
    error_codes={
        "101": {
            "description": "Attempting to read header from non-existant Xml file"
        },
        "102": {
            "description": "Error loading file to read header"
        },
        "103": {
            "description": "Error finding <ccp4i2_header> in file"
        },
        "104": {
            "description": "Error interpreting header from file"
        },
        "105": {
            "description": "File does not have <ccp4i2> root node"
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
    content_qualifiers={
        "projectId": {'allowUnfound': True},
    },
)
class CI2XmlHeaderStub(CData):
    """
    Container for header info from XML file

    This is a pure data class stub. Extend it in core/CI2XmlHeader.py
    to add methods and implementation-specific functionality.
    """

    function: Optional[CFileFunctionStub] = None
    userId: Optional[CUserIdStub] = None
    hostName: Optional[CHostNameStub] = None
    creationTime: Optional[CTimeStub] = None
    ccp4iVersion: Optional[CVersionStub] = None
    pluginName: Optional[CString] = None
    pluginVersion: Optional[CVersionStub] = None
    pluginTitle: Optional[CString] = None
    projectName: Optional[CProjectNameStub] = None
    projectId: Optional[CProjectIdStub] = None
    jobId: Optional[CUUIDStub] = None
    jobNumber: Optional[CString] = None
    comment: Optional[CString] = None
    OS: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CI2XmlHeaderStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/grace',
        "fileExtensions": ['xmgr'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CXmgrDataFileStub(CDataFile):
    """
    An xmgr format file. This is the input format for xmgrace, as output by scala or aimless

    This is a pure data class stub. Extend it in core/CXmgrDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXmgrDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "mimeTypeName": '"text/plain"',
        "mimeTypeDescription": 'Standard plain text',
        "fileLabel": None,
        "fileExtensions": ['txt', 'log'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CTextDataFileStub(CDataFile):
    """
    A text data file

    This is a pure data class stub. Extend it in core/CTextDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTextDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "guiLabel": 'Reflections from DIALS',
        "fileExtensions": ['refl'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CDataReflFileStub(CDataFile):
    """
    Reflection file from DIALS

    This is a pure data class stub. Extend it in core/CDataReflFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReflFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/postscript',
        "fileExtensions": ['ps'],
        "guiLabel": 'Postscript file',
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CPostscriptDataFileStub(CDataFile):
    """
    A postscript format file

    This is a pure data class stub. Extend it in core/CPostscriptDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPostscriptDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "mimeTypeName": '"text/plain"',
        "mimeTypeDescription": 'Standard plain text',
        "guiLabel": 'yml file',
        "fileExtensions": ['yml'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CYmlFileStub(CDataFile):
    """
    A yml data file

    This is a pure data class stub. Extend it in core/CYmlFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CYmlFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "1001": {
            "description": "Unknown error reading XML file"
        },
        "1002": {
            "description": "Error trying to find root node in XML"
        },
        "1006": {
            "description": "Attempting to save XML file with incorrect body"
        },
        "1007": {
            "description": "Error creating XML text"
        },
        "1008": {
            "description": "Error saving XML text to file"
        },
        "1009": {
            "description": "Error reading XML file"
        },
        "1010": {
            "description": "XML file does not exist"
        },
        "1011": {
            "description": "No file name given for making I2XMlDataFile"
        },
        "1012": {
            "description": "Error creating I2XMlDataFile object"
        },
        "1013": {
            "description": "Error creating I2XMlDataFile file"
        }
    },
    qualifiers={
        "fileExtensions": ['xml'],
        "saveToDb": False,
        "mimeTypeName": 'application/xml',
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CXmlDataFileStub(CDataFile):
    """
    A reference to an XML file

    This is a pure data class stub. Extend it in core/CXmlDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXmlDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "fileExtensions": ['cif', 'ent'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CMmcifDataFileStub(CDataFile):
    """
    A generic mmCIF format file.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.

    This is a pure data class stub. Extend it in core/CMmcifDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMmcifDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/x-pdf',
        "fileExtensions": ['pdf'],
        "guiLabel": 'PDF file',
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CPDFDataFileStub(CDataFile):
    """
    An PDF format file

    This is a pure data class stub. Extend it in core/CPDFDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPDFDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "File does not exist"
        },
        "102": {
            "description": "No mime type for data file"
        },
        "103": {
            "description": "Attempting to set file content with inappropriate data"
        },
        "104": {
            "description": "There is no file content class specified for this type of file"
        },
        "105": {
            "description": "The file content class specified for this type of file can not be found"
        },
        "300": {
            "description": "Passed",
            "severity": 0
        },
        "305": {
            "description": "Neither original nor test file exists",
            "severity": 0
        },
        "306": {
            "description": "Original file does not exists"
        },
        "307": {
            "description": "Test file does not exist "
        },
        "308": {
            "description": "Files failed checksum comparison"
        },
        "309": {
            "description": "Files failed size comparison"
        },
        "310": {
            "description": "No comparison testing implemented for this file type",
            "severity": 2
        },
        "311": {
            "description": "Failed loading original file for comparison"
        },
        "312": {
            "description": "Failed loading test file for comparison"
        },
        "313": {
            "description": "Files failed simple text diff comparison"
        },
        "320": {
            "description": "Unrecognised error attempting to load file"
        }
    },
    qualifiers={
        "fileLabel": 'scene',
        "mimeTypeName": 'application/CCP4-scene',
        "mimeTypeDescription": 'CCP4mg scene file',
        "guiLabel": 'CCP4mg scene',
        "fileExtensions": ['scene.xml'],
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CSceneDataFileStub(CDataFile):
    """
    An xml format file for defining scene in CCP4mg.

    This is a pure data class stub. Extend it in core/CSceneDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSceneDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "header": attribute(AttributeType.CUSTOM, custom_class="CI2XmlHeaderStub"),
    },
    error_codes={
        "1003": {
            "description": "XML does not have <ccp4i2> root node"
        },
        "1004": {
            "severity": 2,
            "description": "XML does not have <ccp4i2_header> section"
        },
        "1005": {
            "description": "XML does not have <ccp4i2_body> section"
        }
    },
    qualifiers={
        "fileExtensions": ['xml'],
        "autoLoadHeader": True,
    },
    qualifiers_order=['autoLoadHeader'],
    qualifiers_definition={
        "autoLoadHeader": {'type': 'bool'},
    },
)
class CI2XmlDataFileStub(CXmlDataFileStub):
    """
    A reference to an XML file with CCP4i2 Header

    This is a pure data class stub. Extend it in core/CI2XmlDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    header: Optional[CI2XmlHeaderStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CI2XmlDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "project": attribute(AttributeType.CUSTOM, custom_class="CProjectIdStub"),
        "baseName": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "relPath": attribute(AttributeType.CUSTOM, custom_class="CFilePathStub"),
        "annotation": attribute(AttributeType.STRING),
        "dbFileId": attribute(AttributeType.CUSTOM, custom_class="CUUIDStub"),
        "subType": attribute(AttributeType.INT),
        "contentFlag": attribute(AttributeType.INT),
    },
    error_codes={
        "1001": {
            "description": "Unknown error reading XML file"
        },
        "1002": {
            "description": "Error trying to find root node in XML"
        },
        "1006": {
            "description": "Attempting to save XML file with incorrect body"
        },
        "1007": {
            "description": "Error creating XML text"
        },
        "1008": {
            "description": "Error saving XML text to file"
        },
        "1009": {
            "description": "Error reading XML file"
        },
        "1010": {
            "description": "XML file does not exist"
        },
        "1011": {
            "description": "No file name given for making I2XMlDataFile"
        },
        "1012": {
            "description": "Error creating I2XMlDataFile object"
        },
        "1013": {
            "description": "Error creating I2XMlDataFile file"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/EBI-validation-xml',
        "fileExtensions": ['xml'],
        "guiLabel": 'EBI Validation XML',
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
    qualifiers_definition={
        "allowUndefined": {'type': 'bool', 'description': 'Flag if data file can be undefined at run time'},
        "mustExist": {'type': 'bool', 'description': 'Flag if data file must exist at run time'},
        "fromPreviousJob": {'type': 'bool', 'description': 'Flag if input data file can be inferred from preceeding jobs'},
        "jobCombo": {'type': 'bool', 'description': 'Flag if data widget should be a combo box '},
        "mimeTypeName": {'type': 'str', 'description': ''},
        "mimeTypeDescription": {'type': 'str', 'description': ''},
        "fileLabel": {'type': 'str', 'description': 'Label for file'},
        "fileExtensions": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of strings containing allowed file extensions (no dot)'},
        "fileContentClassName": {'type': 'str', 'editable': False, 'description': 'A string containing the name of a class which will hold the file contents'},
        "isDirectory": {'type': 'bool', 'description': 'Flag if the data is a directory'},
        "ifInfo": {'type': 'bool', 'description': 'Flag if gui widget should have info icon'},
        "saveToDb": {'type': 'bool', 'description': 'Save the name of this file in the database'},
        "requiredSubType": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed sub types'},
        "requiredContentFlag": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed content flags'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CEBIValidationXMLDataFileStub(CXmlDataFileStub):
    """
    An XLM file returned from the EBI validation server

    This is a pure data class stub. Extend it in core/CEBIValidationXMLDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEBIValidationXMLDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
