"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.

This is a stub file - extend classes in core/ to add methods.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any

# Metadata system
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType

# Base classes
from ccp4i2.core.base_object.base_classes import CData, CDataFile, CDataFileContent

# Fundamental types
from ccp4i2.core.base_object.fundamental_types import CBoolean, CFloat, CInt, CList, CString

# Cross-file stub class references
from ccp4i2.core.cdata_stubs.CCP4Data import COneWordStub, CRangeSelectionStub, CUUIDStub
from ccp4i2.core.cdata_stubs.CCP4File import CFilePathStub, CMmcifDataStub, CMmcifDataFileStub, CProjectIdStub
from ccp4i2.core.cdata_stubs.CCP4ModelData import CElementStub, CSeqDataFileStub


@cdata_class(
    error_codes={
        "101": {
            "description": "below minimum"
        },
        "102": {
            "description": "above maximum"
        },
        "103": {
            "description": "not one of limited allowed values"
        }
    },
    qualifiers={
        "min": 0.0,
        "default": None,
        "allowUndefined": False,
        "toolTip": 'Cell length in A',
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'],
    qualifiers_definition={
        "default": {'type': 'float'},
        "max": {'description': 'The inclusive maximum value'},
        "min": {'description': 'The inclusive minimum value'},
        "enumerators": {'type': 'list', 'description': 'A Python list of allowed or recommended values - see onlyEnumerators'},
        "menuText": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A Python list of strings, matching items in enumerators list, to appear on GUI menu'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
    },
)
class CCellLengthStub(CFloat):
    """
    A cell length

    This is a pure data class stub. Extend it in core/CCellLength.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCellLengthStub.

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
        "mimeTypeName": 'application/phaser-rfile',
        "mimeTypeDescription": 'Phaser rotation solution file',
        "fileExtensions": ['phaser_rlist.pkl'],
        "fileContentClassName": None,
        "fileLabel": 'phaser_rfile',
        "guiLabel": 'Phaser rotation solution',
        "toolTip": 'Phaser rfile solutions for rotation search',
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
class CPhaserRFileDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CPhaserRFileDataFile.py
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
        Initialize CPhaserRFileDataFileStub.

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
class CColumnGroupListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CColumnGroupList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CColumnGroupListStub.

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
        "guiLabel": 'mmCIF reflection data',
        "mimeTypeName": 'chemical/x-cif',
        "toolTip": 'A reflection file in mmCIF format',
        "fileContentClassName": 'CMmcifReflData',
        "helpFile": 'data_files#mmCIF',
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
class CMmcifReflDataFileStub(CMmcifDataFileStub):
    """
    A reflection file in mmCIF format

    This is a pure data class stub. Extend it in core/CMmcifReflDataFile.py
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
        Initialize CMmcifReflDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "columnGroupType": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "contentFlag": attribute(AttributeType.INT),
        "dataset": attribute(AttributeType.STRING),
        "columnList": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "selected": attribute(AttributeType.BOOLEAN),
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
    content_qualifiers={
        "columnGroupType": {'onlyEnumerators': True, 'enumerators': ['Obs', 'Phs', 'MapCoeffs', 'FreeR']},
    },
)
class CColumnGroupStub(CData):
    """
    Groups of columns in MTZ - probably from analysis by hklfile

    This is a pure data class stub. Extend it in core/CColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    columnGroupType: Optional[COneWordStub] = None
    contentFlag: Optional[CInt] = None
    dataset: Optional[CString] = None
    columnList: Optional[CList] = None
    selected: Optional[CBoolean] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CColumnGroupStub.

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
        "151": {
            "description": "Failed converting MTZ file to alternative format"
        },
        "152": {
            "description": "Failed merging MTZ file - invalid input"
        },
        "153": {
            "description": "Failed merging MTZ files - error running cmtzjoin - see log"
        },
        "154": {
            "description": "Failed merging MTZ files - error running cad - see log"
        },
        "401": {
            "description": "MTZ file header data differs"
        },
        "402": {
            "description": "MTZ file columns differ"
        },
        "403": {
            "description": "Error trying to access number of reflections",
            "severity": 2
        },
        "404": {
            "description": "MTZ files have different number of reflections"
        },
        "405": {
            "description": "MTZ column mean value differs"
        },
        "406": {
            "description": "MTZ file header data differs - may be autogenerated names",
            "severity": 2
        },
        "407": {
            "description": "Error splitting MTZ file - failed creating input command to cmtzsplit"
        },
        "408": {
            "description": "Error splitting MTZ file - output file missing"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz',
        "mimeTypeDescription": 'MTZ experimental data',
        "fileExtensions": ['mtz'],
        "fileContentClassName": 'CMtzData',
        "guiLabel": 'Experimental data',
        "toolTip": "Experimental data in CCP4's MTZ format",
        "helpFile": 'data_files#MTZ',
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
        "sameCrystalAs": {'type': 'str', 'description': 'Name of CMtzDataFile object that crystal parameters should match - probably the observed data'},
        "sameCrystalLevel": {'type': 'int', 'description': 'Rigour of same crystal test'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CMtzDataFileStub(CDataFile):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CMtzDataFile.py
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
        Initialize CMtzDataFileStub.

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
        "onlyEnumerators": True,
        "enumerators": ['native', 'derivative', 'SAD', 'peak', 'inflection', 'high_remote', 'low_remote', ''],
        "default": 'SAD',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CExperimentalDataTypeStub(CString):
    """
    Experimental data type e.g. native or peak

    This is a pure data class stub. Extend it in core/CExperimentalDataType.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExperimentalDataTypeStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "moleculeType": attribute(AttributeType.STRING),
        "seqFile": attribute(AttributeType.CUSTOM, custom_class="CSeqDataFileStub"),
        "numberOfCopies": attribute(AttributeType.INT),
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
    content_qualifiers={
        "moleculeType": {'onlyEnumerators': True, 'enumerators': ['PROTEIN', 'NUCLEIC'], 'menuText': ['protein', 'nucleic acid'], 'default': 'PROTEIN', 'toolTip': 'Molecule type'},
        "seqFile": {'jobCombo': False, 'mustExist': True, 'allowUndefined': False},
        "numberOfCopies": {'allowUndefined': False, 'toolTip': 'Number of copies of sequence', 'min': 0, 'max': 999, 'default': 1, 'enumerators': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]},
    },
)
class CAsuComponentStub(CData):
    """
    A component of the asymmetric unit. This is for use in MR, defining
what we are searching for.

    This is a pure data class stub. Extend it in core/CAsuComponent.py
    to add methods and implementation-specific functionality.
    """

    moleculeType: Optional[CString] = None
    seqFile: Optional[CSeqDataFileStub] = None
    numberOfCopies: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuComponentStub.

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
        "mimeTypeName": 'application/dials-pfile',
        "mimeTypeDescription": 'Dials pickle data file',
        "fileExtensions": ['pickle', 'refl'],
        "fileContentClassName": None,
        "fileLabel": 'dials_pdata',
        "guiLabel": 'Xia2/Dials pickle data',
        "toolTip": 'Xia2/Dials pickle data files',
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
class CDialsPickleFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CDialsPickleFile.py
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
        Initialize CDialsPickleFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "below minimum"
        },
        "102": {
            "description": "above maximum"
        },
        "103": {
            "description": "not one of limited allowed values"
        }
    },
    qualifiers={
        "min": 0.0,
        "max": 180.0,
        "default": None,
        "allowUndefined": True,
        "toolTip": 'Cell angle in degrees',
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'],
    qualifiers_definition={
        "default": {'type': 'float'},
        "max": {'description': 'The inclusive maximum value'},
        "min": {'description': 'The inclusive minimum value'},
        "enumerators": {'type': 'list', 'description': 'A Python list of allowed or recommended values - see onlyEnumerators'},
        "menuText": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A Python list of strings, matching items in enumerators list, to appear on GUI menu'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
    },
)
class CCellAngleStub(CFloat):
    """
    A cell angle

    This is a pure data class stub. Extend it in core/CCellAngle.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCellAngleStub.

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
        "listMinLength": 2,
        "saveToDb": True,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
    qualifiers_definition={
        "default": {'type': 'list'},
        "listMaxLength": {'type': 'int', 'description': 'Inclusive maximum length of list'},
        "listMinLength": {'type': 'int', 'description': 'Inclusive minimum length of list'},
        "listCompare": {'type': 'int', 'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
    },
)
class CMergeMiniMtzListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CMergeMiniMtzList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMergeMiniMtzListStub.

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
        "mimeTypeName": 'application/CCP4-unmerged-experimental',
        "mimeTypeDescription": 'Unmerged experimental data',
        "fileExtensions": ['mtz', 'hkl', 'HKL', 'sca', 'SCA', 'ent', 'cif'],
        "fileContentClassName": 'CUnmergedDataContent',
        "guiLabel": 'Unmerged reflections',
        "toolTip": 'Unmerged experimental data in any format',
        "helpFile": 'data_files#unmerged_data',
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
class CUnmergedDataFileStub(CDataFile):
    """
    Handle MTZ, XDS and scalepack files. Allow wildcard filename

    This is a pure data class stub. Extend it in core/CUnmergedDataFile.py
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
        Initialize CUnmergedDataFileStub.

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
        "mimeTypeName": 'application/dials-jfile',
        "mimeTypeDescription": 'Dials json data file',
        "fileExtensions": ['json', 'expt', 'jsn'],
        "fileContentClassName": None,
        "fileLabel": 'dials_jdata',
        "guiLabel": 'json data',
        "toolTip": 'json data files',
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
class CDialsJsonFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CDialsJsonFile.py
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
        Initialize CDialsJsonFileStub.

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
        "onlyEnumerators": True,
        "default": 'UNDEFINED',
        "enumerators": ['UNDEFINED', 'HREM', 'LREM', 'PEAK', 'INFL', 'NAT', 'DERI'],
        "menuText": ['undefined', 'high remote', 'low remote', 'peak', 'inflection', 'native', 'derivative'],
        "toolTip": 'Hint to Shelx for the use of the dataset',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CShelxLabelStub(CString):
    """
    A string

    This is a pure data class stub. Extend it in core/CShelxLabel.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CShelxLabelStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "h": attribute(AttributeType.STRING),
        "k": attribute(AttributeType.STRING),
        "l": attribute(AttributeType.STRING),
    },
    error_codes={
        "201": {
            "description": "Operator has bad syntax (needs three comma-separated fields)"
        },
        "202": {
            "description": "Operator contains invalid characters"
        },
        "203": {
            "description": "Operator is not set"
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
    contents_order=['h', 'k', 'l'],
    content_qualifiers={
        "h": {'default': 'h'},
        "k": {'default': 'k'},
        "l": {'default': 'l'},
    },
)
class CReindexOperatorStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CReindexOperator.py
    to add methods and implementation-specific functionality.
    """

    h: Optional[CString] = None
    k: Optional[CString] = None
    l: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CReindexOperatorStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "below minimum"
        },
        "102": {
            "description": "above maximum"
        },
        "103": {
            "description": "not one of limited allowed values"
        }
    },
    qualifiers={
        "min": 0.0,
        "toolTip": 'Data collection wavelength in Angstrom',
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'],
    qualifiers_definition={
        "default": {'type': 'float'},
        "max": {'description': 'The inclusive maximum value'},
        "min": {'description': 'The inclusive minimum value'},
        "enumerators": {'type': 'list', 'description': 'A Python list of allowed or recommended values - see onlyEnumerators'},
        "menuText": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A Python list of strings, matching items in enumerators list, to appear on GUI menu'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
    },
)
class CWavelengthStub(CFloat):
    """
    Wavelength in Angstrom

    This is a pure data class stub. Extend it in core/CWavelength.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CWavelengthStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "runNumber": attribute(AttributeType.INT),
        "batchRange0": attribute(AttributeType.INT),
        "batchRange1": attribute(AttributeType.INT),
        "resolution": attribute(AttributeType.FLOAT),
        "fileNumber": attribute(AttributeType.INT),
    },
    error_codes={
        "101": {
            "description": "End of batch range less than start"
        },
        "102": {
            "description": "All items must be set"
        }
    },
    qualifiers={
        "toolTip": 'Specify range of reflections to treat as one run',
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
        "runNumber": {'allowUndefined': True, 'min': 1},
        "batchRange0": {'allowUndefined': True, 'min': 1},
        "batchRange1": {'allowUndefined': True, 'min': 1},
        "resolution": {'min': 0.0, 'allowUndefined': True},
        "fileNumber": {'allowUndefined': True, 'min': 1},
    },
)
class CRunBatchRangeStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CRunBatchRange.py
    to add methods and implementation-specific functionality.
    """

    runNumber: Optional[CInt] = None
    batchRange0: Optional[CInt] = None
    batchRange1: Optional[CInt] = None
    resolution: Optional[CFloat] = None
    fileNumber: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRunBatchRangeStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Invalid space group"
        },
        "102": {
            "description": "Space group is not chiral",
            "severity": 2
        },
        "103": {
            "description": "Space group is not Hermann-Mauguin standard"
        },
        "104": {
            "description": "Space group is not a chiral Hermann-Mauguin standard. Full syminfo.lib information not loaded."
        },
        "105": {
            "description": "Space group is not Hermann-Mauguin standard - has wrong number of spaces?"
        },
        "106": {
            "description": "Space group is undefined",
            "severity": 1
        },
        "107": {
            "description": "Space group is undefined"
        },
        "108": {
            "description": "Space group is incomplete",
            "severity": 2
        }
    },
    qualifiers={
        "allowUndefined": True,
        "toolTip": 'Hermann-Mauguin space group name',
        "helpFile": 'crystal_data#space_group',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CSpaceGroupStub(CString):
    """
    A string holding the space group

    This is a pure data class stub. Extend it in core/CSpaceGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSpaceGroupStub.

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
        "allowUndefined": False,
        "allowedChars": 1,
        "minLength": 1,
        "toolTip": 'Unique identifier for dataset (one word)',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CDatasetNameStub(CString):
    """
    A string

    This is a pure data class stub. Extend it in core/CDatasetName.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDatasetNameStub.

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
class CMiniMtzDataFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CMiniMtzDataFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMiniMtzDataFileListStub.

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
class CUnmergedDataFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CUnmergedDataFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUnmergedDataFileListStub.

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
        "enumerators": ['H', 'J', 'F', 'D', 'Q', 'G', 'L', 'K', 'M', 'E', 'P', 'W', 'A', 'B', 'Y', 'I', 'R'],
        "onlyEnumerators": True,
        "default": 'F',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CColumnTypeStub(CString):
    """
    A list of recognised MTZ column types

    This is a pure data class stub. Extend it in core/CColumnType.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CColumnTypeStub.

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
class CColumnTypeListStub(CList):
    """
    A list of acceptable MTZ column types

    This is a pure data class stub. Extend it in core/CColumnTypeList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CColumnTypeListStub.

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
        "mimeTypeName": 'application/phaser-sol',
        "mimeTypeDescription": 'Phaser solution file',
        "fileExtensions": ['phaser_sol.pkl'],
        "fileContentClassName": None,
        "fileLabel": 'phaser_sol',
        "guiLabel": 'Phaser solutions',
        "toolTip": 'Possible solutions passed between runs of the Phaser program',
        "helpFile": 'data_files#phasersol',
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
class CPhaserSolDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CPhaserSolDataFile.py
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
        Initialize CPhaserSolDataFileStub.

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
class CImportUnmergedListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CImportUnmergedList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CImportUnmergedListStub.

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
        "mimeTypeName": 'application/CCP4-shelx-FA',
        "mimeTypeDescription": 'Shelx FA',
        "fileExtensions": ['hkl'],
        "fileContentClassName": None,
        "fileLabel": 'shelx_FA',
        "guiLabel": 'Shelx FA',
        "toolTip": 'Data used by Shelx programs',
        "helpFile": 'data_files#shelxfa',
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
class CShelxFADataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CShelxFADataFile.py
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
        Initialize CShelxFADataFileStub.

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
class CAltSpaceGroupListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CAltSpaceGroupList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAltSpaceGroupListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "name": attribute(AttributeType.STRING),
        "columnGroups": attribute(AttributeType.CUSTOM, custom_class="CList"),
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
class CMtzDatasetStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CMtzDataset.py
    to add methods and implementation-specific functionality.
    """

    name: Optional[CString] = None
    columnGroups: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMtzDatasetStub.

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
class CImageFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CImageFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CImageFileListStub.

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
        "guiLabel": 'Contents of asymmetric unit',
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
    qualifiers_definition={
        "default": {'type': 'list'},
        "listMaxLength": {'type': 'int', 'description': 'Inclusive maximum length of list'},
        "listMinLength": {'type': 'int', 'description': 'Inclusive minimum length of list'},
        "listCompare": {'type': 'int', 'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
    },
)
class CAsuComponentListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CAsuComponentList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuComponentListStub.

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
        "guiLabel": 'Reflection data',
        "mimeTypeName": 'application/CCP4-generic-reflections',
        "toolTip": 'A reflection data file in MTZ or a non-CCP4 format',
        "fileContentClassName": 'CUnmergedDataContent',
        "fileExtensions": ['mtz', 'hkl', 'HKL', 'sca', 'SCA', 'mmcif', 'cif', 'ent'],
        "downloadModes": ['ebiSFs'],
        "helpFile": 'import_merged#file_formats',
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
class CGenericReflDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CGenericReflDataFile.py
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
        Initialize CGenericReflDataFileStub.

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
class CXia2ImageSelectionListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CXia2ImageSelectionList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXia2ImageSelectionListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "Fp": attribute(AttributeType.FLOAT),
        "Fpp": attribute(AttributeType.FLOAT),
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
    contents_order=['Fp', 'Fpp'],
    content_qualifiers={
        "Fp": {'toolTip': "Form factor F' for element at given wavelength"},
        "Fpp": {'toolTip': "Form factor F'' for element at given wavelength"},
    },
)
class CFormFactorStub(CData):
    """
    The for factor (Fp and Fpp) for a giving element and wavelength

    This is a pure data class stub. Extend it in core/CFormFactor.py
    to add methods and implementation-specific functionality.
    """

    Fp: Optional[CFloat] = None
    Fpp: Optional[CFloat] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFormFactorStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "201": {
            "description": "Word contains white space item"
        }
    },
    qualifiers={
        "onlyEnumerators": False,
        "enumerators": ['Br', 'Fe', 'Pt', 'Se'],
        "charWidth": 4,
        "default": 'Se',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CAnomalousScatteringElementStub(CElementStub):
    """
    Definition of a anomalous scattering element

    This is a pure data class stub. Extend it in core/CAnomalousScatteringElement.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnomalousScatteringElementStub.

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
        "mimeTypeName": 'application/CCP4-image',
        "mimeTypeDescription": 'Image file',
        "fileExtensions": ['img', 'cbf', 'mccd', 'mar1600', 'h5', 'nxs'],
        "fileContentClassName": None,
        "guiLabel": 'Image file',
        "toolTip": 'First image file in a directory',
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
class CImageFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CImageFile.py
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
        Initialize CImageFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "mustExist": False,
        "mtzFileKey": '',
        "toolTipList": [],
        "default": [],
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CProgramColumnGroupStub(CData):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CProgramColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProgramColumnGroupStub.

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
        "mimeTypeName": 'application/refmac-keywords',
        "mimeTypeDescription": 'Refmac keyword file',
        "fileExtensions": ['txt'],
        "fileContentClassName": None,
        "fileLabel": 'refmac_keywords',
        "guiLabel": 'Refmac keyword file',
        "toolTip": 'A file containing keywords as they are meant to be read by refmac5',
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
class CRefmacKeywordFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CRefmacKeywordFile.py
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
        Initialize CRefmacKeywordFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "low": attribute(AttributeType.FLOAT),
        "high": attribute(AttributeType.FLOAT),
    },
    error_codes={
        "201": {
            "description": "High/low resolution wrong way round?"
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
        "low": {'min': 0.0, 'allowUndefined': True},
        "high": {'min': 0.0, 'allowUndefined': True},
    },
)
class CResolutionRangeStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CResolutionRange.py
    to add methods and implementation-specific functionality.
    """

    low: Optional[CFloat] = None
    high: Optional[CFloat] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CResolutionRangeStub.

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
        "mimeTypeName": 'application/CCP4-map',
        "mimeTypeDescription": 'Map',
        "fileExtensions": ['map', 'mrc'],
        "fileContentClassName": None,
        "guiLabel": 'Map',
        "toolTip": 'A map in CCP4/MRC format',
        "helpFile": 'data_files#map_files',
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
class CMapDataFileStub(CDataFile):
    """
    A CCP4 Map file

    This is a pure data class stub. Extend it in core/CMapDataFile.py
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
        Initialize CMapDataFileStub.

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
class CRunBatchRangeListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CRunBatchRangeList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRunBatchRangeListStub.

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
class CDatasetListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CDatasetList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDatasetListStub.

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
        "fileLabel": 'imosflm',
        "mimeTypeName": 'application/iMosflm-xml',
        "mimeTypeDescription": 'iMosflm data',
        "guiLabel": 'iMosflm data',
        "fileExtensions": ['imosflm.xml'],
        "fileContentClassName": None,
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
class CImosflmXmlDataFileStub(CDataFile):
    """
    An iMosflm data file

    This is a pure data class stub. Extend it in core/CImosflmXmlDataFile.py
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
        Initialize CImosflmXmlDataFileStub.

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
        "allowUndefined": False,
        "minLength": 1,
        "allowedChars": 1,
        "toolTip": 'Unique identifier for crystal (one word)',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CCrystalNameStub(CString):
    """
    A string

    This is a pure data class stub. Extend it in core/CCrystalName.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCrystalNameStub.

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
        "201": {
            "description": "Wrong number of columns"
        },
        "202": {
            "description": "Wrong column types"
        },
        "203": {
            "description": "No correct column types found in file"
        },
        "204": {
            "description": "Duplicate or additional column types found in file"
        },
        "205": {
            "description": "Columns in file have non-standard labels"
        },
        "206": {
            "description": "File contains unmerged data"
        },
        "210": {
            "description": "Failed creating mini-MTZ"
        },
        "211": {
            "description": "Insufficient columns selected from imported MTZ"
        },
        "212": {
            "description": "Data already imported as",
            "severity": 2
        },
        "220": {
            "description": "Can not convert file content, file does not exist"
        },
        "221": {
            "description": "Can not convert file content, existing content insufficiently rich"
        },
        "222": {
            "description": "Can not convert file content, bad input for target content"
        },
        "223": {
            "description": "Can not recognise file content"
        },
        "224": {
            "description": "Not possible to convert to required content - no mechanism implemented"
        },
        "225": {
            "description": "Failed importing from an mmcif file - failed running cif2mtz"
        },
        "226": {
            "description": "Failed importing from an mmcif file - no output from cif2mtz"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-mini',
        "fileExtensions": ['mtz', 'cif', 'ent'],
        "fileContentClassName": 'CMtzData',
        "saveToDb": True,
        "correctColumns": ['FQ', 'JQ', 'GLGL', 'KMKM', 'AAAA', 'PW', 'FP', 'I'],
        "toolTip": 'Mini-MTZ file containing reflection,phases,FreeR set or map coefficients',
        "helpFile": 'data_files#MTZ',
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag',
        'correctColumns',
        'columnGroupClassList',
        'sameCrystalAs'],
    qualifiers_definition={
        "correctColumns": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of coloumn data types expected in the file'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CMiniMtzDataFileStub(CMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CMiniMtzDataFile.py
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
        Initialize CMiniMtzDataFileStub.

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
        "151": {
            "description": "Failed converting MTZ file to alternative format"
        },
        "152": {
            "description": "Failed merging MTZ file - invalid input"
        },
        "153": {
            "description": "Failed merging MTZ files - error running cmtzjoin - see log"
        },
        "154": {
            "description": "Failed merging MTZ files - error running cad - see log"
        },
        "401": {
            "description": "MTZ file header data differs"
        },
        "402": {
            "description": "MTZ file columns differ"
        },
        "403": {
            "description": "Error trying to access number of reflections",
            "severity": 2
        },
        "404": {
            "description": "MTZ files have different number of reflections"
        },
        "405": {
            "description": "MTZ column mean value differs"
        },
        "406": {
            "description": "MTZ file header data differs - may be autogenerated names",
            "severity": 2
        },
        "407": {
            "description": "Error splitting MTZ file - failed creating input command to cmtzsplit"
        },
        "408": {
            "description": "Error splitting MTZ file - output file missing"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-unmerged',
        "mimeTypeDescription": 'MTZ unmerged experimental data',
        "fileExtensions": ['mtz'],
        "fileContentClassName": None,
        "guiLabel": 'Unmerged MTZ reflections',
        "toolTip": "Unmerged experimental data in CCP4's MTZ format",
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
        "sameCrystalAs": {'type': 'str', 'description': 'Name of CMtzDataFile object that crystal parameters should match - probably the observed data'},
        "sameCrystalLevel": {'type': 'int', 'description': 'Rigour of same crystal test'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CUnmergedMtzDataFileStub(CMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CUnmergedMtzDataFile.py
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
        Initialize CUnmergedMtzDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "a": attribute(AttributeType.CUSTOM, custom_class="CCellLengthStub"),
        "b": attribute(AttributeType.CUSTOM, custom_class="CCellLengthStub"),
        "c": attribute(AttributeType.CUSTOM, custom_class="CCellLengthStub"),
        "alpha": attribute(AttributeType.CUSTOM, custom_class="CCellAngleStub"),
        "beta": attribute(AttributeType.CUSTOM, custom_class="CCellAngleStub"),
        "gamma": attribute(AttributeType.CUSTOM, custom_class="CCellAngleStub"),
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
        "toolTip": 'Cell lengths and angles',
        "helpFile": 'crystal_data#cell',
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
    contents_order=['a', 'b', 'c', 'alpha', 'beta', 'gamma'],
    content_qualifiers={
        "a": {'toolTip': 'Cell length a in A', 'guiLabel': 'a'},
        "b": {'toolTip': 'Cell length b in A', 'guiLabel': 'b'},
        "c": {'toolTip': 'Cell length c in A', 'guiLabel': 'c'},
        "alpha": {'toolTip': 'Cell angle alpha in degrees', 'guiLabel': 'alpha'},
        "beta": {'toolTip': 'Cell angle beta in degrees', 'guiLabel': 'beta'},
        "gamma": {'toolTip': 'Cell angle gamma in degrees', 'guiLabel': 'gamma'},
    },
)
class CCellStub(CData):
    """
    A unit cell

    This is a pure data class stub. Extend it in core/CCell.py
    to add methods and implementation-specific functionality.
    """

    a: Optional[CCellLengthStub] = None
    b: Optional[CCellLengthStub] = None
    c: Optional[CCellLengthStub] = None
    alpha: Optional[CCellAngleStub] = None
    beta: Optional[CCellAngleStub] = None
    gamma: Optional[CCellAngleStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CCellStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Invalid space group"
        },
        "102": {
            "description": "Space group is not chiral",
            "severity": 2
        },
        "103": {
            "description": "Space group is not Hermann-Mauguin standard"
        },
        "104": {
            "description": "Space group is not a chiral Hermann-Mauguin standard. Full syminfo.lib information not loaded."
        },
        "105": {
            "description": "Space group is not Hermann-Mauguin standard - has wrong number of spaces?"
        },
        "106": {
            "description": "Space group is undefined",
            "severity": 1
        },
        "107": {
            "description": "Space group is undefined"
        },
        "108": {
            "description": "Space group is incomplete",
            "severity": 2
        }
    },
    qualifiers={
        "allowUndefined": True,
        "toolTip": 'Hermann-Mauguin space group name',
        "helpFile": 'crystal_data#space_group',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CAltSpaceGroupStub(CSpaceGroupStub):
    """
    A string holding the space group

    This is a pure data class stub. Extend it in core/CAltSpaceGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAltSpaceGroupStub.

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
        "enumerators": ['H', 'J', 'F', 'D', 'Q', 'G', 'L', 'K', 'M', 'E', 'P', 'W', 'A', 'B', 'Y', 'I', 'R'],
        "onlyEnumerators": True,
        "default": 'F',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
    qualifiers_definition={
        "default": {'type': 'str'},
        "maxLength": {'type': 'int', 'description': 'Maximum length of string'},
        "minLength": {'type': 'int', 'description': 'Minimum length of string'},
        "enumerators": {'type': 'list', 'description': 'A list of allowed or recommended values for string'},
        "menuText": {'type': 'list', 'description': 'A list of strings equivalent to the enumerators that will appear in the GUI'},
        "onlyEnumerators": {'type': 'bool', 'description': 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
        "allowedCharsCode": {'type': 'int', 'description': 'Flag if the text is limited to set of allowed characters'},
    },
)
class CMtzColumnGroupTypeStub(CColumnTypeStub):
    """
    A list of recognised MTZ column types

    This is a pure data class stub. Extend it in core/CMtzColumnGroupType.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMtzColumnGroupTypeStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "columnLabel": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "columnType": attribute(AttributeType.CUSTOM, custom_class="CColumnTypeStub"),
        "dataset": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "groupIndex": attribute(AttributeType.INT),
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
    content_qualifiers={
        "columnLabel": {'allowUndefined': True},
    },
)
class CMtzColumnStub(CData):
    """
    An MTZ column with column label and column type

    This is a pure data class stub. Extend it in core/CMtzColumn.py
    to add methods and implementation-specific functionality.
    """

    columnLabel: Optional[COneWordStub] = None
    columnType: Optional[CColumnTypeStub] = None
    dataset: Optional[COneWordStub] = None
    groupIndex: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMtzColumnStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "columnName": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "defaultList": attribute(AttributeType.STRING),
        "columnType": attribute(AttributeType.CUSTOM, custom_class="CColumnTypeListStub"),
        "partnerTo": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "partnerOffset": attribute(AttributeType.INT),
    },
    error_codes={
        "1": {
            "description": "Attempting to change immutable object"
        },
        "2": {
            "description": "Attempting to access unknown attribute"
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
class CColumnGroupItemStub(CData):
    """
    Definition of set of columns that form a 'group'

    This is a pure data class stub. Extend it in core/CColumnGroupItem.py
    to add methods and implementation-specific functionality.
    """

    columnName: Optional[COneWordStub] = None
    defaultList: Optional[CString] = None
    columnType: Optional[CColumnTypeListStub] = None
    partnerTo: Optional[COneWordStub] = None
    partnerOffset: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CColumnGroupItemStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "imageFile": attribute(AttributeType.CUSTOM, custom_class="CImageFileStub"),
        "imageStart": attribute(AttributeType.INT),
        "imageEnd": attribute(AttributeType.INT),
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
        "toolTip": 'select an image file and an optional range of files to define a dataset',
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
    contents_order=['imageFile', 'imageStart', 'imageEnd'],
    content_qualifiers={
        "imageStart": {'allowUndefined': True, 'min': 0},
        "imageEnd": {'allowUndefined': True, 'min': 0},
    },
)
class CXia2ImageSelectionStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CXia2ImageSelection.py
    to add methods and implementation-specific functionality.
    """

    imageFile: Optional[CImageFileStub] = None
    imageStart: Optional[CInt] = None
    imageEnd: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXia2ImageSelectionStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Anomalous structure factors and sigma',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CFPairColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CFPairColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFPairColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Anomalous intensities and sigma',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CIPairColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CIPairColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CIPairColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Structure factor and phase to define a map',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CMapColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CMapColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMapColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Hendrickson-Lattmann coefficients',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CHLColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CHLColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHLColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Set of FreeR flags',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CFreeRColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CFreeRColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFreeRColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Structure factor and sigma',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CFSigFColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CFSigFColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFSigFColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "toolTipList": ['The real part of the experimental intensity', 'The anomalous part of the experimental intensity'],
        "guiLabel": 'Intensity and anomalous intensity',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CAnomalousIntensityColumnGroupStub(CProgramColumnGroupStub):
    """
    Selection of I and AnomI columns from MTZ.
Expected to be part of ab initio phasing dataset ( CDataset)

    This is a pure data class stub. Extend it in core/CAnomalousIntensityColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnomalousIntensityColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Intensity and sigma',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CISigIColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CISigIColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CISigIColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "toolTipList": ['The real part of the experimental structure factors', 'The anomalous part of the experimental structure factors'],
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CAnomalousColumnGroupStub(CProgramColumnGroupStub):
    """
    Selection of F/I and AnomF/I columns from MTZ.
Expected to be part of ab initio phasing dataset ( CDataset)

    This is a pure data class stub. Extend it in core/CAnomalousColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnomalousColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "Error setting columnGroup qualifier"
        },
        "104": {
            "description": "Missing column selection"
        },
        "105": {
            "description": "Specified column not found in MTZ file"
        },
        "106": {
            "description": "Specified column has wrong type in MTZ file"
        },
        "107": {
            "description": "Error reading columnGroup qualifier from XML file"
        },
        "108": {
            "description": "No columnGroup qualifier"
        }
    },
    qualifiers={
        "guiLabel": 'Phase and figure of merit',
    },
    qualifiers_order=['mtzFileKey', 'mustExist', 'toolTipList', 'default'],
    qualifiers_definition={
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
        "toolTipList": {'type': 'list', 'description': 'Tooltips for columns in group'},
        "default": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'Preferred values for column names'},
    },
)
class CPhiFomColumnGroupStub(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input

    This is a pure data class stub. Extend it in core/CPhiFomColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPhiFomColumnGroupStub.

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
        "201": {
            "description": "Wrong number of columns"
        },
        "202": {
            "description": "Wrong column types"
        },
        "203": {
            "description": "No correct column types found in file"
        },
        "204": {
            "description": "Duplicate or additional column types found in file"
        },
        "205": {
            "description": "Columns in file have non-standard labels"
        },
        "206": {
            "description": "File contains unmerged data"
        },
        "210": {
            "description": "Failed creating mini-MTZ"
        },
        "211": {
            "description": "Insufficient columns selected from imported MTZ"
        },
        "212": {
            "description": "Data already imported as",
            "severity": 2
        },
        "220": {
            "description": "Can not convert file content, file does not exist"
        },
        "221": {
            "description": "Can not convert file content, existing content insufficiently rich"
        },
        "222": {
            "description": "Can not convert file content, bad input for target content"
        },
        "223": {
            "description": "Can not recognise file content"
        },
        "224": {
            "description": "Not possible to convert to required content - no mechanism implemented"
        },
        "225": {
            "description": "Failed importing from an mmcif file - failed running cif2mtz"
        },
        "226": {
            "description": "Failed importing from an mmcif file - no output from cif2mtz"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-map',
        "mimeTypeDescription": 'MTZ F-phi',
        "fileExtensions": ['mtz', 'cif', 'ent'],
        "fileContentClassName": 'CMtzData',
        "fileLabel": 'map_coefficients',
        "guiLabel": 'Map coefficients',
        "toolTip": 'Electron density map coefficients: F,Phi',
        "correctColumns": ['FP', 'FQP'],
        "columnGroupClassList": ["<class 'ccp4x.data_scan.CCP4XtalData.CMapColumnGroup'>"],
        "downloadModes": ['Uppsala-EDS'],
        "helpFile": 'data_files#MapCoeffs',
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag',
        'correctColumns',
        'columnGroupClassList',
        'sameCrystalAs'],
    qualifiers_definition={
        "correctColumns": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of coloumn data types expected in the file'},
    },
    content_qualifiers={
        "subType": {'default': 1, 'enumerators': [1, 2, 3], 'onlyEnumerators': True, 'menuText': ['normal map', 'difference map', 'anomalous difference map']},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CMapCoeffsDataFileStub(CMiniMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CMapCoeffsDataFile.py
    to add methods and implementation-specific functionality.
    """

    # Subtype constants
    SUBTYPE_NORMAL = 1  # normal map
    SUBTYPE_DIFFERENCE = 2  # difference map
    SUBTYPE_ANOM_DIFFERENCE = 3  # anomalous difference map

    # Content flag constants
    CONTENT_FLAG_FPHI = 1  # FPhi

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['FPhi']

    # Column signatures for each content flag (indexed by contentFlag - 1)
    CONTENT_SIGNATURE_LIST = [['F', 'PHI']]

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMapCoeffsDataFileStub.

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
        "301": {
            "description": "Running ctruncate failed"
        },
        "302": {
            "description": "Running cmtzsplit to convert observed data type failed"
        },
        "303": {
            "description": "Running sftools failed"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-observed',
        "mimeTypeDescription": 'MTZ observed',
        "fileExtensions": ['mtz', 'cif', 'ent'],
        "fileContentClassName": 'CMtzData',
        "fileLabel": 'observed_data',
        "guiLabel": 'Reflections',
        "toolTip": 'Observed structure factors or intensities',
        "correctColumns": ['KMKM', 'GLGL', 'JQ', 'FQ'],
        "columnGroupClassList": ["<class 'ccp4x.data_scan.CCP4XtalData.CIPairColumnGroup'>", "<class 'ccp4x.data_scan.CCP4XtalData.CFPairColumnGroup'>", "<class 'ccp4x.data_scan.CCP4XtalData.CISigIColumnGroup'>", "<class 'ccp4x.data_scan.CCP4XtalData.CFSigFColumnGroup'>"],
        "downloadModes": ['ebiSFs'],
        "helpFile": 'data_files#Obs',
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag',
        'correctColumns',
        'columnGroupClassList',
        'sameCrystalAs'],
    qualifiers_definition={
        "correctColumns": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of coloumn data types expected in the file'},
    },
    content_qualifiers={
        "subType": {'default': 1, 'enumerators': [1, 2, 3], 'onlyEnumerators': True, 'menuText': ['observed data', 'derived data', 'reference data']},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CObsDataFileStub(CMiniMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CObsDataFile.py
    to add methods and implementation-specific functionality.
    """

    # Subtype constants
    SUBTYPE_OBSERVED = 1  # observed data
    SUBTYPE_DERIVED = 2  # derived data
    SUBTYPE_REFERENCE = 3  # reference data

    # Content flag constants
    CONTENT_FLAG_IPAIR = 1  # Anomalous Is
    CONTENT_FLAG_FPAIR = 2  # Anomalous SFs
    CONTENT_FLAG_IMEAN = 3  # Mean Is
    CONTENT_FLAG_FMEAN = 4  # Mean SFs

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = [
        'Anomalous Is',
        'Anomalous SFs',
        'Mean Is',
        'Mean SFs']

    # Column signatures for each content flag (indexed by contentFlag - 1)
    CONTENT_SIGNATURE_LIST = [['Iplus', 'SIGIplus', 'Iminus', 'SIGIminus'], [
        'Fplus', 'SIGFplus', 'Fminus', 'SIGFminus'], ['I', 'SIGI'], ['F', 'SIGF']]

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CObsDataFileStub.

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
        "201": {
            "description": "Wrong number of columns"
        },
        "202": {
            "description": "Wrong column types"
        },
        "203": {
            "description": "No correct column types found in file"
        },
        "204": {
            "description": "Duplicate or additional column types found in file"
        },
        "205": {
            "description": "Columns in file have non-standard labels"
        },
        "206": {
            "description": "File contains unmerged data"
        },
        "210": {
            "description": "Failed creating mini-MTZ"
        },
        "211": {
            "description": "Insufficient columns selected from imported MTZ"
        },
        "212": {
            "description": "Data already imported as",
            "severity": 2
        },
        "220": {
            "description": "Can not convert file content, file does not exist"
        },
        "221": {
            "description": "Can not convert file content, existing content insufficiently rich"
        },
        "222": {
            "description": "Can not convert file content, bad input for target content"
        },
        "223": {
            "description": "Can not recognise file content"
        },
        "224": {
            "description": "Not possible to convert to required content - no mechanism implemented"
        },
        "225": {
            "description": "Failed importing from an mmcif file - failed running cif2mtz"
        },
        "226": {
            "description": "Failed importing from an mmcif file - no output from cif2mtz"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-phases',
        "mimeTypeDescription": 'MTZ phases',
        "fileExtensions": ['mtz', 'cif', 'ent'],
        "fileContentClassName": 'CMtzData',
        "guiLabel": 'Phases',
        "fileLabel": 'phases',
        "toolTip": 'Phases in Hendrickson-Lattmann or Phi/FOM form',
        "correctColumns": ['AAAA', 'PW'],
        "columnGroupClassList": ["<class 'ccp4x.data_scan.CCP4XtalData.CHLColumnGroup'>", "<class 'ccp4x.data_scan.CCP4XtalData.CPhiFomColumnGroup'>"],
        "helpFile": 'data_files#Phs',
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag',
        'correctColumns',
        'columnGroupClassList',
        'sameCrystalAs'],
    qualifiers_definition={
        "correctColumns": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of coloumn data types expected in the file'},
    },
    content_qualifiers={
        "subType": {'default': 1, 'enumerators': [1, 2], 'onlyEnumerators': True, 'menuText': ['unbiased data', 'biased data']},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CPhsDataFileStub(CMiniMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CPhsDataFile.py
    to add methods and implementation-specific functionality.
    """

    # Subtype constants
    SUBTYPE_UNBIASED = 1  # unbiased data
    SUBTYPE_BIASED = 2  # biased data

    # Content flag constants
    CONTENT_FLAG_HL = 1  # Hendrickson-Lattmann coeffs
    CONTENT_FLAG_PHIFOM = 2  # Phi,FOM

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['Hendrickson-Lattmann coeffs', 'Phi,FOM']

    # Column signatures for each content flag (indexed by contentFlag - 1)
    CONTENT_SIGNATURE_LIST = [['HLA', 'HLB', 'HLC', 'HLD'], ['PHI', 'FOM']]

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPhsDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "fileName": attribute(AttributeType.CUSTOM, custom_class="CMiniMtzDataFileStub"),
        "columnTag": attribute(AttributeType.STRING),
        "columnNames": attribute(AttributeType.STRING),
    },
    error_codes={
        "201": {
            "description": "Selected file is not a suitable 'mini' MTZ containing experimental data object"
        },
        "202": {
            "description": "Output column name list does not have correct number of names"
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
    contents_order=['fileName', 'columnTag', 'columnNames'],
    content_qualifiers={
        "fileName": {'fromPreviousJob': False},
    },
)
class CMergeMiniMtzStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CMergeMiniMtz.py
    to add methods and implementation-specific functionality.
    """

    fileName: Optional[CMiniMtzDataFileStub] = None
    columnTag: Optional[CString] = None
    columnNames: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMergeMiniMtzStub.

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
        "201": {
            "description": "Wrong number of columns"
        },
        "202": {
            "description": "Wrong column types"
        },
        "203": {
            "description": "No correct column types found in file"
        },
        "204": {
            "description": "Duplicate or additional column types found in file"
        },
        "205": {
            "description": "Columns in file have non-standard labels"
        },
        "206": {
            "description": "File contains unmerged data"
        },
        "210": {
            "description": "Failed creating mini-MTZ"
        },
        "211": {
            "description": "Insufficient columns selected from imported MTZ"
        },
        "212": {
            "description": "Data already imported as",
            "severity": 2
        },
        "220": {
            "description": "Can not convert file content, file does not exist"
        },
        "221": {
            "description": "Can not convert file content, existing content insufficiently rich"
        },
        "222": {
            "description": "Can not convert file content, bad input for target content"
        },
        "223": {
            "description": "Can not recognise file content"
        },
        "224": {
            "description": "Not possible to convert to required content - no mechanism implemented"
        },
        "225": {
            "description": "Failed importing from an mmcif file - failed running cif2mtz"
        },
        "226": {
            "description": "Failed importing from an mmcif file - no output from cif2mtz"
        }
    },
    qualifiers={
        "mimeTypeName": 'application/CCP4-mtz-freerflag',
        "mimeTypeDescription": 'FreeR flag',
        "fileExtensions": ['mtz', 'cif', 'ent'],
        "fileContentClassName": 'CMtzData',
        "fileLabel": 'freeRflag',
        "guiLabel": 'Free R set',
        "toolTip": 'Set of reflections used for FreeR calculation',
        "correctColumns": ['I'],
        "columnGroupClassList": ["<class 'ccp4x.data_scan.CCP4XtalData.CFreeRColumnGroup'>"],
        "helpFile": 'data_files#FreeR',
    },
    qualifiers_order=[
        'fileExtensions',
        'mimeTypeName',
        'mimeTypeDescription',
        'allowUndefined',
        'mustExist',
        'fromPreviousJob',
        'jobCombo',
        'fileContentClassName',
        'isDirectory',
        'saveToDb',
        'requiredSubType',
        'requiredContentFlag',
        'correctColumns',
        'columnGroupClassList',
        'sameCrystalAs'],
    qualifiers_definition={
        "correctColumns": {'type': 'list', 'listItemType': "<class 'str'>", 'description': 'A list of coloumn data types expected in the file'},
    },
    content_qualifiers={
        "subType": {'enumerators': [], 'onlyEnumerators': True},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CFreeRDataFileStub(CMiniMtzDataFileStub):
    """
    An MTZ experimental data file

    This is a pure data class stub. Extend it in core/CFreeRDataFile.py
    to add methods and implementation-specific functionality.
    """

    # Content flag constants
    CONTENT_FLAG_FREER = 1  # FreeR

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['FreeR']

    # Column signatures for each content flag (indexed by contentFlag - 1)
    CONTENT_SIGNATURE_LIST = [['FREER']]

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFreeRDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "spaceGroup": attribute(AttributeType.CUSTOM, custom_class="CSpaceGroupStub"),
        "cell": attribute(AttributeType.CUSTOM, custom_class="CCellStub"),
    },
    error_codes={
        "101": {
            "description": "Cell lengths should NOT be identical"
        },
        "102": {
            "description": "Cell angles should NOT be identical"
        },
        "103": {
            "description": "Cell angle should be 90"
        },
        "104": {
            "description": "Cell angle should NOT be 90"
        },
        "105": {
            "description": "Cell lengths should be identical"
        },
        "106": {
            "description": "Cell angle should be 120"
        },
        "107": {
            "description": "Cell angle should be identical"
        }
    },
    qualifiers={
        "toolTip": 'Space group and cell length and angles',
        "helpFile": 'crystal_data#cell_space_group',
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
    contents_order=['spaceGroup', 'cell'],
    content_qualifiers={
        "spaceGroup": {'guilabel': 'space group'},
        "cell": {'guilabel': 'cell'},
    },
)
class CSpaceGroupCellStub(CData):
    """
    Cell space group and parameters

    This is a pure data class stub. Extend it in core/CSpaceGroupCell.py
    to add methods and implementation-specific functionality.
    """

    spaceGroup: Optional[CSpaceGroupStub] = None
    cell: Optional[CCellStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSpaceGroupCellStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "cell": attribute(AttributeType.CUSTOM, custom_class="CCellStub"),
        "spaceGroup": attribute(AttributeType.CUSTOM, custom_class="CSpaceGroupStub"),
        "wavelength": attribute(AttributeType.CUSTOM, custom_class="CWavelengthStub"),
        "haveFreeRColumn": attribute(AttributeType.BOOLEAN),
        "haveFobsColumn": attribute(AttributeType.BOOLEAN),
        "haveFpmObsColumn": attribute(AttributeType.BOOLEAN),
        "haveIobsColumn": attribute(AttributeType.BOOLEAN),
        "haveIpmObsColumn": attribute(AttributeType.BOOLEAN),
    },
    error_codes={
        "101": {
            "description": "Attempting to load mmCIF data from non-existant/broken file"
        },
        "102": {
            "description": "Error reading interpreting line in cif file"
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
class CMmcifReflDataStub(CMmcifDataStub):
    """
    Generic mmCIF data.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.

    This is a pure data class stub. Extend it in core/CMmcifReflData.py
    to add methods and implementation-specific functionality.
    """

    cell: Optional[CCellStub] = None
    spaceGroup: Optional[CSpaceGroupStub] = None
    wavelength: Optional[CWavelengthStub] = None
    haveFreeRColumn: Optional[CBoolean] = None
    haveFobsColumn: Optional[CBoolean] = None
    haveFpmObsColumn: Optional[CBoolean] = None
    haveIobsColumn: Optional[CBoolean] = None
    haveIpmObsColumn: Optional[CBoolean] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMmcifReflDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "file": attribute(AttributeType.CUSTOM, custom_class="CUnmergedDataFileStub"),
        "cell": attribute(AttributeType.CUSTOM, custom_class="CCellStub"),
        "wavelength": attribute(AttributeType.CUSTOM, custom_class="CWavelengthStub"),
        "crystalName": attribute(AttributeType.STRING),
        "dataset": attribute(AttributeType.STRING),
        "excludeSelection": attribute(AttributeType.CUSTOM, custom_class="CRangeSelectionStub"),
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
        "toolTip": 'Imported data file, cell parameters and crystal/dataset identifiers',
        "helpFile": 'import_merged#file_formats',
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
    contents_order=['file', 'crystalName', 'dataset', 'excludeSelection'],
    content_qualifiers={
        "file": {'allowUndefined': False, 'mustExist': True, 'fromPreviousJob': True},
        "crystalName": {'allowUndefined': True, 'minLength': 1, 'guiLabel': 'crystal name', 'allowedCharsCode': 1},
        "dataset": {'allowUndefined': True, 'minLength': 1, 'guiLabel': 'dataset name', 'allowedCharsCode': 1},
        "excludeSelection": {'allowUndefined': True},
    },
)
class CImportUnmergedStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CImportUnmerged.py
    to add methods and implementation-specific functionality.
    """

    file: Optional[CUnmergedDataFileStub] = None
    cell: Optional[CCellStub] = None
    wavelength: Optional[CWavelengthStub] = None
    crystalName: Optional[CString] = None
    dataset: Optional[CString] = None
    excludeSelection: Optional[CRangeSelectionStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CImportUnmergedStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "cell": attribute(AttributeType.CUSTOM, custom_class="CCellStub"),
        "spaceGroup": attribute(AttributeType.CUSTOM, custom_class="CSpaceGroupStub"),
        "resolutionRange": attribute(AttributeType.CUSTOM, custom_class="CResolutionRangeStub"),
        "listOfColumns": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "datasets": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "crystalNames": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "wavelengths": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "datasetCells": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "merged": attribute(AttributeType.BOOLEAN),
    },
    error_codes={
        "101": {
            "description": "Attempting to load MTZ data from non-existant/broken file"
        },
        "102": {
            "description": "Error creating command file for mtzdump"
        },
        "103": {
            "description": "No log file found from mtzdump"
        },
        "104": {
            "description": "Error reading log file from mtzdump"
        },
        "105": {
            "severity": 2,
            "description": "Different spacegroup"
        },
        "106": {
            "severity": 2,
            "description": "Different cell parameter"
        },
        "107": {
            "severity": 2,
            "description": "Different cell parameters"
        },
        "108": {
            "severity": 4,
            "description": "Different Laue group"
        },
        "109": {
            "severity": 4,
            "description": "Different point group"
        },
        "410": {
            "description": "Invalid CSeqDataFile passed to matthewCoeff"
        },
        "411": {
            "description": "Failed to run matthewCoeff"
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
class CMtzDataStub(CDataFileContent):
    """
    Base class for classes holding file contents

    This is a pure data class stub. Extend it in core/CMtzData.py
    to add methods and implementation-specific functionality.
    """

    cell: Optional[CCellStub] = None
    spaceGroup: Optional[CSpaceGroupStub] = None
    resolutionRange: Optional[CResolutionRangeStub] = None
    listOfColumns: Optional[CList] = None
    datasets: Optional[CList] = None
    crystalNames: Optional[CList] = None
    wavelengths: Optional[CList] = None
    datasetCells: Optional[CList] = None
    merged: Optional[CBoolean] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMtzDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "format": attribute(AttributeType.STRING),
        "merged": attribute(AttributeType.STRING),
        "crystalName": attribute(AttributeType.CUSTOM, custom_class="CCrystalNameStub"),
        "datasetName": attribute(AttributeType.CUSTOM, custom_class="CDatasetNameStub"),
        "cell": attribute(AttributeType.CUSTOM, custom_class="CCellStub"),
        "spaceGroup": attribute(AttributeType.CUSTOM, custom_class="CSpaceGroupStub"),
        "batchs": attribute(AttributeType.STRING),
        "lowRes": attribute(AttributeType.FLOAT),
        "highRes": attribute(AttributeType.FLOAT),
        "knowncell": attribute(AttributeType.BOOLEAN),
        "knownwavelength": attribute(AttributeType.BOOLEAN),
        "numberLattices": attribute(AttributeType.INT),
        "wavelength": attribute(AttributeType.CUSTOM, custom_class="CWavelengthStub"),
        "numberofdatasets": attribute(AttributeType.INT),
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
    content_qualifiers={
        "format": {'onlyEnumerators': True, 'enumerators': ['unk', 'mtz', 'xds', 'sca', 'saint', 'shelx', 'mmcif'], 'default': 'unk'},
        "merged": {'onlyEnumerators': True, 'enumerators': ['unk', 'merged', 'unmerged'], 'default': 'unk'},
    },
)
class CUnmergedDataContentStub(CDataFileContent):
    """
    Base class for classes holding file contents

    This is a pure data class stub. Extend it in core/CUnmergedDataContent.py
    to add methods and implementation-specific functionality.
    """

    format: Optional[CString] = None
    merged: Optional[CString] = None
    crystalName: Optional[CCrystalNameStub] = None
    datasetName: Optional[CDatasetNameStub] = None
    cell: Optional[CCellStub] = None
    spaceGroup: Optional[CSpaceGroupStub] = None
    batchs: Optional[CString] = None
    lowRes: Optional[CFloat] = None
    highRes: Optional[CFloat] = None
    knowncell: Optional[CBoolean] = None
    knownwavelength: Optional[CBoolean] = None
    numberLattices: Optional[CInt] = None
    wavelength: Optional[CWavelengthStub] = None
    numberofdatasets: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUnmergedDataContentStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "groupType": attribute(AttributeType.CUSTOM, custom_class="CMtzColumnGroupTypeStub"),
        "columns": attribute(AttributeType.CUSTOM, custom_class="CList"),
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
class CMtzColumnGroupStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CMtzColumnGroup.py
    to add methods and implementation-specific functionality.
    """

    groupType: Optional[CMtzColumnGroupTypeStub] = None
    columns: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMtzColumnGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "selected": attribute(AttributeType.BOOLEAN),
        "obsDataFile": attribute(AttributeType.CUSTOM, custom_class="CObsDataFileStub"),
        "crystalName": attribute(AttributeType.CUSTOM, custom_class="CCrystalNameStub"),
        "datasetName": attribute(AttributeType.CUSTOM, custom_class="CDatasetNameStub"),
        "formFactors": attribute(AttributeType.CUSTOM, custom_class="CFormFactorStub"),
        "formFactorSource": attribute(AttributeType.STRING),
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
    content_qualifiers={
        "formFactorSource": {'onlyEnumerators': True, 'enumerators': ['no', 'composition', 'xia2'], 'menuText': ['user input', 'atomic composition', 'from XIA2'], 'default': 'no'},
    },
)
class CDatasetStub(CData):
    """
    The experimental data model for ab initio phasing

    This is a pure data class stub. Extend it in core/CDataset.py
    to add methods and implementation-specific functionality.
    """

    selected: Optional[CBoolean] = None
    obsDataFile: Optional[CObsDataFileStub] = None
    crystalName: Optional[CCrystalNameStub] = None
    datasetName: Optional[CDatasetNameStub] = None
    formFactors: Optional[CFormFactorStub] = None
    formFactorSource: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDatasetStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "columnGroup": attribute(AttributeType.CUSTOM, custom_class="CMtzColumnGroupStub"),
        "datasetName": attribute(AttributeType.STRING),
    },
    error_codes={
        "101": {
            "description": "Column not in MTZ file"
        },
        "102": {
            "description": "Column wrong type"
        },
        "103": {
            "description": "MTZ file is not defined",
            "severity": 2
        },
        "104": {
            "description": "No column group selected"
        },
        "105": {
            "description": "No column group selected",
            "severity": 2
        }
    },
    qualifiers={
        "mustExist": False,
        "mtzFileKey": '',
        "groupTypes": [],
    },
    qualifiers_order=['groupTypes', 'mtzFileKey', 'mustExist'],
    qualifiers_definition={
        "groupTypes": {'type': 'list', 'description': 'Type of columnGroup required by program'},
        "mtzFileKey": {'type': 'str', 'description': 'The key for a CMtxDataFile in the same CContainer'},
        "mustExist": {'type': 'bool', 'description': 'Flag if the parameter must be set at run time'},
    },
)
class CProgramColumnGroup0Stub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CProgramColumnGroup0.py
    to add methods and implementation-specific functionality.
    """

    columnGroup: Optional[CMtzColumnGroupStub] = None
    datasetName: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProgramColumnGroup0Stub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
