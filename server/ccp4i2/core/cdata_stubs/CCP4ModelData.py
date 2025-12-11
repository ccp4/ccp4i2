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
from ccp4i2.core.cdata_stubs.CCP4Data import CDictStub, COneWordStub, CUUIDStub
from ccp4i2.core.cdata_stubs.CCP4File import CFilePathStub, CI2XmlDataFileStub, CProjectIdStub


@cdata_class(
    error_codes={
        "150": {
            "description": "No file content information"
        },
        "151": {
            "description": "Two sequences have the same identifier"
        },
        "152": {
            "description": "Failed in merging sequence files to read sequence file"
        },
        "153": {
            "description": "Failed in merging sequence files to write merged file"
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
class CSeqDataFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CSeqDataFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSeqDataFileListStub.

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
        "selection": attribute(AttributeType.CUSTOM, custom_class="CDictStub"),
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
        "mimeTypeName": 'application/CCP4-asu-content',
        "mimeTypeDescription": 'AU content',
        "fileExtensions": ['asu.xml'],
        "fileContentClassName": 'CAsuContent',
        "fileLabel": 'AU contents',
        "guiLabel": 'AU contents',
        "toolTip": 'A CCP4i2 file specifying AU contents',
        "helpFile": 'model_data#sequences',
        "saveToDb": True,
        "selectionMode": 0,
    },
    qualifiers_order=['autoLoadHeader'],
    qualifiers_definition={
        "selectionMode": {'type': 'int', 'description': 'Chain selection options'},
    },
    contents_order=['selection'],
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CAsuDataFileStub(CI2XmlDataFileStub):
    """
    A reference to an XML file with CCP4i2 Header

    This is a pure data class stub. Extend it in core/CAsuDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None
    selection: Optional[CDictStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuDataFileStub.

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
class CAtomRefmacSelectionListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CAtomRefmacSelectionList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomRefmacSelectionListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "groupId": attribute(AttributeType.INT),
        "chainIds": attribute(AttributeType.STRING),
        "firstRes": attribute(AttributeType.INT),
        "lastRes": attribute(AttributeType.INT),
        "atoms": attribute(AttributeType.STRING),
        "alt": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
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
    contents_order=[
        'groupId',
        'chainIds',
        'firstRes',
        'lastRes',
        'atoms',
        'alt'],
)
class CAtomRefmacSelectionOccupancyStub(CData):
    """
    A residue range selection for occupancy groups

    This is a pure data class stub. Extend it in core/CAtomRefmacSelectionOccupancy.py
    to add methods and implementation-specific functionality.
    """

    groupId: Optional[CInt] = None
    chainIds: Optional[CString] = None
    firstRes: Optional[CInt] = None
    lastRes: Optional[CInt] = None
    atoms: Optional[CString] = None
    alt: Optional[COneWordStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomRefmacSelectionOccupancyStub.

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
class CResidueRangeListStub(CList):
    """
    A list of residue range selections

    This is a pure data class stub. Extend it in core/CResidueRangeList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CResidueRangeListStub.

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
        "fileLabel": 'tls',
        "mimeTypeName": 'application/refmac-TLS',
        "mimeTypeDescription": 'Refmac TLS file',
        "guiLabel": 'TLS coefficients',
        "toolTip": 'Definition of model domains for TLS refinement',
        "fileExtensions": ['tls'],
        "fileContentClassName": None,
        "helpFile": 'model_data#tls_file',
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
class CTLSDataFileStub(CDataFile):
    """
    A refmac TLS file

    This is a pure data class stub. Extend it in core/CTLSDataFile.py
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
        Initialize CTLSDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "401": {
            "description": "Non-alphabet character removed from sequence",
            "severity": 2
        },
        "402": {
            "description": "Invalid characters (BJOXZ) in sequence"
        },
        "403": {
            "description": "Sequence undefined",
            "severity": 2
        }
    },
    qualifiers={
        "minLength": None,
        "maxLength": None,
        "enumerators": [],
        "menuText": [],
        "onlyEnumerators": False,
        "charWidth": -1,
        "allowedCharsCode": 0,
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
class CSequenceStringStub(CString):
    """
    A string

    This is a pure data class stub. Extend it in core/CSequenceString.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSequenceStringStub.

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
        "fileLabel": 'HHPred sequence search',
        "mimeTypeName": 'application/HHPred-alignments',
        "mimeTypeDescription": 'HHPred sequence search results',
        "guiLabel": 'HHPred results',
        "tooltip": 'Output from HHPred search',
        "fileExtensions": ['hhr'],
        "fileContentClassName": 'CHhpredData',
        "helpFile": 'model_data#ali',
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
class CHhpredDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CHhpredDataFile.py
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
        Initialize CHhpredDataFileStub.

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
        "onlyEnumerators": True,
        "enumerators": ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'],
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
class CElementStub(COneWordStub):
    """
    Chemical element

    This is a pure data class stub. Extend it in core/CElement.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CElementStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "identifier": attribute(AttributeType.STRING),
        "referenceDb": attribute(AttributeType.STRING),
        "reference": attribute(AttributeType.STRING),
        "name": attribute(AttributeType.STRING),
        "description": attribute(AttributeType.STRING),
        "sequence": attribute(AttributeType.STRING),
        "moleculeType": attribute(AttributeType.STRING),
    },
    error_codes={
        "201": {
            "description": "Sequence undefined",
            "severity": 1
        },
        "202": {
            "description": "error reading from file"
        },
        "203": {
            "description": "Comparing sequences: Sequence item different"
        },
        "204": {
            "description": "Comparing sequences: One item set - the other is unset"
        },
        "401": {
            "description": "Attempting to load from non-existent file"
        },
        "402": {
            "description": "Error reading from file"
        },
        "403": {
            "description": "Unknown sequence file format"
        },
        "405": {
            "description": "Error reading identifiers from multi-record file"
        },
        "406": {
            "description": "Error opening file"
        },
        "407": {
            "description": "The 'PIR' file did not have the correct format"
        },
        "408": {
            "severity": 2,
            "description": "The 'PIR' file format was corrected"
        },
        "409": {
            "description": "Error opening file to write"
        },
        "410": {
            "description": "Error attempting to write out sequence file"
        },
        "411": {
            "description": "Error attempting to create a temporary sequence file"
        },
        "412": {
            "description": "Sequence file is empty"
        },
        "413": {
            "description": "Unable to read BLAST format file"
        },
        "414": {
            "description": "Unable to read hhpred format file"
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
    contents_order=[
        'identifier',
        'name',
        'description',
        'referenceDb',
        'reference',
        'moleculeType',
        'sequence'],
    content_qualifiers={
        "identifier": {'toolTip': 'Description of sequence', 'minlength': 4},
        "referenceDb": {'onlyEnumerators': False, 'default': 'unk', 'enumerators': ['unk', 'sp', 'tr', 'pdb'], 'menuText': ['Unknown', 'UniProt/Swiss-Prot', 'UniProt/TrEMBL', 'ProteinDatabank']},
        "reference": {'toolTip': 'Optional reference for sequence'},
        "name": {'toolTip': 'User friendly name of sequence'},
        "description": {'toolTip': 'User friendly description of sequence'},
        "sequence": {'toolTip': 'Single letter sequence (white space and dash ignored)'},
        "moleculeType": {'onlyEnumerators': True, 'enumerators': ['PROTEIN', 'NUCLEIC'], 'menuText': ['protein', 'nucleic acid'], 'default': 'PROTEIN', 'toolTip': 'Molecule type'},
    },
)
class CSequenceStub(CData):
    """
    A string of sequence one-letter codes
Need to be able to parse common seq file formats
Do we need to support alternative residues
What about nucleic/polysach?

    This is a pure data class stub. Extend it in core/CSequence.py
    to add methods and implementation-specific functionality.
    """

    identifier: Optional[CString] = None
    referenceDb: Optional[CString] = None
    reference: Optional[CString] = None
    name: Optional[CString] = None
    description: Optional[CString] = None
    sequence: Optional[CString] = None
    moleculeType: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSequenceStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "identifier": attribute(AttributeType.STRING),
        "moleculeType": attribute(AttributeType.STRING),
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
    contents_order=['identifier', 'moleculeType'],
    content_qualifiers={
        "identifier": {'toolTip': 'Optional convenient name for sequence alignment'},
        "moleculeType": {'onlyEnumerators': True, 'enumerators': ['PROTEIN', 'NUCLEIC'], 'menuText': ['protein', 'nucleic acid'], 'default': 'PROTEIN', 'toolTip': 'Molecule type'},
    },
)
class CSequenceAlignmentStub(CData):
    """
    An alignment of two or more sequences.
Each sequence is obviously related to class CSequence, but
will also contain gaps relevant to the alignment. We could
implement the contents as a list of CSequence objects?
The alignment is typically formatted in a file as consecutive
or interleaved sequences.

    This is a pure data class stub. Extend it in core/CSequenceAlignment.py
    to add methods and implementation-specific functionality.
    """

    identifier: Optional[CString] = None
    moleculeType: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSequenceAlignmentStub.

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
        "fileLabel": 'Blast sequence search',
        "mimeTypeName": 'application/Blast-alignments',
        "mimeTypeDescription": 'Blast sequence search results',
        "guiLabel": 'Blast results',
        "tooltip": 'Output from Blast search',
        "fileExtensions": ['bla', 'blast', 'xml'],
        "fileContentClassName": 'CBlastData',
        "helpFile": 'model_data#ali',
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
class CBlastDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CBlastDataFile.py
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
        Initialize CBlastDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "monomerList": attribute(AttributeType.CUSTOM, custom_class="CList"),
    },
    error_codes={
        "101": {
            "description": "Error opening MMCIF format file"
        },
        "102": {
            "description": "Error merging data - monomer already in geometry file"
        },
        "103": {
            "severity": 2,
            "description": "Warning merging data - overwriting geometry for monomer with same id"
        },
        "104": {
            "description": "Error reading geometry cif file - does not contain expected data"
        },
        "105": {
            "description": "Unknown error reading geometry file"
        },
        "106": {
            "description": "_chem_comp section not found in geometry file"
        },
        "110": {
            "description": "Attemting to delete unrecognised chem_comp.id"
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
class CDictDataStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CDictData.py
    to add methods and implementation-specific functionality.
    """

    monomerList: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDictDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "identifier": attribute(AttributeType.STRING),
        "formula": attribute(AttributeType.STRING),
        "dictionaryName": attribute(AttributeType.STRING),
        "smiles": attribute(AttributeType.STRING),
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
    contents_order=['identifier', 'formula', 'dictionaryName', 'smiles'],
    content_qualifiers={
        "identifier": {'toolTip': 'The name you use for the monomer'},
        "formula": {'toolTip': 'The formula for the monomer'},
        "dictionaryName": {'toolTip': 'The REFMAC dictionary name if not the same as the name'},
        "smiles": {'toolTip': 'The smiles string for the monomer'},
    },
)
class CMonomerStub(CData):
    """
    A monomer compound. ?smiles

    This is a pure data class stub. Extend it in core/CMonomer.py
    to add methods and implementation-specific functionality.
    """

    identifier: Optional[CString] = None
    formula: Optional[CString] = None
    dictionaryName: Optional[CString] = None
    smiles: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMonomerStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "hitId": attribute(AttributeType.STRING),
        "querySequence": attribute(AttributeType.STRING),
        "hitSequence": attribute(AttributeType.STRING),
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
class CBlastItemStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CBlastItem.py
    to add methods and implementation-specific functionality.
    """

    hitId: Optional[CString] = None
    querySequence: Optional[CString] = None
    hitSequence: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBlastItemStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "groupIds": attribute(AttributeType.STRING),
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
    contents_order=['groupIds'],
)
class CAtomRefmacSelectionGroupsStub(CData):
    """
    A group selection for occupancy groups

    This is a pure data class stub. Extend it in core/CAtomRefmacSelectionGroups.py
    to add methods and implementation-specific functionality.
    """

    groupIds: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomRefmacSelectionGroupsStub.

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
class COccRelationRefmacListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/COccRelationRefmacList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize COccRelationRefmacListStub.

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
            "description": "Error attempting to merge geometry files - no libcheck script"
        },
        "202": {
            "description": "Error attempting to merge geometry files - failed creating working directory"
        },
        "203": {
            "description": "Error attempting to merge geometry files - setting libcheck parameters"
        },
        "204": {
            "description": "Error attempting to merge geometry files - running libcheck"
        },
        "205": {
            "description": "Error attempting to merge geometry files - failed to run libcheck"
        }
    },
    qualifiers={
        "fileLabel": 'dictionary',
        "mimeTypeName": 'application/refmac-dictionary',
        "mimeTypeDescription": 'Geometry file',
        "guiLabel": 'Geometry dictionary',
        "toolTip": 'Idealised geometry of ligands for refinement',
        "fileExtensions": ['cif'],
        "fileContentClassName": 'CDictData',
        "helpFile": 'model_data#ligand_geometry',
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
class CDictDataFileStub(CDataFile):
    """
    A refmac dictionary file

    This is a pure data class stub. Extend it in core/CDictDataFile.py
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
        Initialize CDictDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "label": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "number": attribute(AttributeType.INT),
        "use": attribute(AttributeType.BOOLEAN),
        "pdbItemList": attribute(AttributeType.CUSTOM, custom_class="CList"),
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
        "guiLabel": 'Ensemble',
        "allowUndefined": False,
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
        "number": {'min': 0, 'default': 1, 'enumerators': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'menuText': ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']},
        "use": {'default': True},
        "pdbItemList": {'listMinLength': 1},
    },
)
class CEnsembleStub(CData):
    """
    An ensemble of models. Typically, this would be a set of related
PDB files, but models could also be xtal or EM maps. This should
be indicated by the types entry.
A single ensemble is a CList of structures.

    This is a pure data class stub. Extend it in core/CEnsemble.py
    to add methods and implementation-specific functionality.
    """

    label: Optional[COneWordStub] = None
    number: Optional[CInt] = None
    use: Optional[CBoolean] = None
    pdbItemList: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEnsembleStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "chainId": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "firstRes": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "lastRes": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
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
        "pdbFileKey": None,
    },
    qualifiers_order=['pdbFileKey'],
    qualifiers_definition={
        "pdbFileKey": {'type': 'str', 'description': 'The key for a CPdbDataFile in the same CContainer'},
    },
    contents_order=['chainId', 'firstRes', 'lastRes'],
    content_qualifiers={
        "chainId": {'default': ''},
    },
)
class CResidueRangeStub(CData):
    """
    A residue range selection

    This is a pure data class stub. Extend it in core/CResidueRange.py
    to add methods and implementation-specific functionality.
    """

    chainId: Optional[COneWordStub] = None
    firstRes: Optional[COneWordStub] = None
    lastRes: Optional[COneWordStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CResidueRangeStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "text": attribute(AttributeType.STRING),
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
        "pdbFileKey": '',
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
        "pdbFileKey": {'type': 'str', 'description': 'The key for a CPdbDataFile in the same CContainer'},
    },
)
class CAtomSelectionStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CAtomSelection.py
    to add methods and implementation-specific functionality.
    """

    text: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomSelectionStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "queryId": attribute(AttributeType.STRING),
        "alignmentList": attribute(AttributeType.CUSTOM, custom_class="CList"),
    },
    error_codes={
        "201": {
            "description": "Failed reading blast file"
        },
        "202": {
            "description": "Blast file contains results of more than one query - only the first is read",
            "severity": 2
        },
        "203": {
            "description": "Failed parsing Blast file"
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
class CBlastDataStub(CDataFileContent):
    """
    Base class for classes holding file contents

    This is a pure data class stub. Extend it in core/CBlastData.py
    to add methods and implementation-specific functionality.
    """

    queryId: Optional[CString] = None
    alignmentList: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBlastDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "alignmentList": attribute(AttributeType.CUSTOM, custom_class="CList"),
    },
    error_codes={
        "201": {
            "description": "Failed to read HHPred file"
        },
        "202": {
            "description": "Failed to load iotbx software to read HHPred file"
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
class CHhpredDataStub(CDataFileContent):
    """
    Base class for classes holding file contents

    This is a pure data class stub. Extend it in core/CHhpredData.py
    to add methods and implementation-specific functionality.
    """

    alignmentList: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHhpredDataStub.

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
class CPdbDataFileListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CPdbDataFileList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPdbDataFileListStub.

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
        "fileLabel": 'mol2',
        "mimeTypeName": 'chemical/x-mol2',
        "mimeTypeDescription": 'MOL2 file',
        "guiLabel": 'MOL2 file',
        "toolTip": 'Structure geometry of ligands for refinement in MOL2 format',
        "fileExtensions": ['mol2'],
        "fileContentClassName": None,
        "helpFile": 'model_data#mol2_file',
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
class CMol2DataFileStub(CDataFile):
    """
    A molecule definition file (MOL2)

    This is a pure data class stub. Extend it in core/CMol2DataFile.py
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
        Initialize CMol2DataFileStub.

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
class COccRefmacSelectionListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/COccRefmacSelectionList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize COccRefmacSelectionListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "uniprotId": attribute(AttributeType.STRING),
        "organism": attribute(AttributeType.STRING),
        "expressionSystem": attribute(AttributeType.STRING),
    },
    error_codes={
        "401": {
            "description": "No uniprot id available"
        },
        "402": {
            "description": "No uniprot xml file available to read"
        },
        "403": {
            "description": "No project id provided to determine uniprot xml filename"
        },
        "404": {
            "description": "Reading uniprot xml file failed"
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
class CSequenceMetaStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CSequenceMeta.py
    to add methods and implementation-specific functionality.
    """

    uniprotId: Optional[CString] = None
    organism: Optional[CString] = None
    expressionSystem: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSequenceMetaStub.

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
            "description": "Error reading sequence file"
        },
        "202": {
            "description": "Error in BioPython attempting to identify file type"
        }
    },
    qualifiers={
        "fileLabel": 'sequence',
        "mimeTypeName": 'application/CCP4-seq',
        "mimeTypeDescription": 'Sequence file',
        "guiLabel": 'Sequence',
        "tooltip": 'Sequence in any of the common formats (pir,fasta..)',
        "fileExtensions": ['seq', 'pir', 'fasta'],
        "fileContentClassName": 'CSequence',
        "downloadModes": ['uniprotFasta'],
        "helpFile": 'model_data#sequences',
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
class CSeqDataFileStub(CDataFile):
    """
    A sequence file

    This is a pure data class stub. Extend it in core/CSeqDataFile.py
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
        Initialize CSeqDataFileStub.

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
        "202": {
            "description": "Error reading from file"
        },
        "203": {
            "description": "Unknown alignment file format"
        },
        "204": {
            "description": "Can not read Blast or HHPred file format"
        },
        "205": {
            "description": "Error reading identifiers from multi-record file"
        },
        "206": {
            "description": "Error attempting to identify file format"
        },
        "250": {
            "description": "Alignment file format not recognised - can not convert"
        },
        "251": {
            "description": "Alignment file conversion failed to overwrite existing file"
        },
        "252": {
            "description": "Alignment file conversion failed writing file"
        },
        "260": {
            "description": "Alignment file does not contain required number of sequences"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "mustExist": False,
        "fromPreviousJob": False,
        "jobCombo": True,
        "mimeTypeName": 'application/CCP4-seqalign',
        "mimeTypeDescription": 'Sequence alignment file',
        "fileLabel": None,
        "fileExtensions": ['aln', 'pir', 'fasta', 'msf', 'phy'],
        "fileContentClassName": 'CSequenceAlignment',
        "isDirectory": False,
        "saveToDb": True,
        "requiredSubType": None,
        "requiredContentFlag": None,
        "guiLabel": 'Aligned sequences',
        "toolTip": 'Multiple sequence alignment in any of the common formats (pir,fasta..)',
        "helpFile": 'model_data#alignments',
    },
    qualifiers_order=['requiredSequences'],
    qualifiers_definition={
        "requiredSequences": {'type': 'list', 'listItemType': "<class 'int'>", 'description': 'A list of allowed numbers of sequences in file (usually [2])'},
    },
    content_qualifiers={
        "subType": {'default': None},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CSeqAlignDataFileStub(CDataFile):
    """
    A (multiple) sequence alignment file

    This is a pure data class stub. Extend it in core/CSeqAlignDataFile.py
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
        Initialize CSeqAlignDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "annotation": attribute(AttributeType.STRING),
        "identifier": attribute(AttributeType.STRING),
        "chain": attribute(AttributeType.STRING),
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
class CHhpredItemStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CHhpredItem.py
    to add methods and implementation-specific functionality.
    """

    annotation: Optional[CString] = None
    identifier: Optional[CString] = None
    chain: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHhpredItemStub.

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
class CEnsembleListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CEnsembleList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEnsembleListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "101": {
            "description": "Unable to load mmdb - ensure LD_LIBRARY_PATH is set"
        },
        "102": {
            "description": "Error reading PDB file into MMDB object"
        },
        "103": {
            "description": "Residue range selection does not specify chain"
        },
        "104": {
            "description": "Residue range selection specifies non-existant chain id"
        },
        "105": {
            "description": "Residue range selection - no residues selected"
        },
        "106": {
            "description": "Residue range selection - residue number is not an integer"
        },
        "112": {
            "description": "Atom selection failed. Failed creating CMMDBManager object"
        },
        "113": {
            "description": "Atom selection failed. Faied reading coordinate file."
        },
        "114": {
            "description": "Atom selection failed. Failed parsing command"
        },
        "115": {
            "description": "Atom selection failed. Error creating PPCAtom"
        },
        "116": {
            "description": "Atom selection failed. Error in GetSelIndex"
        },
        "117": {
            "description": "Atom selection failed. Error loading selection tree"
        },
        "118": {
            "description": "Atom selection failed. Error applying selection tree"
        },
        "119": {
            "description": "Creating new PDB file failed on writing file"
        },
        "120": {
            "description": "Creating new PDB file failed converting from fractional coordinates"
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
class CPdbDataStub(CDataFileContent):
    """
    Contents of a PDB file - a subset with functionality for GUI

    This is a pure data class stub. Extend it in core/CPdbData.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPdbDataStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "groupId": attribute(AttributeType.INT),
        "chainId": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "firstRes": attribute(AttributeType.INT),
        "lastRes": attribute(AttributeType.INT),
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
    contents_order=['groupId', 'chainId', 'firstRes', 'lastRes'],
)
class CAtomRefmacSelectionStub(CData):
    """
    A residue range selection for rigid body groups

    This is a pure data class stub. Extend it in core/CAtomRefmacSelection.py
    to add methods and implementation-specific functionality.
    """

    groupId: Optional[CInt] = None
    chainId: Optional[COneWordStub] = None
    firstRes: Optional[CInt] = None
    lastRes: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomRefmacSelectionStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "id": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "three_letter_code": attribute(AttributeType.CUSTOM, custom_class="COneWordStub"),
        "name": attribute(AttributeType.STRING),
        "group": attribute(AttributeType.STRING),
        "number_atoms_all": attribute(AttributeType.INT),
        "number_atoms_nh": attribute(AttributeType.INT),
        "desc_level": attribute(AttributeType.INT),
    },
    error_codes={
        "201": {
            "description": "Error reading monomer id and name"
        },
        "202": {
            "description": "Error writing monomer id and name"
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
class CChemCompStub(CData):
    """
    Component of CDictDataFile contents

    This is a pure data class stub. Extend it in core/CChemComp.py
    to add methods and implementation-specific functionality.
    """

    id: Optional[COneWordStub] = None
    three_letter_code: Optional[COneWordStub] = None
    name: Optional[CString] = None
    group: Optional[CString] = None
    number_atoms_all: Optional[CInt] = None
    number_atoms_nh: Optional[CInt] = None
    desc_level: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CChemCompStub.

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
        "fileLabel": 'mol',
        "mimeTypeName": 'chemical/x-mdl-molfile',
        "mimeTypeDescription": 'MDL Molfile',
        "guiLabel": 'Mol file',
        "toolTip": 'Structure geometry of ligands for refinement in MDL mol format',
        "fileExtensions": ['mol', 'sdf'],
        "fileContentClassName": None,
        "helpFile": 'model_data#mol_file',
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
class CMDLMolDataFileStub(CDataFile):
    """
    A molecule definition file (MDL)

    This is a pure data class stub. Extend it in core/CMDLMolDataFile.py
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
        Initialize CMDLMolDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    error_codes={
        "401": {
            "description": "Sequence the same as a sequence that is already loaded"
        },
        "402": {
            "description": "Sequence names are not unique: "
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
class CAsuContentSeqListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CAsuContentSeqList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuContentSeqListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "sequence": attribute(AttributeType.CUSTOM, custom_class="CSequenceStringStub"),
        "nCopies": attribute(AttributeType.INT),
        "polymerType": attribute(AttributeType.STRING),
        "name": attribute(AttributeType.STRING),
        "description": attribute(AttributeType.STRING),
        "source": attribute(AttributeType.CUSTOM, custom_class="CSeqDataFileStub"),
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
        "sequence": {'allowUndefined': False, 'minLength': 1},
        "nCopies": {'enumerators': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'default': 1, 'min': 0},
        "polymerType": {'onlyEnumerators': True, 'enumerators': ['PROTEIN', 'RNA', 'DNA'], 'default': 'PROTEIN'},
        "name": {'allowUndefined': False, 'minLength': 1, 'allowedCharsCode': 1},
        "description": {'allowUndefined': True},
    },
)
class CAsuContentSeqStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CAsuContentSeq.py
    to add methods and implementation-specific functionality.
    """

    sequence: Optional[CSequenceStringStub] = None
    nCopies: Optional[CInt] = None
    polymerType: Optional[CString] = None
    name: Optional[CString] = None
    description: Optional[CString] = None
    source: Optional[CSeqDataFileStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuContentSeqStub.

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
        "selection": attribute(AttributeType.CUSTOM, custom_class="CAtomSelectionStub"),
    },
    error_codes={
        "401": {
            "description": "Failed running coord_format to fix coordinate file - is it a PDB file?"
        },
        "402": {
            "severity": 2,
            "description": "Badly formated PDB file fixed"
        },
        "403": {
            "severity": 2,
            "description": "Fixed by removing text"
        },
        "404": {
            "severity": 2,
            "description": "Fixed by adding text"
        },
        "405": {
            "description": "There are no ATOM or HETATM lines in the PDB file"
        },
        "410": {
            "description": "No file loaded - can not convert coordinate file format"
        },
        "411": {
            "description": "Failed loading file - can not convert coordinate file format"
        },
        "412": {
            "description": "Can not overwrite existing file - can not convert coordinate file format"
        },
        "413": {
            "description": "Failed writing coordinate file"
        },
        "414": {
            "description": "Failed to identify coordinate file format"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "mustExist": False,
        "fromPreviousJob": False,
        "jobCombo": True,
        "mimeTypeName": 'chemical/x-pdb',
        "mimeTypeDescription": 'Model coordinates',
        "fileLabel": 'coordinates',
        "fileExtensions": ['pdb', 'cif', 'mmcif', 'ent'],
        "fileContentClassName": 'CPdbData',
        "isDirectory": False,
        "saveToDb": True,
        "requiredSubType": None,
        "requiredContentFlag": None,
        "guiLabel": 'Atomic model',
        "toolTip": 'A model coordinate file in PDB or mmCIF format',
        "ifInfo": True,
        "ifAtomSelection": False,
        "downloadModes": ['ebiPdb', 'rcsbPdb', 'uniprotAFPdb'],
        "helpFile": 'model_data#coordinate_files',
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
        "ifAtomSelection": {'type': 'bool', 'description': 'Atom selection option enabled'},
    },
    content_qualifiers={
        "subType": {'default': 0, 'enumerators': [0, 1, 2, 3, 4], 'onlyEnumerators': True, 'menuText': ['unknown', 'model', 'homolog', 'fragment', 'heavy atoms']},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CPdbDataFileStub(CDataFile):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CPdbDataFile.py
    to add methods and implementation-specific functionality.
    """

    # Subtype constants
    SUBTYPE_UNKNOWN = 0  # unknown
    SUBTYPE_MODEL = 1  # model
    SUBTYPE_HOMOLOG = 2  # homolog
    SUBTYPE_FRAGMENT = 3  # fragment
    SUBTYPE_HEAVY_ATOMS = 4  # heavy atoms

    # Content flag constants
    CONTENT_FLAG_PDB = 1  # PDB format
    CONTENT_FLAG_MMCIF = 2  # mmCIF format

    # Content annotations (indexed by contentFlag - 1)
    CONTENT_ANNOTATION = ['PDB format', 'mmCIF format']

    # Column signatures for each content flag (indexed by contentFlag - 1)
    CONTENT_SIGNATURE_LIST = [None, None]

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None
    selection: Optional[CAtomSelectionStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPdbDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "seqList": attribute(AttributeType.CUSTOM, custom_class="CAsuContentSeqListStub"),
    },
    error_codes={
        "101": {
            "description": "Failed reading file - is it correct file type?"
        },
        "102": {
            "description": "Failed reading file - it is not AU contents file"
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
class CAsuContentStub(CDataFileContent):
    """
    Base class for classes holding file contents

    This is a pure data class stub. Extend it in core/CAsuContent.py
    to add methods and implementation-specific functionality.
    """

    seqList: Optional[CAsuContentSeqListStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAsuContentStub.

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
        "selection": attribute(AttributeType.CUSTOM, custom_class="CAtomSelectionStub"),
    },
    error_codes={
        "401": {
            "description": "Failed running coord_format to fix coordinate file - is it a PDB file?"
        },
        "402": {
            "severity": 2,
            "description": "Badly formated PDB file fixed"
        },
        "403": {
            "severity": 2,
            "description": "Fixed by removing text"
        },
        "404": {
            "severity": 2,
            "description": "Fixed by adding text"
        },
        "405": {
            "description": "There are no ATOM or HETATM lines in the PDB file"
        },
        "410": {
            "description": "No file loaded - can not convert coordinate file format"
        },
        "411": {
            "description": "Failed loading file - can not convert coordinate file format"
        },
        "412": {
            "description": "Can not overwrite existing file - can not convert coordinate file format"
        },
        "413": {
            "description": "Failed writing coordinate file"
        },
        "414": {
            "description": "Failed to identify coordinate file format"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "mustExist": False,
        "fromPreviousJob": False,
        "jobCombo": True,
        "mimeTypeName": 'chemical/x-pdb',
        "mimeTypeDescription": 'Model coordinates',
        "fileLabel": 'ensemble coordinates',
        "fileExtensions": ['pdb', 'cif', 'mmcif', 'ent'],
        "fileContentClassName": 'CPdbData',
        "isDirectory": False,
        "saveToDb": True,
        "requiredSubType": None,
        "requiredContentFlag": None,
        "guiLabel": 'Model ensemble',
        "toolTip": 'An ensemble of model coordinates in PDB or mmCIF format',
        "ifInfo": True,
        "ifAtomSelection": False,
        "downloadModes": [],
        "helpFile": 'model_data#ensemble_coordinate_files',
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
        "ifAtomSelection": {'type': 'bool', 'description': 'Atom selection option enabled'},
    },
    content_qualifiers={
        "subType": {'default': 0, 'enumerators': [0, 1, 2, 3, 4], 'onlyEnumerators': True, 'menuText': ['unknown', 'model', 'homolog', 'fragment', 'heavy atoms']},
        "contentFlag": {'min': 0, 'default': None},
    },
)
class CEnsemblePdbDataFileStub(CPdbDataFileStub):
    """
    A PDB coordinate file containing ensemble of structures as 'NMR' models

    This is a pure data class stub. Extend it in core/CEnsemblePdbDataFile.py
    to add methods and implementation-specific functionality.
    """

    project: Optional[CProjectIdStub] = None
    baseName: Optional[CFilePathStub] = None
    relPath: Optional[CFilePathStub] = None
    annotation: Optional[CString] = None
    dbFileId: Optional[CUUIDStub] = None
    subType: Optional[CInt] = None
    contentFlag: Optional[CInt] = None
    selection: Optional[CAtomSelectionStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEnsemblePdbDataFileStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "structure": attribute(AttributeType.CUSTOM, custom_class="CPdbDataFileStub"),
        "identity_to_target": attribute(AttributeType.FLOAT),
        "rms_to_target": attribute(AttributeType.FLOAT),
    },
    error_codes={
        "101": {
            "description": "No sequence identity or structure RMS to target set"
        }
    },
    qualifiers={
        "guiLabel": 'Structure in ensemble',
        "toolTip": 'Homologous model and its similarity to the target structure',
        "allowUndefined": False,
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
    contents_order=['structure', 'identity_to_target', 'rms_to_target'],
    content_qualifiers={
        "structure": {'allowUndefined': False, 'mustExist': True, 'fromPreviousJob': True, 'ifAtomSelection': True},
        "identity_to_target": {'min': 0.0, 'max': 1.0},
        "rms_to_target": {'min': 0.0, 'max': 100.0},
    },
)
class CPdbEnsembleItemStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CPdbEnsembleItem.py
    to add methods and implementation-specific functionality.
    """

    structure: Optional[CPdbDataFileStub] = None
    identity_to_target: Optional[CFloat] = None
    rms_to_target: Optional[CFloat] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPdbEnsembleItemStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
