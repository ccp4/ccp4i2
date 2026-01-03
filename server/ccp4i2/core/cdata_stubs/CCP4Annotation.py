"""Auto-generated from CCP4i2 metadata. DO NOT EDIT.

This is a stub file - extend classes in core/ to add methods.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Optional, Any

# Metadata system
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType

# Base classes
from ccp4i2.core.base_object.base_classes import CData, CDataFile

# Fundamental types
from ccp4i2.core.base_object.fundamental_types import CBoolean, CFloat, CInt, CList, CString


@cdata_class(
    attributes={
        "pmid": attribute(AttributeType.INT),
        "title": attribute(AttributeType.STRING),
        "authorList": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "source": attribute(AttributeType.STRING),
        "url": attribute(AttributeType.STRING),
        "selected": attribute(AttributeType.BOOLEAN),
    },
    error_codes={
        "101": {
            "description": "Failed to load Medline data"
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
        'pmid',
        'title',
        'authorList',
        'source',
        'url',
        'selected'],
)
class CBibReferenceStub(CData):
    """
    Bibliographic reference

    This is a pure data class stub. Extend it in core/CBibReference.py
    to add methods and implementation-specific functionality.
    """

    pmid: Optional[CInt] = None
    title: Optional[CString] = None
    authorList: Optional[CList] = None
    source: Optional[CString] = None
    url: Optional[CString] = None
    selected: Optional[CBoolean] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReferenceStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "taskName": attribute(AttributeType.STRING),
        "version": attribute(AttributeType.STRING),
        "title": attribute(AttributeType.STRING),
        "references": attribute(AttributeType.CUSTOM, custom_class="CList"),
    },
    error_codes={
        "100": {
            "description": "Failed attempting to load MedLine file - file not found"
        },
        "101": {
            "description": "Failed attempting to find references file"
        },
        "102": {
            "description": "Error copying file"
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
    contents_order=['taskName', 'version', 'title', 'references'],
)
class CBibReferenceGroupStub(CData):
    """
    Set of bibliographic references for a task

    This is a pure data class stub. Extend it in core/CBibReferenceGroup.py
    to add methods and implementation-specific functionality.
    """

    taskName: Optional[CString] = None
    version: Optional[CString] = None
    title: Optional[CString] = None
    references: Optional[CList] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReferenceGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "family": attribute(AttributeType.STRING),
        "style": attribute(AttributeType.INT),
        "pointSize": attribute(AttributeType.INT),
        "weight": attribute(AttributeType.INT),
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
        "family": {'default': 'Helvetica'},
        "style": {'onlyEnumerators': True, 'default': 0, 'enumerators': [0, 1, 2], 'menuText': ['normal', 'italic', 'oblique']},
        "pointSize": {'min': 1, 'default': 12},
        "weight": {'min': 0, 'max': 99, 'default': 50, 'allowUndefined': False, 'enumerators': [25, 50, 63, 75, 87], 'menuText': ['light', 'normal', 'demi-bold', 'bold', 'black']},
    },
)
class CFontStub(CData):
    """
    Simplified Qt font options

    This is a pure data class stub. Extend it in core/CFont.py
    to add methods and implementation-specific functionality.
    """

    family: Optional[CString] = None
    style: Optional[CInt] = None
    pointSize: Optional[CInt] = None
    weight: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFontStub.

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
        "label": 'User id',
        "toolTip": 'User id as me@myplace.ac.uk',
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'],
)
class CUserIdStub(CString):
    """
    A user ID

    This is a pure data class stub. Extend it in core/CUserId.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUserIdStub.

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
        "min": 0,
        "label": 'Time',
        "toolTip": 'Time and date as hh:mm dd/mm/yyyy',
        "format": '%H:%M %d/%b/%y',
    },
    qualifiers_order=['format'],
    qualifiers_definition={
        "format": {'type': 'str', 'description': 'Argument to Python time.strftime to display time in human readable format'},
    },
)
class CTimeStub(CInt):
    """
    The time. Uses Python time module

    This is a pure data class stub. Extend it in core/CTime.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTimeStub.

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
class CAnnotationListStub(CList):
    """
    A list of annotation

    This is a pure data class stub. Extend it in core/CAnnotationList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnnotationListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "tag": attribute(AttributeType.STRING),
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
        "enumeratorsFunction": None,
        "addEnumeratorFunction": None,
    },
    qualifiers_order=['enumeratorsFunction', 'addEnumeratorFunction'],
    qualifiers_definition={
        "enumeratorsFunction": {'type': '"method"', 'definition': 'Function returning list of enumerators'},
        "addEnumeratorFunction": {'type': '"method"', 'definition': 'Function to add to list of enumerators'},
    },
    contents_order=['tag'],
)
class CMetaDataTagStub(CData):
    """
    This class will extend list of enumerators if new value for string is entered

    This is a pure data class stub. Extend it in core/CMetaDataTag.py
    to add methods and implementation-specific functionality.
    """

    tag: Optional[CString] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMetaDataTagStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "year": attribute(AttributeType.INT),
        "month": attribute(AttributeType.STRING),
        "day": attribute(AttributeType.INT),
        "yearRange": attribute(AttributeType.INT),
        "monthRange": attribute(AttributeType.INT),
        "dayRange": attribute(AttributeType.INT),
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
        'year',
        'month',
        'day',
        'yearRange',
        'monthRange',
        'dayRange'],
    content_qualifiers={
        "year": {'enumerators': []},
        "month": {'onlyEnumerators': True, 'enumerators': ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'], 'default': 'January'},
        "day": {'default': 1, 'min': 1, 'max': 31},
        "yearRange": {'default': 0, 'min': 0, 'max': 100},
        "monthRange": {'default': 0, 'min': 0, 'max': 12},
        "dayRange": {'default': 0, 'min': 0, 'max': 30},
    },
)
class CDateRangeStub(CData):
    """
    A date range - may be on a scale of years,months or days

    This is a pure data class stub. Extend it in core/CDateRange.py
    to add methods and implementation-specific functionality.
    """

    year: Optional[CInt] = None
    month: Optional[CString] = None
    day: Optional[CInt] = None
    yearRange: Optional[CInt] = None
    monthRange: Optional[CInt] = None
    dayRange: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDateRangeStub.

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
class CAuthorStub(CString):
    """
    Placeholder for bibliographic author

    This is a pure data class stub. Extend it in core/CAuthor.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAuthorStub.

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
        "label": 'Machine name',
        "toolTip": 'Hostname as mycomputer.myplace.ac.uk',
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
class CHostNameStub(CString):
    """
    Computer name

    This is a pure data class stub. Extend it in core/CHostName.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHostNameStub.

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
class CMetaDataTagListStub(CList):
    """
    A list with all items of one CData sub-class

    This is a pure data class stub. Extend it in core/CMetaDataTagList.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMetaDataTagListStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "name": attribute(AttributeType.STRING),
        "mechanism": attribute(AttributeType.STRING),
        "serverList": attribute(AttributeType.CUSTOM, custom_class="CList"),
        "userExtensible": attribute(AttributeType.BOOLEAN),
        "customCodeFile": attribute(AttributeType.CUSTOM, custom_class="CDataFile"),
        "queueOptionsFile": attribute(AttributeType.CUSTOM, custom_class="CDataFile"),
        "ccp4Dir": attribute(AttributeType.STRING),
        "tempDir": attribute(AttributeType.STRING),
        "sge_root": attribute(AttributeType.STRING),
        "keyFilename": attribute(AttributeType.STRING),
        "validate": attribute(AttributeType.STRING),
        "timeout": attribute(AttributeType.FLOAT),
        "maxTries": attribute(AttributeType.INT),
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
        'name',
        'mechanism',
        'serverList',
        'userExtensible',
        'ccp4Dir',
        'tempDir',
        'sge_root',
        'keyFilename',
        'validate',
        'customCodeFile',
        'queueOptionsFile'],
    content_qualifiers={
        "mechanism": {'enumerators': ['ssh', 'ssh_shared', 'qsub_local', 'qsub_remote', 'slurm_remote', 'custom'], 'menuText': ['ssh', 'ssh with shared filesystem', 'qsub queue', 'qsub on another machine', 'Slurm on another machine', 'custom'], 'onlyEnumerators': True, 'default': 'ssh'},
        "serverList": {'minLength': 1},
        "userExtensible": {'default': False},
        "customCodeFile": {'allowUndefind': True, 'fileExtensions': ['py']},
        "queueOptionsFile": {'allowUndefind': True},
        "ccp4Dir": {'allowUndefind': True},
        "tempDir": {'allowUndefind': True},
        "sge_root": {'allowUndefind': True},
        "keyFilename": {'allowUndefind': True},
        "validate": {'onlyEnumerators': True, 'default': 'password', 'enumerators': ['password', 'key_filename', 'pass_key_filename']},
        "maxTries": {'default': 2},
    },
)
class CServerGroupStub(CData):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    This is a pure data class stub. Extend it in core/CServerGroup.py
    to add methods and implementation-specific functionality.
    """

    name: Optional[CString] = None
    mechanism: Optional[CString] = None
    serverList: Optional[CList] = None
    userExtensible: Optional[CBoolean] = None
    customCodeFile: Optional[CDataFile] = None
    queueOptionsFile: Optional[CDataFile] = None
    ccp4Dir: Optional[CString] = None
    tempDir: Optional[CString] = None
    sge_root: Optional[CString] = None
    keyFilename: Optional[CString] = None
    validate: Optional[CString] = None
    timeout: Optional[CFloat] = None
    maxTries: Optional[CInt] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CServerGroupStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "platformNode": attribute(AttributeType.STRING),
        "userId": attribute(AttributeType.CUSTOM, custom_class="CUserIdStub"),
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
        "label": 'User id and current machine',
        "toolTip": 'User id as me@myplace.ac.uk and machine name',
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
    contents_order=['platformNode', 'userId'],
)
class CUserAddressStub(CData):
    """
    User id and platform node

    This is a pure data class stub. Extend it in core/CUserAddress.py
    to add methods and implementation-specific functionality.
    """

    platformNode: Optional[CString] = None
    userId: Optional[CUserIdStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUserAddressStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)


@cdata_class(
    attributes={
        "text": attribute(AttributeType.STRING),
        "time": attribute(AttributeType.CUSTOM, custom_class="CTimeStub"),
        "author": attribute(AttributeType.CUSTOM, custom_class="CUserIdStub"),
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
        "label": 'Annotation',
        "toolTip": 'Enter your comments',
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
        "text": {'allowUndefined': True, 'charWidth': -1},
        "time": {'allowUndefined': True, 'default': None},
        "author": {'allowUndefined': True, 'default': None},
    },
)
class CAnnotationStub(CData):
    """
    Annotation text with user id and time

    This is a pure data class stub. Extend it in core/CAnnotation.py
    to add methods and implementation-specific functionality.
    """

    text: Optional[CString] = None
    time: Optional[CTimeStub] = None
    author: Optional[CUserIdStub] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnnotationStub.

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
        "label": 'Machine name',
        "toolTip": 'Hostname as mycomputer.myplace.ac.uk',
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
class CHostnameStub(CHostNameStub):
    """
    Computer name

    This is a pure data class stub. Extend it in core/CHostname.py
    to add methods and implementation-specific functionality.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHostnameStub.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
