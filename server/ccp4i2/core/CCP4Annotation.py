from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType
from ccp4i2.core.base_object.base_classes import CData, CDataFile
from ccp4i2.core.base_object.fundamental_types import CBoolean, CFloat, CInt, CList, CString
"""
Implementation classes for CCP4Annotation.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""




@cdata_class(
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
    contents_order=[
        'pmid',
        'title',
        'authorList',
        'source',
        'url',
        'selected'],
)
class CBibReference(CData):

    pmid: Optional["CInt"] = None
    title: Optional["CString"] = None
    authorList: Optional["CList"] = None
    source: Optional["CString"] = None
    url: Optional["CString"] = None
    selected: Optional["CBoolean"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReference.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Bibliographic reference
    
    Extends CBibReference with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
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
    contents_order=['taskName', 'version', 'title', 'references'],
)
class CBibReferenceGroup(CData):

    taskName: Optional["CString"] = None
    version: Optional["CString"] = None
    title: Optional["CString"] = None
    references: Optional["CList"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReferenceGroup.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Set of bibliographic references for a task
    
    Extends CBibReferenceGroup with implementation-specific methods.
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
        "family": {'default': 'Helvetica'},
        "style": {'onlyEnumerators': True, 'default': 0, 'enumerators': [0, 1, 2], 'menuText': ['normal', 'italic', 'oblique']},
        "pointSize": {'min': 1, 'default': 12},
        "weight": {'min': 0, 'max': 99, 'default': 50, 'allowUndefined': False, 'enumerators': [25, 50, 63, 75, 87], 'menuText': ['light', 'normal', 'demi-bold', 'bold', 'black']},
    },
)
class CFont(CData):

    family: Optional["CString"] = None
    style: Optional["CInt"] = None
    pointSize: Optional["CInt"] = None
    weight: Optional["CInt"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFont.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Simplified Qt font options
    
    Extends CFont with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
class CUserId(CString):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUserId.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A user ID
    
    Extends CUserId with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "min": 0,
        "label": 'Time',
        "toolTip": 'Time and date as hh:mm dd/mm/yyyy',
        "format": '%H:%M %d/%b/%y',
    },
    qualifiers_order=['format'],
)
class CTime(CInt):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTime.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    The time. Uses Python time module
    
    Extends CTime with implementation-specific methods.
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
class CAnnotationList(CList):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnnotationList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A list of annotation
    
    Extends CAnnotationList with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "enumeratorsFunction": None,
        "addEnumeratorFunction": None,
    },
    qualifiers_order=['enumeratorsFunction', 'addEnumeratorFunction'],
    contents_order=['tag'],
)
class CMetaDataTag(CData):

    tag: Optional["CString"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMetaDataTag.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    This class will extend list of enumerators if new value for string is entered
    
    Extends CMetaDataTag with implementation-specific methods.
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
class CDateRange(CData):

    year: Optional["CInt"] = None
    month: Optional["CString"] = None
    day: Optional["CInt"] = None
    yearRange: Optional["CInt"] = None
    monthRange: Optional["CInt"] = None
    dayRange: Optional["CInt"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDateRange.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A date range - may be on a scale of years,months or days
    
    Extends CDateRange with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
)
class CAuthor(CString):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAuthor.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Placeholder for bibliographic author
    
    Extends CAuthor with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
)
class CHostName(CString):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHostName.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Computer name
    
    Extends CHostName with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "listMinLength": 1,
    },
    qualifiers_order=['listMinLength', 'listMaxLength', 'listCompare'],
)
class CMetaDataTagList(CList):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMetaDataTagList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    A list with all items of one CData sub-class
    
    Extends CMetaDataTagList with implementation-specific methods.
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
class CServerGroup(CData):

    name: Optional["CString"] = None
    mechanism: Optional["CString"] = None
    serverList: Optional["CList"] = None
    userExtensible: Optional["CBoolean"] = None
    customCodeFile: Optional["CDataFile"] = None
    queueOptionsFile: Optional["CDataFile"] = None
    ccp4Dir: Optional["CString"] = None
    tempDir: Optional["CString"] = None
    sge_root: Optional["CString"] = None
    keyFilename: Optional["CString"] = None
    validate: Optional["CString"] = None
    timeout: Optional["CFloat"] = None
    maxTries: Optional["CInt"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CServerGroup.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CServerGroup with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
    contents_order=['platformNode', 'userId'],
)
class CUserAddress(CData):

    platformNode: Optional["CString"] = None
    userId: Optional["CUserId"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUserAddress.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    User id and platform node
    
    Extends CUserAddress with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
    content_qualifiers={
        "text": {'allowUndefined': True, 'charWidth': -1},
        "time": {'allowUndefined': True, 'default': None},
        "author": {'allowUndefined': True, 'default': None},
    },
)
class CAnnotation(CData):

    text: Optional["CString"] = None
    time: Optional["CTime"] = None
    author: Optional["CUserId"] = None

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnnotation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Annotation text with user id and time
    
    Extends CAnnotation with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
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
)
class CHostname(CHostName):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CHostname.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    
    Inherits from:
    - CHostname: Metadata and structure
    - CHostName: Shared full-fat methods
    Computer name
    
    Extends CHostname with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

