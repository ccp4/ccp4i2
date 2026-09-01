from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType, content
from ccp4i2.core.base_object.base_classes import CData, CDataFile
from ccp4i2.core.base_object.fundamental_types import CBoolean, CFloat, CInt, CList, CString
"""
Implementation classes for CCP4Annotation.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""




class CBibReference(CData):


    class Meta:
        error_codes = {
            "101": {
                "description": "Failed to load Medline data"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    pmid = content("CInt")
    title = content("CString")
    authorList = content("CList")
    source = content("CString")
    url = content("CString")
    selected = content("CBoolean")

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


class CBibReferenceGroup(CData):


    class Meta:
        error_codes = {
            "100": {
                "description": "Failed attempting to load MedLine file - file not found"
            },
            "101": {
                "description": "Failed attempting to find references file"
            },
            "102": {
                "description": "Error copying file"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    taskName = content("CString")
    version = content("CString")
    title = content("CString")
    references = content("CList")

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


class CFont(CData):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    family = content("CString", default='Helvetica')
    style = content("CInt",
                    onlyEnumerators=True,
                    default=0,
                    enumerators=[0, 1, 2],
                    menuText=['normal', 'italic', 'oblique'])
    pointSize = content("CInt", min=1, default=12)
    weight = content("CInt",
                     min=0,
                     max=99,
                     default=50,
                     allowUndefined=False,
                     enumerators=[25, 50, 63, 75, 87],
                     menuText=['light', 'normal', 'demi-bold', 'bold', 'black'])

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


class CUserId(CString):


    class Meta:
        qualifiers = {
            "label": 'User id',
            "toolTip": 'User id as me@myplace.ac.uk',
        }
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


class CTime(CInt):


    class Meta:
        qualifiers = {
            "min": 0,
            "label": 'Time',
            "toolTip": 'Time and date as hh:mm dd/mm/yyyy',
            "format": '%H:%M %d/%b/%y',
        }
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


class CAnnotationList(CList):


    class Meta:
        qualifiers = {
            "listMinLength": 0,
        }
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


class CMetaDataTag(CData):


    class Meta:
        qualifiers = {
            "enumeratorsFunction": None,
            "addEnumeratorFunction": None,
        }
    tag = content("CString")

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


class CDateRange(CData):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    year = content("CInt", enumerators=[])
    month = content("CString",
                    onlyEnumerators=True,
                    enumerators=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'],
                    default='January')
    day = content("CInt", default=1, min=1, max=31)
    yearRange = content("CInt", default=0, min=0, max=100)
    monthRange = content("CInt", default=0, min=0, max=12)
    dayRange = content("CInt", default=0, min=0, max=30)

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


class CAuthor(CString):


    class Meta:
        qualifiers = {
            "minLength": None,
            "maxLength": None,
            "enumerators": [],
            "menuText": [],
            "onlyEnumerators": False,
            "charWidth": -1,
            "allowedCharsCode": 0,
        }
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


class CHostName(CString):


    class Meta:
        qualifiers = {
            "label": 'Machine name',
            "toolTip": 'Hostname as mycomputer.myplace.ac.uk',
        }
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


class CMetaDataTagList(CList):


    class Meta:
        qualifiers = {
            "listMinLength": 1,
        }
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


class CServerGroup(CData):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    name = content("CString")
    mechanism = content("CString",
                        enumerators=['ssh', 'ssh_shared', 'qsub_local', 'qsub_remote', 'slurm_remote', 'custom'],
                        menuText=['ssh', 'ssh with shared filesystem', 'qsub queue', 'qsub on another machine', 'Slurm on another machine', 'custom'],
                        onlyEnumerators=True,
                        default='ssh')
    serverList = content("CList", minLength=1)
    userExtensible = content("CBoolean", default=False)
    ccp4Dir = content("CString", allowUndefind=True)
    tempDir = content("CString", allowUndefind=True)
    sge_root = content("CString", allowUndefind=True)
    keyFilename = content("CString", allowUndefind=True)
    validate = content("CString",
                       onlyEnumerators=True,
                       default='password',
                       enumerators=['password', 'key_filename', 'pass_key_filename'])
    customCodeFile = content("CDataFile", allowUndefind=True, fileExtensions=['py'])
    queueOptionsFile = content("CDataFile", allowUndefind=True)
    timeout = content("CFloat")
    maxTries = content("CInt", default=2)

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


class CUserAddress(CData):


    class Meta:
        qualifiers = {
            "label": 'User id and current machine',
            "toolTip": 'User id as me@myplace.ac.uk and machine name',
        }
    platformNode = content("CString")
    userId = content("CUserId")

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


class CAnnotation(CData):


    class Meta:
        qualifiers = {
            "label": 'Annotation',
            "toolTip": 'Enter your comments',
        }
    text = content("CString", allowUndefined=True, charWidth=-1)
    time = content("CTime", allowUndefined=True, default=None)
    author = content("CUserId", allowUndefined=True, default=None)

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


class CHostname(CHostName):


    class Meta:
        qualifiers = {
            "label": 'Machine name',
            "toolTip": 'Hostname as mycomputer.myplace.ac.uk',
        }
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

