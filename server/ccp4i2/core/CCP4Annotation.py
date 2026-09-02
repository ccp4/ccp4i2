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


    """
    Bibliographic reference
    
    Extends CBibReference with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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
    pmid = content("CInt", guiLabel='PMID', toolTip='Pubmed ID')
    title = content("CString", guiLabel='Title', toolTip='Title of the referenced work')
    authorList = content("CList", guiLabel='Author list', toolTip='List of authors')
    source = content(
        "CString",
        guiLabel='Source',
        toolTip='Journal or other publication the work appeared in')
    url = content(
        "CString",
        guiLabel='URL',
        toolTip='Web address of the referenced work')
    selected = content(
        "CBoolean",
        guiLabel='Selected',
        toolTip='Whether this reference is cited in the job report')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReference.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CBibReferenceGroup(CData):


    """
    Set of bibliographic references for a task
    
    Extends CBibReferenceGroup with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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
    taskName = content(
        "CString",
        guiLabel='Task',
        toolTip='Task these references should be cited for')
    version = content(
        "CString",
        guiLabel='Version',
        toolTip='Version of the task the references apply to')
    title = content(
        "CString",
        guiLabel='Title',
        toolTip='Heading for this group of references')
    references = content(
        "CList",
        guiLabel='References',
        toolTip='The bibliographic references in this group')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CBibReferenceGroup.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CFont(CData):


    """
    Simplified Qt font options
    
    Extends CFont with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    family = content(
        "CString",
        default='Helvetica',
        guiLabel='Font',
        toolTip='Name of the font family, e.g. Helvetica')
    style = content(
        "CInt",
        onlyEnumerators=True,
        default=0,
        enumerators=[0, 1, 2],
        menuText=['normal', 'italic', 'oblique'],
        guiLabel='Style',
        toolTip='Normal, italic or oblique')
    pointSize = content(
        "CInt",
        min=1,
        default=12,
        guiLabel='Size',
        toolTip='Font size in points')
    weight = content(
        "CInt",
        min=0,
        max=99,
        default=50,
        allowUndefined=False,
        enumerators=[25, 50, 63, 75, 87],
        menuText=['light', 'normal', 'demi-bold', 'bold', 'black'],
        guiLabel='Weight',
        toolTip='Stroke weight, from light through normal to bold')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFont.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CUserId(CString):


    """
    A user ID
    
    Extends CUserId with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CTime(CInt):


    """
    The time. Uses Python time module
    
    Extends CTime with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CAnnotationList(CList):


    """
    A list of annotation
    
    Extends CAnnotationList with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CMetaDataTag(CData):


    """
    This class will extend list of enumerators if new value for string is entered
    
    Extends CMetaDataTag with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "enumeratorsFunction": None,
            "addEnumeratorFunction": None,
        }
    tag = content(
        "CString",
        guiLabel='Tag',
        toolTip='Metadata keyword recorded against the item')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMetaDataTag.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CDateRange(CData):


    """
    A date range - may be on a scale of years,months or days
    
    Extends CDateRange with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    year = content("CInt", enumerators=[], guiLabel='Year', toolTip='The year')
    month = content(
        "CString",
        onlyEnumerators=True,
        enumerators=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'],
        default='January',
        guiLabel='Month',
        toolTip='Month of the year')
    day = content(
        "CInt",
        default=1,
        min=1,
        max=31,
        guiLabel='Day',
        toolTip='Day of the month')
    yearRange = content(
        "CInt",
        default=0,
        min=0,
        max=100,
        guiLabel='Years',
        toolTip='Number of years the range spans')
    monthRange = content(
        "CInt",
        default=0,
        min=0,
        max=12,
        guiLabel='Months',
        toolTip='Number of months the range spans')
    dayRange = content(
        "CInt",
        default=0,
        min=0,
        max=30,
        guiLabel='Days',
        toolTip='Number of days the range spans')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDateRange.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CAuthor(CString):


    """
    Placeholder for bibliographic author
    
    Extends CAuthor with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CHostName(CString):


    """
    Computer name
    
    Extends CHostName with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CMetaDataTagList(CList):


    """
    A list with all items of one CData sub-class
    
    Extends CMetaDataTagList with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass


class CServerGroup(CData):


    """
    Extends CServerGroup with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    name = content(
        "CString",
        guiLabel='Name',
        toolTip='Name of this server configuration')
    mechanism = content(
        "CString",
        enumerators=['ssh', 'ssh_shared', 'qsub_local', 'qsub_remote', 'slurm_remote', 'custom'],
        menuText=['ssh', 'ssh with shared filesystem', 'qsub queue', 'qsub on another machine', 'Slurm on another machine', 'custom'],
        onlyEnumerators=True,
        default='ssh',
        guiLabel='Mechanism',
        toolTip='How jobs are submitted, e.g. by ssh or to a queue')
    serverList = content(
        "CList",
        minLength=1,
        guiLabel='Servers',
        toolTip='Machines jobs may be sent to')
    userExtensible = content(
        "CBoolean",
        default=False,
        guiLabel='User extensible',
        toolTip='Whether a user may add servers to this group')
    ccp4Dir = content(
        "CString",
        allowUndefind=True,
        guiLabel='CCP4 directory',
        toolTip='Location of the CCP4 installation on the remote machine')
    tempDir = content(
        "CString",
        allowUndefind=True,
        guiLabel='Temporary directory',
        toolTip="Directory for a job's scratch files on the remote machine")
    sge_root = content(
        "CString",
        allowUndefind=True,
        guiLabel='SGE root',
        toolTip='Root of the Grid Engine installation, for queue submission')
    keyFilename = content(
        "CString",
        allowUndefind=True,
        guiLabel='Key file',
        toolTip='Private key used to authenticate to the remote machine')
    validate = content(
        "CString",
        onlyEnumerators=True,
        default='password',
        enumerators=['password', 'key_filename', 'pass_key_filename'],
        guiLabel='Validate',
        toolTip='Command used to check the remote machine is usable')
    customCodeFile = content(
        "CDataFile",
        allowUndefind=True,
        fileExtensions=['py'],
        guiLabel='Custom code',
        toolTip='Script supplying site-specific submission logic')
    queueOptionsFile = content(
        "CDataFile",
        allowUndefind=True,
        guiLabel='Queue options',
        toolTip='File of options passed to the queueing system')
    timeout = content(
        "CFloat",
        guiLabel='Timeout',
        toolTip='How long to wait for the remote machine before giving up, in seconds')
    maxTries = content(
        "CInt",
        default=2,
        guiLabel='Max. tries',
        toolTip='How many times to retry a failed submission')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CServerGroup.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CUserAddress(CData):


    """
    User id and platform node
    
    Extends CUserAddress with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "label": 'User id and current machine',
            "toolTip": 'User id as me@myplace.ac.uk and machine name',
        }
    platformNode = content(
        "CString",
        guiLabel='Machine',
        toolTip='Network name of the machine')
    userId = content("CUserId", guiLabel='User', toolTip='Login name of the user')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CUserAddress.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CAnnotation(CData):


    """
    Annotation text with user id and time
    
    Extends CAnnotation with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "label": 'Annotation',
            "toolTip": 'Enter your comments',
        }
    text = content(
        "CString",
        allowUndefined=True,
        charWidth=-1,
        guiLabel='Text',
        toolTip='Text of annotation')
    time = content(
        "CTime",
        allowUndefined=True,
        default=None,
        guiLabel='Timestamp',
        toolTip='Time of annotation')
    author = content(
        "CUserId",
        allowUndefined=True,
        default=None,
        guiLabel='Author',
        toolTip='Annotation author')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAnnotation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CHostname(CHostName):


    """
    
    Inherits from:
    - CHostname: Metadata and structure
    - CHostName: Shared full-fat methods
    Computer name
    
    Extends CHostname with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
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

    # Add your methods here
    pass

