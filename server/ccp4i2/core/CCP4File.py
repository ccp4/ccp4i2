"""
CCP4File.py --- CData classes.

These were once two classes each: a generated stub carrying the data
model, and an implementation carrying the methods. The generator is no
longer run and the split cost more than it saved --- the two halves
interleaved in the MRO, so an implementation could drop out of its own
subclass's ancestry and `isinstance` would say no to a file that plainly
was one. They are one class now.
"""

import platform
import socket
import time


# Re-export CDataFile for legacy code compatibility
# Many legacy files use "CCP4File.CDataFile" which is actually in base_object
from ccp4i2.core.base_object.cdata_file import CDataFile
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType, content
from ccp4i2.core.base_object.base_classes import CData, CDataFile, CDataFileContent
from ccp4i2.core.base_object.fundamental_types import CInt, CList, CString
from ccp4i2.core.CCP4Annotation import CHostName, CTime, CUserId
from ccp4i2.core.CCP4Data import CUUID


class CProjectId(CUUID):


    """The CCP4i2 database project id - a global unique id"""
    class Meta:
        error_codes = {
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
        }
        qualifiers = {
            "allowUndefined": True,
            "allowUnfound": True,
            "default": None,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProjectId.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CVersion(CString):


    """A (string) version number of the form n.m.i"""
    class Meta:
        error_codes = {
            "101": {
                "description": "Version is not of form n.m or n.m.i"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "default": None,
            "charWidth": 10,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CVersion.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CMmcifData(CDataFileContent):


    """
    Generic mmCIF data.
    This is intended to be a base class for other classes
    specific to coordinates, reflections or geometry data.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMmcifData.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CExePath(CData):

    """A program and where to find it, when it is not on PATH."""

    exeName = content("CString", guiLabel='Program', toolTip='Name of the executable')
    exePath = content(
        "CDataFile",
        mustExist=True,
        allowUndefined=False,
        guiLabel='Path',
        toolTip='Where the executable is, when it is not on PATH')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExePath.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CProjectName(CString):


    """The name of a CCP4i project or directory alias"""
    class Meta:
        error_codes = {
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
        }
        qualifiers = {
            "allowUndefined": True,
            "allowAlias": True,
            "allowUnfound": True,
            "default": None,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CProjectName.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CFileFunction(CString):


    """List of recognised XML file functions"""
    class Meta:
        qualifiers = {
            "enumerators": ['DEF', 'PARAMS', 'LOG', 'PROJECTDIRECTORIES', 'COM', 'REFMAC', 'OUTPUT', 'STATUS', 'PROJECTDATABASE', 'MGSCENE', 'JOBSERVERSTATUS', 'WORKFLOW', 'COMFILEPATCH', 'CUSTOMTASK', 'IMPORTEDJOB', 'I1SUPPLEMENT', 'ASUCONTENT', 'UNKNOWN'],
            "onlyEnumerators": True,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFileFunction.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CExportedFile(CData):


    """A file written out of a project by an export."""

    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    exportId = content(
        "CUUID",
        guiLabel='Export id',
        toolTip='Identifier of the export this file was written for')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExportedFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CExportedFileList(CList):


    """A list with all items of one CData sub-class"""
    class Meta:
        qualifiers = {
            "listMinLength": 0,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExportedFileList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CExePathList(CList):


    """A list with all items of one CData sub-class"""
    class Meta:
        qualifiers = {
            "listMinLength": 1,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExePathList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CSearchPath(CData):


    """A directory searched for files."""

    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    name = content("CString", guiLabel='Name', toolTip='Name of this item')
    path = content("CDataFile", guiLabel='Path', toolTip='Directory searched for files')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSearchPath.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CSearchPathList(CList):


    """A list with all items of one CData sub-class"""
    class Meta:
        qualifiers = {
            "listMinLength": 0,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSearchPathList.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CFilePath(CString):


    """A file path"""
    class Meta:
        error_codes = {
            "101": {
                "description": "Invalid characters in file name"
            },
            "102": {
                "severity": 2,
                "description": "Invalid characters in file name"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "allowedCharacters": '',
            "allowedCharactersMode": 1,
            "default": None,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CFilePath.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CI2XmlHeader(CData):


    """Container for header info from XML file"""
    class Meta:
        error_codes = {
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
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    function = content(
        "CFileFunction",
        guiLabel='Function',
        toolTip='What the file is for, e.g. the parameters of a job or its report')
    userId = content(
        "CUserId",
        guiLabel='User',
        toolTip='Login name of the user who created the file')
    hostName = content(
        "CHostName",
        guiLabel='Host',
        toolTip='Machine on which the file was created')
    creationTime = content(
        "CTime",
        guiLabel='Created',
        toolTip='When the file was written')
    pluginName = content(
        "CString",
        guiLabel='Task',
        toolTip='Name of the task that wrote the file')
    pluginVersion = content(
        "CVersion",
        guiLabel='Task version',
        toolTip='Version of that task, so a file can be read by the code that wrote it')
    pluginTitle = content(
        "CString",
        guiLabel='Task title',
        toolTip='Human-readable title of the task')
    projectName = content(
        "CProjectName",
        guiLabel='Project',
        toolTip='Name of the project the file belongs to')
    projectId = content(
        "CProjectId",
        allowUnfound=True,
        guiLabel='Project id',
        toolTip='Identifier of that project')
    jobId = content(
        "CUUID",
        guiLabel='Job id',
        toolTip='Identifier of the job that produced the file')
    jobNumber = content(
        "CString",
        guiLabel='Job number',
        toolTip="The job's number within its project, as shown in the job list")
    comment = content(
        "CString",
        guiLabel='Comment',
        toolTip='Free text recorded with the file')
    OS = content("CString", guiLabel='OS', toolTip='Platform the file was created on')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CI2XmlHeader.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    def setCurrent(self):
        """
        Set header fields to current values: time, hostname, OS, CCP4 version.

        This populates standard header metadata that should be set when creating
        a new XML file.
        """

        # Set creation time to now (Unix timestamp as integer)
        if hasattr(self, 'creationTime'):
            self.creationTime.set(int(time.time()))

        # Set hostname
        if hasattr(self, 'hostName'):
            self.hostName.set(socket.gethostname())

        # Set OS
        if hasattr(self, 'OS'):
            self.OS.set(f"{platform.system()} {platform.release()}")


class CXmgrDataFile(CDataFile):



    """An xmgr format file. This is the input format for xmgrace, as output by scala or aimless"""
    class Meta:
        qualifiers = {
            "mimeTypeName": 'application/grace',
            "fileExtensions": ['xmgr'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXmgrDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CTextDataFile(CDataFile):



    """A text data file"""
    class Meta:
        qualifiers = {
            "mimeTypeName": '"text/plain"',
            "mimeTypeDescription": 'Standard plain text',
            "fileLabel": None,
            "fileExtensions": ['txt', 'log'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTextDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CDataReflFile(CDataFile):



    """Reflection file from DIALS"""
    class Meta:
        qualifiers = {
            "guiLabel": 'Reflections from DIALS',
            "fileExtensions": ['refl'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReflFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CPostscriptDataFile(CDataFile):



    """A postscript format file"""
    class Meta:
        qualifiers = {
            "mimeTypeName": 'application/postscript',
            "fileExtensions": ['ps'],
            "guiLabel": 'Postscript file',
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPostscriptDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CYmlFile(CDataFile):



    """A yml data file"""
    class Meta:
        qualifiers = {
            "mimeTypeName": '"text/plain"',
            "mimeTypeDescription": 'Standard plain text',
            "guiLabel": 'yml file',
            "fileExtensions": ['yml'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CYmlFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CXmlDataFile(CDataFile):



    """A reference to an XML file"""
    class Meta:
        error_codes = {
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
        }
        qualifiers = {
            "fileExtensions": ['xml'],
            "saveToDb": False,
            "mimeTypeName": 'application/xml',
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXmlDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    def saveFile(self, bodyEtree=None):
        """
        Save XML content to the file path specified by this object.

        Args:
            bodyEtree: Optional ElementTree.Element to use as the body content.
                      If None, creates an empty file.

        Returns:
            bool: True if save was successful
        """
        import xml.etree.ElementTree as ET
        from pathlib import Path

        # If bodyEtree is provided, write it directly
        if bodyEtree is not None:
            tree = ET.ElementTree(bodyEtree)
            file_path = Path(self.getFullPath())

            # Ensure directory exists
            file_path.parent.mkdir(parents=True, exist_ok=True)

            # Write with pretty formatting
            ET.indent(tree, space='  ')
            tree.write(file_path, encoding='utf-8', xml_declaration=True)

            return True
        else:
            # Create empty root if no body provided
            root = ET.Element('data')
            tree = ET.ElementTree(root)
            file_path = Path(self.getFullPath())

            file_path.parent.mkdir(parents=True, exist_ok=True)
            tree.write(file_path, encoding='utf-8', xml_declaration=True)

            return True


class CMmcifDataFile(CDataFile):



    """
    A generic mmCIF format file.
    This is intended to be a base class for other classes
    specific to coordinates, reflections or geometry data.
    """
    class Meta:
        qualifiers = {
            "fileExtensions": ['cif', 'ent'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMmcifDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CPDFDataFile(CDataFile):



    """An PDF format file"""
    class Meta:
        qualifiers = {
            "mimeTypeName": 'application/x-pdf',
            "fileExtensions": ['pdf'],
            "guiLabel": 'PDF file',
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPDFDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CSceneDataFile(CDataFile):



    """An xml format file for defining scene in CCP4mg."""
    class Meta:
        qualifiers = {
            "fileLabel": 'scene',
            "mimeTypeName": 'application/CCP4-scene',
            "mimeTypeDescription": 'CCP4mg scene file',
            "guiLabel": 'CCP4mg scene',
            "fileExtensions": ['scene.xml'],
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSceneDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CI2XmlDataFile(CXmlDataFile):


    """A reference to an XML file with CCP4i2 Header"""
    class Meta:
        error_codes = {
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
        }
        qualifiers = {
            "fileExtensions": ['xml'],
            "autoLoadHeader": True,
        }
    header = content(
        "CI2XmlHeader",
        guiLabel='Header',
        toolTip='Provenance recorded at the top of a CCP4i2 XML file')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CI2XmlDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    def saveFile(self, bodyEtree=None, useLXML=None):
        """
        Save the XML file with header and body structure.

        CCP4i2 XML files have a standard structure:
        <ccp4i2>
          <header>...</header>
          <body>...</body>
        </ccp4i2>

        Args:
            bodyEtree: Optional ElementTree element for the body content.
                      If not provided, an empty body will be created.
            useLXML: Ignored - kept for backward compatibility with legacy code.
        """
        import xml.etree.ElementTree as ET
        import traceback
        from pathlib import Path
        import sys
        from datetime import datetime
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Called at {traceback.format_stack()}")
        # Create root element
        root = ET.Element('ccp4i2')
        #traceback.print_stack(file=sys.stdout)
        # Add header
        if hasattr(self, 'header') and self.header is not None:
            header_elem = self.header.getEtree()
            if header_elem is not None:
                root.append(header_elem)
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Header etree: {ET.tostring(header_elem, encoding='unicode')}")
        # Add body
        if bodyEtree is not None:
            # If bodyEtree is provided, use it as the body
            if bodyEtree.tag in ('body', 'ccp4i2_body'):
                root.append(bodyEtree)
            else:
                # Create ccp4i2_body and unwrap the container (CCP4i2 format has no <container> wrapper)
                # The container's children (inputData, outputData, etc.) go directly into ccp4i2_body
                body = ET.Element('ccp4i2_body')
                # If bodyEtree is a container, unwrap it and append its children
                if bodyEtree.tag.lower() == 'container':
                    for child in bodyEtree:
                        body.append(child)
                else:
                    # For non-container elements, append as-is
                    body.append(bodyEtree)
                root.append(body)
        else:
            # Create empty body
            body = ET.Element('ccp4i2_body')
            root.append(body)
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Root etree: {ET.tostring(root, encoding='unicode') if root is not None else 'None'}")

        # Create tree and write to file
        tree = ET.ElementTree(root)
        full_path_str = self.getFullPath()
        #print(f"[DEBUG CI2XmlDataFile.saveFile] getFullPath() returned: '{full_path_str}'")

        if not full_path_str or full_path_str.strip() == '':
            #pass  # DEBUG: print(f"[DEBUG CI2XmlDataFile.saveFile] ERROR: getFullPath() returned empty string!")
            #pass  # DEBUG: print(f"[DEBUG CI2XmlDataFile.saveFile] baseName.value: {self.baseName.value if hasattr(self, 'baseName') else 'N/A'}")
            return False

        file_path = Path(full_path_str)
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Writing XML file to: {file_path}")
        
        # Ensure directory exists
        file_path.parent.mkdir(parents=True, exist_ok=True)

        # Write with pretty formatting
        ET.indent(root, space='  ')
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Root etree: {ET.tostring(root, encoding='unicode') if root is not None else 'None'}")
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Etree that will be written: {ET.tostring(root, encoding='unicode')}")
        # Instead of tree.write(file_path, encoding='utf-8', xml_declaration=True)
        xml_string = ET.tostring(root, encoding='unicode')
        #print(f"[DEBUG] Type of xml_string: {type(xml_string)}")
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Final XML string to write:\n{xml_string}")
        with open(file_path, 'w', encoding='utf-8') as f:  # Note: 'wb' for bytes
            f.write(xml_string)
            f.flush()
        #with open(file_path, 'r', encoding='utf-8') as f:
        #    content = f.read()
        #    #print(f"[DEBUG CI2XmlDataFile.saveFile] Written file content:\n{content}")
        #tree.write(file_path, encoding='utf-8', xml_declaration=True)
        #print(f"[DEBUG CI2XmlDataFile.saveFile] Successfully wrote file to: {file_path}")

        return True


class CEBIValidationXMLDataFile(CXmlDataFile):



    """An XLM file returned from the EBI validation server"""
    class Meta:
        qualifiers = {
            "mimeTypeName": 'application/EBI-validation-xml',
            "fileExtensions": ['xml'],
            "guiLabel": 'EBI Validation XML',
        }
        content_qualifiers = {
            "subType": {'default': None},
            "contentFlag": {'min': 0, 'default': None},
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEBIValidationXMLDataFile.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass

