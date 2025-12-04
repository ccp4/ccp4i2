"""
Implementation classes for CCP4File.py

Extends stub classes from core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from __future__ import annotations
from typing import Optional, Any

from core.cdata_stubs.CCP4File import CDataReflFileStub, CEBIValidationXMLDataFileStub, CExePathStub, CExePathListStub, CExportedFileStub, CExportedFileListStub, CFileFunctionStub, CFilePathStub, CI2XmlDataFileStub, CI2XmlHeaderStub, CMmcifDataStub, CMmcifDataFileStub, CPDFDataFileStub, CPostscriptDataFileStub, CProjectIdStub, CProjectNameStub, CSceneDataFileStub, CSearchPathStub, CSearchPathListStub, CTextDataFileStub, CVersionStub, CXmgrDataFileStub, CXmlDataFileStub, CYmlFileStub

# Re-export CDataFile for legacy code compatibility
# Many legacy files use "CCP4File.CDataFile" which is actually in base_object
from core.base_object.cdata_file import CDataFile


class CDataReflFile(CDataReflFileStub):
    """
    Reflection file from DIALS
    
    Extends CDataReflFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CEBIValidationXMLDataFile(CEBIValidationXMLDataFileStub):
    """
    An XLM file returned from the EBI validation server 
    
    Extends CEBIValidationXMLDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExePath(CExePathStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CExePathStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExePathList(CExePathListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CExePathListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExportedFile(CExportedFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CExportedFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExportedFileList(CExportedFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CExportedFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFileFunction(CFileFunctionStub):
    """
    List of recognised XML file functions
    
    Extends CFileFunctionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFilePath(CFilePathStub):
    """
    A file path
    
    Extends CFilePathStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CI2XmlDataFile(CI2XmlDataFileStub):
    """
    A reference to an XML file with CCP4i2 Header

    Extends CI2XmlDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

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


class CI2XmlHeader(CI2XmlHeaderStub):
    """
    Container for header info from XML file

    Extends CI2XmlHeaderStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def setCurrent(self):
        """
        Set header fields to current values: time, hostname, OS, CCP4 version.

        This populates standard header metadata that should be set when creating
        a new XML file.
        """
        import time
        import socket
        import platform
        import sys

        # Set creation time to now (Unix timestamp as integer)
        if hasattr(self, 'creationTime'):
            self.creationTime.set(int(time.time()))

        # Set hostname
        if hasattr(self, 'hostName'):
            self.hostName.set(socket.gethostname())

        # Set OS
        if hasattr(self, 'OS'):
            self.OS.set(f"{platform.system()} {platform.release()}")

        # Set CCP4i version (Python version as proxy for now)
        if hasattr(self, 'ccp4iVersion'):
            python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            self.ccp4iVersion.set(python_version)


class CMmcifData(CMmcifDataStub):
    """
    Generic mmCIF data.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.
    
    Extends CMmcifDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMmcifDataFile(CMmcifDataFileStub):
    """
    A generic mmCIF format file.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.
    
    Extends CMmcifDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPDFDataFile(CPDFDataFileStub):
    """
    An PDF format file
    
    Extends CPDFDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPostscriptDataFile(CPostscriptDataFileStub):
    """
    A postscript format file
    
    Extends CPostscriptDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CProjectId(CProjectIdStub):
    """
    The CCP4i2 database project id - a global unique id
    
    Extends CProjectIdStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CProjectName(CProjectNameStub):
    """
    The name of a CCP4i project or directory alias
    
    Extends CProjectNameStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSceneDataFile(CSceneDataFileStub):
    """
    An xml format file for defining scene in CCP4mg.
    
    Extends CSceneDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSearchPath(CSearchPathStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CSearchPathStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSearchPathList(CSearchPathListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CSearchPathListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CTextDataFile(CTextDataFileStub):
    """
    A text data file
    
    Extends CTextDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CVersion(CVersionStub):
    """
    A (string) version number of the form n.m.i
    
    Extends CVersionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CXmgrDataFile(CXmgrDataFileStub):
    """
    An xmgr format file. This is the input format for xmgrace, as output by scala or aimless
    
    Extends CXmgrDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CXmlDataFile(CXmlDataFileStub):
    """
    A reference to an XML file

    Extends CXmlDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

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


class CYmlFile(CYmlFileStub):
    """
    A yml data file
    
    Extends CYmlFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

