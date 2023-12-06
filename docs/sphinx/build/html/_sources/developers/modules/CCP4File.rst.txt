CDataFile: A data file - expected to have associated class for file contents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  |                  |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CDataFileContent: Base class for classes holding file contents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__

.. _CCP4File-CCP4File-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-1:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CExePath
~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CExePath <./CCP4File.html#CExePath>`__

.. _CCP4File-CCP4File-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

======= =========================================
exePath `CDataFile <./CCP4File.html#CDataFile>`__
exeName `CString <./CCP4Data.html#CString>`__
======= =========================================

.. _CCP4File-CCP4File-qualifiers-for-class-2:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CExePathList
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4File.CExePathList <./CCP4File.html#CExePathList>`__

.. _CCP4File-CCP4File-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CExePath <./CCP4File.html#CExePath>`__

.. _CCP4File-CCP4File-qualifiers-for-class-3:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| listMinLength  | *int*  | **1**            | Inclusive        |
|                |        |                  | minimum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listMaxLength  | *int*  | **               | Inclusive        |
|                |        | NotImplemented** | maximum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listCompare    | *int*  | **               | If has value     |
|                |        | NotImplemented** | 1/-1 consecutive |
|                |        |                  | items in list    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than preceeding  |
|                |        |                  | item. The list   |
|                |        |                  | item class must  |
|                |        |                  | have a           |
|                |        |                  | \__cmp__()       |
|                |        |                  | method.          |
+----------------+--------+------------------+------------------+

CExportedFile
~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CExportedFile <./CCP4File.html#CExportedFile>`__

.. _CCP4File-CCP4File-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

======== =================================
exportId `CUUID <./CCP4Data.html#CUUID>`__
======== =================================

.. _CCP4File-CCP4File-qualifiers-for-class-4:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CExportedFileList
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4File.CExportedFileList <./CCP4File.html#CExportedFileList>`__

.. _CCP4File-CCP4File-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CExportedFile <./CCP4File.html#CExportedFile>`__

.. _CCP4File-CCP4File-qualifiers-for-class-5:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| listMinLength  | *int*  | **0**            | Inclusive        |
|                |        |                  | minimum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listMaxLength  | *int*  | **               | Inclusive        |
|                |        | NotImplemented** | maximum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listCompare    | *int*  | **               | If has value     |
|                |        | NotImplemented** | 1/-1 consecutive |
|                |        |                  | items in list    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than preceeding  |
|                |        |                  | item. The list   |
|                |        |                  | item class must  |
|                |        |                  | have a           |
|                |        |                  | \__cmp__()       |
|                |        |                  | method.          |
+----------------+--------+------------------+------------------+

CFileFunction: List of recognised XML file functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4File.CFileFunction <./CCP4File.html#CFileFunction>`__

.. _CCP4File-CCP4File-contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-6:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| charWidth        | *int*  | **-1**           |                  |
+------------------+--------+------------------+------------------+
| minLength        | *int*  | **None**         | Minimum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| maxLength        | *int*  | **None**         | Maximum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| onlyEnumerators  | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **['DEF',        |                  |
|                  |        | 'PARAMS', 'LOG', |                  |
|                  |        | 'PROJ            |                  |
|                  |        | ECTDIRECTORIES', |                  |
|                  |        | 'COM', 'REFMAC', |                  |
|                  |        | 'OUTPUT',        |                  |
|                  |        | 'STATUS',        |                  |
|                  |        | 'P               |                  |
|                  |        | ROJECTDATABASE', |                  |
|                  |        | 'MGSCENE',       |                  |
|                  |        | 'J               |                  |
|                  |        | OBSERVERSTATUS', |                  |
|                  |        | 'WORKFLOW',      |                  |
|                  |        | 'COMFILEPATCH',  |                  |
|                  |        | 'CUSTOMTASK',    |                  |
|                  |        | 'IMPORTEDJOB',   |                  |
|                  |        | 'I1SUPPLEMENT',  |                  |
|                  |        | 'UNKNOWN']**     |                  |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | equivalent to    |
|                  |        |                  | the enumerators  |
|                  |        |                  | that will appear |
|                  |        |                  | in the GUI       |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+

CFilePath: A file path
~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4File.CFilePath <./CCP4File.html#CFilePath>`__

.. _CCP4File-CCP4File-contents-of-class-7:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-7:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if allow    |
|                  |        |                  | undefined value  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| default          | *str*  | **None**         | Default file     |
|                  |        |                  | path             |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| charWidth        | *int*  | **-1**           |                  |
+------------------+--------+------------------+------------------+
| minLength        | *int*  | **None**         | Minimum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| maxLength        | *int*  | **None**         | Maximum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| onlyEnumerators  | *bool* | **False**        | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **[]**           | A list of        |
|                  |        |                  | allowed or       |
|                  |        |                  | recommended      |
|                  |        |                  | values for       |
|                  |        |                  | string           |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | equivalent to    |
|                  |        |                  | the enumerators  |
|                  |        |                  | that will appear |
|                  |        |                  | in the GUI       |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+
| a                | *str*  |                  | Set of           |
| llowedCharacters |        |                  | characters       |
|                  |        |                  | allowed in file  |
|                  |        |                  | name             |
+------------------+--------+------------------+------------------+
| allow            | *int*  | **1**            | Handling of      |
| edCharactersMode |        |                  | violation of     |
|                  |        |                  | allowed          |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+

CI2XmlDataFile: A reference to an XML file with CCP4i2 Header
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CXmlDataFile <./CCP4File.html#CXmlDataFile>`__ ->
`CCP4File.CI2XmlDataFile <./CCP4File.html#CI2XmlDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-8:

Contents of class:
^^^^^^^^^^^^^^^^^^

======== ===============================================
project  `CProjectId <./CCP4File.html#CProjectId>`__
header   `CI2XmlHeader <./CCP4File.html#CI2XmlHeader>`__
baseName `CFilePath <./CCP4File.html#CFilePath>`__
relPath  `CFilePath <./CCP4File.html#CFilePath>`__
======== ===============================================

.. _CCP4File-CCP4File-qualifiers-for-class-8:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['xml']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  |                  |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+
| autoLoadHeader   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+

CI2XmlHeader: Container for header info from XML file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CI2XmlHeader <./CCP4File.html#CI2XmlHeader>`__

.. _CCP4File-CCP4File-contents-of-class-9:

Contents of class:
^^^^^^^^^^^^^^^^^^

============= =================================================
function      `CFileFunction <./CCP4File.html#CFileFunction>`__
comment       `CString <./CCP4Data.html#CString>`__
projectName   `CProjectName <./CCP4File.html#CProjectName>`__
hostName      `CHostName <./CCP4Annotation.html#CHostName>`__
userId        `CUserId <./CCP4Annotation.html#CUserId>`__
ccp4iVersion  `CVersion <./CCP4File.html#CVersion>`__
jobId         `CUUID <./CCP4Data.html#CUUID>`__
pluginVersion `CVersion <./CCP4File.html#CVersion>`__
jobNumber     `CString <./CCP4Data.html#CString>`__
projectId     `CProjectId <./CCP4File.html#CProjectId>`__
creationTime  `CTime <./CCP4Annotation.html#CTime>`__
pluginName    `CString <./CCP4Data.html#CString>`__
pluginTitle   `CString <./CCP4Data.html#CString>`__
OS            `CString <./CCP4Data.html#CString>`__
============= =================================================

.. _CCP4File-CCP4File-qualifiers-for-class-9:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CMmcifData: Generic mmCIF data. This is intended to be a base class for other classes specific to coordinates, reflections or geometry data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4File.CMmcifData <./CCP4File.html#CMmcifData>`__

.. _CCP4File-CCP4File-contents-of-class-10:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-10:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CMmcifDataFile: A generic mmCIF format file. This is intended to be a base class for other classes specific to coordinates, reflections or geometry data.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CMmcifDataFile <./CCP4File.html#CMmcifDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-11:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-11:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['cif',        | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  |                  |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CPostscriptDataFile: An postscript format file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CPostscriptDataFile <./CCP4File.html#CPostscriptDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-12:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-12:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['ps']**       | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/postscript** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CProjectId: The CCP4i2 database project id - a global unique id
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CUUID <./CCP4Data.html#CUUID>`__ ->
`CCP4File.CProjectId <./CCP4File.html#CProjectId>`__

.. _CCP4File-CCP4File-contents-of-class-13:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-13:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if allow    |
|                  |        |                  | undefined value  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| default          | *str*  | **None**         |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| charWidth        | *int*  | **-1**           |                  |
+------------------+--------+------------------+------------------+
| minLength        | *int*  | **None**         | Minimum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| maxLength        | *int*  | **None**         | Maximum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| onlyEnumerators  | *bool* | **False**        | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **[]**           | A list of        |
|                  |        |                  | allowed or       |
|                  |        |                  | recommended      |
|                  |        |                  | values for       |
|                  |        |                  | string           |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | equivalent to    |
|                  |        |                  | the enumerators  |
|                  |        |                  | that will appear |
|                  |        |                  | in the GUI       |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+
| allowUnfound     | *bool* | **True**         | Flag if allow    |
|                  |        |                  | unfound project  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+

CProjectName: The name of a CCP4i project or directory alias
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4File.CProjectName <./CCP4File.html#CProjectName>`__

.. _CCP4File-CCP4File-contents-of-class-14:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-14:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if allow    |
|                  |        |                  | undefined value  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| default          | *str*  | **None**         |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| charWidth        | *int*  | **-1**           |                  |
+------------------+--------+------------------+------------------+
| minLength        | *int*  | **None**         | Minimum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| maxLength        | *int*  | **None**         | Maximum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| onlyEnumerators  | *bool* | **False**        | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **[]**           | A list of        |
|                  |        |                  | allowed or       |
|                  |        |                  | recommended      |
|                  |        |                  | values for       |
|                  |        |                  | string           |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | equivalent to    |
|                  |        |                  | the enumerators  |
|                  |        |                  | that will appear |
|                  |        |                  | in the GUI       |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+
| allowAlias       | *bool* | **True**         | Flag if allow    |
|                  |        |                  | project to be    |
|                  |        |                  | directory alias  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| allowUnfound     | *bool* | **True**         | Flag if allow    |
|                  |        |                  | unfound project  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+

CSceneDataFile: An xml format file for defining scene in CCP4mg.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CSceneDataFile <./CCP4File.html#CSceneDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-15:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-15:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **CCP4mg scene** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | *                | A list of        |
|                  |        | *['scene.xml']** | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/CCP4-scene** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **CCP4mg scene   |                  |
| eTypeDescription |        | file**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **scene**        | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **               | A string         |
| ContentClassName |        | NotImplemented** | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CSearchPath
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CSearchPath <./CCP4File.html#CSearchPath>`__

.. _CCP4File-CCP4File-contents-of-class-16:

Contents of class:
^^^^^^^^^^^^^^^^^^

==== =========================================
path `CDataFile <./CCP4File.html#CDataFile>`__
name `CString <./CCP4Data.html#CString>`__
==== =========================================

.. _CCP4File-CCP4File-qualifiers-for-class-16:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CSearchPathList
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4File.CSearchPathList <./CCP4File.html#CSearchPathList>`__

.. _CCP4File-CCP4File-contents-of-class-17:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CSearchPath <./CCP4File.html#CSearchPath>`__

.. _CCP4File-CCP4File-qualifiers-for-class-17:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| listMinLength  | *int*  | **0**            | Inclusive        |
|                |        |                  | minimum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listMaxLength  | *int*  | **               | Inclusive        |
|                |        | NotImplemented** | maximum length   |
|                |        |                  | of list          |
+----------------+--------+------------------+------------------+
| listCompare    | *int*  | **               | If has value     |
|                |        | NotImplemented** | 1/-1 consecutive |
|                |        |                  | items in list    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than preceeding  |
|                |        |                  | item. The list   |
|                |        |                  | item class must  |
|                |        |                  | have a           |
|                |        |                  | \__cmp__()       |
|                |        |                  | method.          |
+----------------+--------+------------------+------------------+

CTextDataFile: A text data file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CTextDataFile <./CCP4File.html#CTextDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-18:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-18:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['txt',        | A list of        |
|                  |        | 'log']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **"text/plain"** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Standard plain |                  |
| eTypeDescription |        | text**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CVersion: A (string) version number of the form n.m.i
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4File.CVersion <./CCP4File.html#CVersion>`__

.. _CCP4File-CCP4File-contents-of-class-19:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4File-CCP4File-qualifiers-for-class-19:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if allow an |
|                  |        |                  | unset value at   |
|                  |        |                  | run time         |
+------------------+--------+------------------+------------------+
| default          | *str*  | **None**         | A default value  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| charWidth        | *int*  | **10**           | Number of        |
|                  |        |                  | characters       |
|                  |        |                  | allowed for      |
|                  |        |                  | widget in GUI    |
+------------------+--------+------------------+------------------+
| minLength        | *int*  | **None**         | Minimum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| maxLength        | *int*  | **None**         | Maximum length   |
|                  |        |                  | of string        |
+------------------+--------+------------------+------------------+
| onlyEnumerators  | *bool* | **False**        | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **[]**           | A list of        |
|                  |        |                  | allowed or       |
|                  |        |                  | recommended      |
|                  |        |                  | values for       |
|                  |        |                  | string           |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **[]**           | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | equivalent to    |
|                  |        |                  | the enumerators  |
|                  |        |                  | that will appear |
|                  |        |                  | in the GUI       |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+

CXmgrDataFile: An xmgr format file. This is the input format for xmgrace, as output by scala or aimless
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CXmgrDataFile <./CCP4File.html#CXmgrDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-20:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-20:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['xmgr']**     | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  |                  |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+

CXmlDataFile: A reference to an XML file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CXmlDataFile <./CCP4File.html#CXmlDataFile>`__

.. _CCP4File-CCP4File-contents-of-class-21:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ===========================================

.. _CCP4File-CCP4File-qualifiers-for-class-21:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **False**        | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['xml']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  |                  |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  |                  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **None**         | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **False**        | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **True**         | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **None**         | A string         |
| ContentClassName |        |                  | containing the   |
|                  |        |                  | name of a class  |
|                  |        |                  | which will hold  |
|                  |        |                  | the file         |
|                  |        |                  | contents         |
+------------------+--------+------------------+------------------+
| isDirectory      | *bool* | **False**        | Flag if the data |
|                  |        |                  | is a directory   |
+------------------+--------+------------------+------------------+
| requiredSubType  | *list* | **None**         | A list of        |
|                  |        |                  | allowed sub      |
|                  |        |                  | types            |
+------------------+--------+------------------+------------------+
| req              | *list* | **None**         | A list of        |
| uiredContentFlag |        |                  | allowed content  |
|                  |        |                  | flags            |
+------------------+--------+------------------+------------------+
