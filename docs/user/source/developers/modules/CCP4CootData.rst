CCootHistoryDataFile
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4CootData.CCootHistoryDataFile <./CCP4CootData.html#CCootHistoryDataFile>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
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
| toolTip          | *str*  | **history.scm or |                  |
|                  |        | 0-state.scm file |                  |
|                  |        | from Coot**      |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Coot history** |                  |
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
| fileExtensions   | *list* | **['scm',        | A list of        |
|                  |        | 'py']**          | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicati      |                  |
|                  |        | on/coot-script** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Coot           |                  |
| eTypeDescription |        | history/script   |                  |
|                  |        | file**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **coot_history** | Label for file   |
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
