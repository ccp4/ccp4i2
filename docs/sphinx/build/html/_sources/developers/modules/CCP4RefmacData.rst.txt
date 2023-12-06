CRefmacAnomalousAtom
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4RefmacData.CRefmacAnomalousAtom <./CCP4RefmacData.html#CRefmacAnomalousAtom>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

======== =====================================
Fp       `CFloat <./CCP4Data.html#CFloat>`__
atomType `CString <./CCP4Data.html#CString>`__
Fpp      `CFloat <./CCP4Data.html#CFloat>`__
======== =====================================

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

CRefmacRestraintsDataFile
~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4RefmacData.CRefmacRestraintsDataFile <./CCP4RefmacData.html#CRefmacRestraintsDataFile>`__

.. _CCP4RefmacData-contents-of-class-1:

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

.. _CCP4RefmacData-qualifiers-for-class-1:

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
| guiLabel         | *str*  | **Additional     |                  |
|                  |        | restraints**     |                  |
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
| fileExtensions   | *list* | **['txt']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/refmac-exter |                  |
|                  |        | nal-restraints** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Refmac         |                  |
| eTypeDescription |        | external         |                  |
|                  |        | restraints**     |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **restraints**   | Label for file   |
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

CRefmacRigidGroupItem
~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4RefmacData.CRefmacRigidGroupItem <./CCP4RefmacData.html#CRefmacRigidGroupItem>`__

.. _CCP4RefmacData-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

============== =====================================
rigid_group_id `CString <./CCP4Data.html#CString>`__
segmentList    `CList <./CCP4Data.html#CList>`__
============== =====================================

.. _CCP4RefmacData-qualifiers-for-class-2:

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

CRefmacRigidGroupList
~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4RefmacData.CRefmacRigidGroupList <./CCP4RefmacData.html#CRefmacRigidGroupList>`__

.. _CCP4RefmacData-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of
`CRefmacRigidGroupItem <./CCP4RefmacData.html#CRefmacRigidGroupItem>`__

.. _CCP4RefmacData-qualifiers-for-class-3:

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

CRefmacRigidGroupSegment
~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4RefmacData.CRefmacRigidGroupSegment <./CCP4RefmacData.html#CRefmacRigidGroupSegment>`__

.. _CCP4RefmacData-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

========= =====================================
chain_id  `CString <./CCP4Data.html#CString>`__
residue_1 `CInt <./CCP4Data.html#CInt>`__
residue_2 `CInt <./CCP4Data.html#CInt>`__
========= =====================================

.. _CCP4RefmacData-qualifiers-for-class-4:

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
