CAtomSelection
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CAtomSelection <./CCP4ModelData.html#CAtomSelection>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

==== =====================================
text `CString <./CCP4Data.html#CString>`__
==== =====================================

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

CBlastData
~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4ModelData.CBlastData <./CCP4ModelData.html#CBlastData>`__

.. _CCP4ModelData-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

============= =====================================
queryId       `CString <./CCP4Data.html#CString>`__
alignmentList `CList <./CCP4Data.html#CList>`__
============= =====================================

.. _CCP4ModelData-qualifiers-for-class-1:

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

CBlastDataFile
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CBlastDataFile <./CCP4ModelData.html#CBlastDataFile>`__

.. _CCP4ModelData-contents-of-class-2:

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

.. _CCP4ModelData-qualifiers-for-class-2:

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
| guiLabel         | *str*  | **Blast          |                  |
|                  |        | results**        |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | model_data#ali** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['bla',        | A list of        |
|                  |        | 'blast',         | strings          |
|                  |        | 'xml']**         | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **application/Bl |                  |
|                  |        | ast-alignments** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Blast sequence |                  |
| eTypeDescription |        | search results** |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **Blast sequence | Label for file   |
|                  |        | search**         |                  |
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
| file             | *str*  | **CBlastData**   | A string         |
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

CBlastItem
~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CBlastItem <./CCP4ModelData.html#CBlastItem>`__

.. _CCP4ModelData-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

============= =====================================
hitSequence   `CString <./CCP4Data.html#CString>`__
hitId         `CString <./CCP4Data.html#CString>`__
querySequence `CString <./CCP4Data.html#CString>`__
============= =====================================

.. _CCP4ModelData-qualifiers-for-class-3:

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

CChemComp
~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CChemComp <./CCP4ModelData.html#CChemComp>`__

.. _CCP4ModelData-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

================= =======================================
number_atoms_all  `CInt <./CCP4Data.html#CInt>`__
group             `CString <./CCP4Data.html#CString>`__
name              `CString <./CCP4Data.html#CString>`__
desc_level        `CInt <./CCP4Data.html#CInt>`__
number_atoms_nh   `CInt <./CCP4Data.html#CInt>`__
three_letter_code `COneWord <./CCP4Data.html#COneWord>`__
id                `COneWord <./CCP4Data.html#COneWord>`__
================= =======================================

.. _CCP4ModelData-qualifiers-for-class-4:

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

CContainsSeMet
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CBoolean <./CCP4Data.html#CBoolean>`__ ->
`CCP4ModelData.CContainsSeMet <./CCP4ModelData.html#CContainsSeMet>`__

.. _CCP4ModelData-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4ModelData-qualifiers-for-class-5:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *bool* | **False**        |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Does structure |                  |
|                |        | contain          |                  |
|                |        | sel              |                  |
|                |        | enomethionine?** |                  |
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
| charWidth      | *int*  | **10**           |                  |
+----------------+--------+------------------+------------------+

CDictData
~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CDictData <./CCP4ModelData.html#CDictData>`__

.. _CCP4ModelData-contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== =================================
monomerList `CList <./CCP4Data.html#CList>`__
=========== =================================

.. _CCP4ModelData-qualifiers-for-class-6:

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

CDictDataFile: A refmac dictionary file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CDictDataFile <./CCP4ModelData.html#CDictDataFile>`__

.. _CCP4ModelData-contents-of-class-7:

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

.. _CCP4ModelData-qualifiers-for-class-7:

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
| toolTip          | *str*  | **Idealised      |                  |
|                  |        | geometry of      |                  |
|                  |        | ligands for      |                  |
|                  |        | refinement**     |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Geometry       |                  |
|                  |        | dictionary**     |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_data#l   |                  |
|                  |        | igand_geometry** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['cif']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | *                |                  |
|                  |        | *application/ref |                  |
|                  |        | mac-dictionary** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Geometry       |                  |
| eTypeDescription |        | file**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **dictionary**   | Label for file   |
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
| file             | *str*  | **CDictData**    | A string         |
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

CElement: Chemical element
~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.COneWord <./CCP4Data.html#COneWord>`__ ->
`CCP4ModelData.CElement <./CCP4ModelData.html#CElement>`__

.. _CCP4ModelData-contents-of-class-8:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4ModelData-qualifiers-for-class-8:

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
| onlyEnumerators  | *bool* | **True**         | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **['H', 'He',    | A list of        |
|                  |        | 'Li', 'Be', 'B', | allowed or       |
|                  |        | 'C', 'N', 'O',   | recommended      |
|                  |        | 'F', 'Ne', 'Na', | values for       |
|                  |        | 'Mg', 'Al',      | string           |
|                  |        | 'Si', 'P', 'S',  |                  |
|                  |        | 'Cl', 'Ar', 'K', |                  |
|                  |        | 'Ca', 'Sc',      |                  |
|                  |        | 'Ti', 'V', 'Cr', |                  |
|                  |        | 'Mn', 'Fe',      |                  |
|                  |        | 'Co', 'Ni',      |                  |
|                  |        | 'Cu', 'Zn',      |                  |
|                  |        | 'Ga', 'Ge',      |                  |
|                  |        | 'As', 'Se',      |                  |
|                  |        | 'Br', 'Kr',      |                  |
|                  |        | 'Rb', 'Sr', 'Y', |                  |
|                  |        | 'Zr', 'Nb',      |                  |
|                  |        | 'Mo', 'Tc',      |                  |
|                  |        | 'Ru', 'Rh',      |                  |
|                  |        | 'Pd', 'Ag',      |                  |
|                  |        | 'Cd', 'In',      |                  |
|                  |        | 'Sn', 'Sb',      |                  |
|                  |        | 'Te', 'I', 'Xe', |                  |
|                  |        | 'Cs', 'Ba',      |                  |
|                  |        | 'La', 'Ce',      |                  |
|                  |        | 'Pr', 'Nd',      |                  |
|                  |        | 'Pm', 'Sm',      |                  |
|                  |        | 'Eu', 'Gd',      |                  |
|                  |        | 'Tb', 'Dy',      |                  |
|                  |        | 'Ho', 'Er',      |                  |
|                  |        | 'Tm', 'Yb',      |                  |
|                  |        | 'Lu', 'Hf',      |                  |
|                  |        | 'Ta', 'W', 'Re', |                  |
|                  |        | 'Os', 'Ir',      |                  |
|                  |        | 'Pt', 'Au',      |                  |
|                  |        | 'Hg', 'Tl',      |                  |
|                  |        | 'Pb', 'Bi',      |                  |
|                  |        | 'Po', 'At',      |                  |
|                  |        | 'Rn', 'Fr',      |                  |
|                  |        | 'Ra', 'Ac',      |                  |
|                  |        | 'Th', 'Pa', 'U', |                  |
|                  |        | 'Np', 'Pu',      |                  |
|                  |        | 'Am', 'Cm',      |                  |
|                  |        | 'Bk', 'Cf',      |                  |
|                  |        | 'Es', 'Fm',      |                  |
|                  |        | 'Md', 'No',      |                  |
|                  |        | 'Lr', 'Rf',      |                  |
|                  |        | 'Db', 'Sg',      |                  |
|                  |        | 'Bh', 'Hs',      |                  |
|                  |        | 'Mt', 'Ds',      |                  |
|                  |        | 'Rg', 'Cn']**    |                  |
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

CEnsemble: An ensemble of models. Typically, this would be a set of related PDB files, but models could also be xtal or EM maps. This should be indicated by the types entry. A single ensemble is a CList of structures.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CEnsemble <./CCP4ModelData.html#CEnsemble>`__

.. _CCP4ModelData-contents-of-class-9:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== =======================================
use         `CBoolean <./CCP4Data.html#CBoolean>`__
pdbItemList `CList <./CCP4Data.html#CList>`__
number      `CInt <./CCP4Data.html#CInt>`__
label       `COneWord <./CCP4Data.html#COneWord>`__
=========== =======================================

.. _CCP4ModelData-qualifiers-for-class-9:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **False**          
default        *dict* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **Ensemble**       
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CEnsembleList
~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4ModelData.CEnsembleList <./CCP4ModelData.html#CEnsembleList>`__

.. _CCP4ModelData-contents-of-class-10:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CEnsemble <./CCP4ModelData.html#CEnsemble>`__

.. _CCP4ModelData-qualifiers-for-class-10:

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

CEnsemblePdbDataFile: A PDB coordinate file containing ensemble of structures as 'NMR' models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__ ->
`CCP4ModelData.CEnsemblePdbDataFile <./CCP4ModelData.html#CEnsemblePdbDataFile>`__

.. _CCP4ModelData-contents-of-class-11:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ========================================================
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
selection   `CAtomSelection <./CCP4ModelData.html#CAtomSelection>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ========================================================

.. _CCP4ModelData-qualifiers-for-class-11:

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
| toolTip          | *str*  | **An ensemble of |                  |
|                  |        | model            |                  |
|                  |        | coordinates in   |                  |
|                  |        | PDB or mmCIF     |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Model          |                  |
|                  |        | ensemble**       |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_         |                  |
|                  |        | data#ensemble_co |                  |
|                  |        | ordinate_files** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['pdb', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **               |                  |
|                  |        | chemical/x-pdb** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Model          |                  |
| eTypeDescription |        | coordinates**    |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **ensemble       | Label for file   |
|                  |        | coordinates**    |                  |
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
| file             | *str*  | **CPdbData**     | A string         |
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

CHhpredData
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4ModelData.CHhpredData <./CCP4ModelData.html#CHhpredData>`__

.. _CCP4ModelData-contents-of-class-12:

Contents of class:
^^^^^^^^^^^^^^^^^^

============= =================================
alignmentList `CList <./CCP4Data.html#CList>`__
============= =================================

.. _CCP4ModelData-qualifiers-for-class-12:

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

CHhpredDataFile
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CHhpredDataFile <./CCP4ModelData.html#CHhpredDataFile>`__

.. _CCP4ModelData-contents-of-class-13:

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

.. _CCP4ModelData-qualifiers-for-class-13:

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
| guiLabel         | *str*  | **HHPred         |                  |
|                  |        | results**        |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | model_data#ali** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['hhr']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | *                |                  |
|                  |        | *application/HHP |                  |
|                  |        | red-alignments** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **HHPred         |                  |
| eTypeDescription |        | sequence search  |                  |
|                  |        | results**        |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **HHPred         | Label for file   |
|                  |        | sequence         |                  |
|                  |        | search**         |                  |
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
| file             | *str*  | **CHhpredData**  | A string         |
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

CHhpredItem
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CHhpredItem <./CCP4ModelData.html#CHhpredItem>`__

.. _CCP4ModelData-contents-of-class-14:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =====================================
identifier `CString <./CCP4Data.html#CString>`__
annotation `CString <./CCP4Data.html#CString>`__
chain      `CString <./CCP4Data.html#CString>`__
========== =====================================

.. _CCP4ModelData-qualifiers-for-class-14:

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

CHomolog
~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__ ->
`CCP4ModelData.CHomolog <./CCP4ModelData.html#CHomolog>`__

.. _CCP4ModelData-contents-of-class-15:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ========================================================
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
selection   `CAtomSelection <./CCP4ModelData.html#CAtomSelection>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ========================================================

.. _CCP4ModelData-qualifiers-for-class-15:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **False**        | Flag if data     |
|                  |        |                  | file can be      |
|                  |        |                  | undefined at run |
|                  |        |                  | time             |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **A model        |                  |
|                  |        | coordinate file  |                  |
|                  |        | in PDB or mmCIF  |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Atomic model** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_data#co  |                  |
|                  |        | ordinate_files** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['pdb', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **               |                  |
|                  |        | chemical/x-pdb** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Model          |                  |
| eTypeDescription |        | coordinates**    |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **coordinates**  | Label for file   |
+------------------+--------+------------------+------------------+
| mustExist        | *bool* | **True**         | Flag if data     |
|                  |        |                  | file must exist  |
|                  |        |                  | at run time      |
+------------------+--------+------------------+------------------+
| fromPreviousJob  | *bool* | **False**        | Flag if input    |
|                  |        |                  | data file can be |
|                  |        |                  | inferred from    |
|                  |        |                  | preceeding jobs  |
+------------------+--------+------------------+------------------+
| jobCombo         | *bool* | **False**        | Flag if data     |
|                  |        |                  | widget should be |
|                  |        |                  | a combo box      |
+------------------+--------+------------------+------------------+
| file             | *str*  | **CPdbData**     | A string         |
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

CHomologList
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4ModelData.CHomologList <./CCP4ModelData.html#CHomologList>`__

.. _CCP4ModelData-contents-of-class-16:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CHomolog <./CCP4ModelData.html#CHomolog>`__

.. _CCP4ModelData-qualifiers-for-class-16:

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

CMDLMolDataFile: A molecule definition file (MDL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CMDLMolDataFile <./CCP4ModelData.html#CMDLMolDataFile>`__

.. _CCP4ModelData-contents-of-class-17:

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

.. _CCP4ModelData-qualifiers-for-class-17:

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
| toolTip          | *str*  | **Structure      |                  |
|                  |        | geometry of      |                  |
|                  |        | ligands for      |                  |
|                  |        | refinement in    |                  |
|                  |        | MDL mol format** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Mol file**     |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model          |                  |
|                  |        | _data#mol_file** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mol']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **chemical       |                  |
|                  |        | /x-mdl-molfile** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MDL Molfile**  |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **mol**          | Label for file   |
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

CMonomer: A monomer compound. ?smiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CMonomer <./CCP4ModelData.html#CMonomer>`__

.. _CCP4ModelData-contents-of-class-18:

Contents of class:
^^^^^^^^^^^^^^^^^^

============== =====================================
formula        `CString <./CCP4Data.html#CString>`__
identifier     `CString <./CCP4Data.html#CString>`__
dictionaryName `CString <./CCP4Data.html#CString>`__
smiles         `CString <./CCP4Data.html#CString>`__
============== =====================================

.. _CCP4ModelData-qualifiers-for-class-18:

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

CPdbData: Contents of a PDB file - a subset with functionality for GUI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4ModelData.CPdbData <./CCP4ModelData.html#CPdbData>`__

.. _CCP4ModelData-contents-of-class-19:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4ModelData-qualifiers-for-class-19:

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

CPdbDataFile: A PDB coordinate file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__

.. _CCP4ModelData-contents-of-class-20:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ========================================================
project     `CProjectId <./CCP4File.html#CProjectId>`__
subType     `CInt <./CCP4Data.html#CInt>`__
dbFileId    `CUUID <./CCP4Data.html#CUUID>`__
contentFlag `CInt <./CCP4Data.html#CInt>`__
selection   `CAtomSelection <./CCP4ModelData.html#CAtomSelection>`__
baseName    `CFilePath <./CCP4File.html#CFilePath>`__
annotation  `CString <./CCP4Data.html#CString>`__
relPath     `CFilePath <./CCP4File.html#CFilePath>`__
=========== ========================================================

.. _CCP4ModelData-qualifiers-for-class-20:

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
| toolTip          | *str*  | **A model        |                  |
|                  |        | coordinate file  |                  |
|                  |        | in PDB or mmCIF  |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Atomic model** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_data#co  |                  |
|                  |        | ordinate_files** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['pdb', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **               |                  |
|                  |        | chemical/x-pdb** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Model          |                  |
| eTypeDescription |        | coordinates**    |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **coordinates**  | Label for file   |
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
| file             | *str*  | **CPdbData**     | A string         |
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

CPdbDataFileList
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4ModelData.CPdbDataFileList <./CCP4ModelData.html#CPdbDataFileList>`__

.. _CCP4ModelData-contents-of-class-21:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__

.. _CCP4ModelData-qualifiers-for-class-21:

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

CPdbEnsembleItem
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CPdbEnsembleItem <./CCP4ModelData.html#CPdbEnsembleItem>`__

.. _CCP4ModelData-contents-of-class-22:

Contents of class:
^^^^^^^^^^^^^^^^^^

================== ====================================================
identity_to_target `CFloat <./CCP4Data.html#CFloat>`__
structure          `CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__
rms_to_target      `CFloat <./CCP4Data.html#CFloat>`__
================== ====================================================

.. _CCP4ModelData-qualifiers-for-class-22:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **False**        |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Homologous     |                  |
|                |        | model and its    |                  |
|                |        | similarity to    |                  |
|                |        | the target       |                  |
|                |        | structure**      |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Structure in   |                  |
|                |        | ensemble**       |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+

CRefmacNcs: Definition of a NCS for Refmac - not implemented yet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CRefmacNcs <./CCP4ModelData.html#CRefmacNcs>`__

.. _CCP4ModelData-contents-of-class-23:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4ModelData-qualifiers-for-class-23:

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

CRefmacRigidDomain: Definition of a rigid domain for Refmac
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CRefmacRigidDomain <./CCP4ModelData.html#CRefmacRigidDomain>`__

.. _CCP4ModelData-contents-of-class-24:

Contents of class:
^^^^^^^^^^^^^^^^^^

+-------------------+-------------------------------------------------+
| excludeSideChains | `CBoolean <./CCP4Data.html#CBoolean>`__         |
+-------------------+-------------------------------------------------+
| selection         | `CResidueRangeL                                 |
|                   | ist <./CCP4ModelData.html#CResidueRangeList>`__ |
+-------------------+-------------------------------------------------+
| applyRotation     | `CBoolean <./CCP4Data.html#CBoolean>`__         |
+-------------------+-------------------------------------------------+
| transformation    | `CTransform                                     |
|                   | ation <./CCP4MathsData.html#CTransformation>`__ |
+-------------------+-------------------------------------------------+
| applyTranslation  | `CBoolean <./CCP4Data.html#CBoolean>`__         |
+-------------------+-------------------------------------------------+

.. _CCP4ModelData-qualifiers-for-class-24:

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

CResidueRange: A residue range selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CResidueRange <./CCP4ModelData.html#CResidueRange>`__

.. _CCP4ModelData-contents-of-class-25:

Contents of class:
^^^^^^^^^^^^^^^^^^

======== =======================================
firstRes `COneWord <./CCP4Data.html#COneWord>`__
chainId  `COneWord <./CCP4Data.html#COneWord>`__
lastRes  `COneWord <./CCP4Data.html#COneWord>`__
======== =======================================

.. _CCP4ModelData-qualifiers-for-class-25:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
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
| pdbFileKey     | *str*  | **None**         | The key for a    |
|                |        |                  | CPdbDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+

CResidueRangeList: A list of residue range selections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4ModelData.CResidueRangeList <./CCP4ModelData.html#CResidueRangeList>`__

.. _CCP4ModelData-contents-of-class-26:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CResidueRange <./CCP4ModelData.html#CResidueRange>`__

.. _CCP4ModelData-qualifiers-for-class-26:

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

CSeqAlignDataFile: A (multiple) sequence alignment file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CSeqAlignDataFile <./CCP4ModelData.html#CSeqAlignDataFile>`__

.. _CCP4ModelData-contents-of-class-27:

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

.. _CCP4ModelData-qualifiers-for-class-27:

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
| toolTip          | *str*  | **Multiple       |                  |
|                  |        | sequence         |                  |
|                  |        | alignment in any |                  |
|                  |        | of the common    |                  |
|                  |        | formats          |                  |
|                  |        | (pir,fasta..)**  |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Aligned        |                  |
|                  |        | sequences**      |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_d        |                  |
|                  |        | ata#alignments** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['aln', 'pir', | A list of        |
|                  |        | 'fasta', 'msf',  | strings          |
|                  |        | 'phy']**         | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **application    |                  |
|                  |        | /CCP4-seqalign** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Sequence       |                  |
| eTypeDescription |        | alignment file** |                  |
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
| file             | *str*  | **CSeq           | A string         |
| ContentClassName |        | uenceAlignment** | containing the   |
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
| r                | *list* | **               | A list of        |
| equiredSequences |        | NotImplemented** | allowed numbers  |
|                  |        |                  | of sequences in  |
|                  |        |                  | file (usually    |
|                  |        |                  | [2])             |
+------------------+--------+------------------+------------------+

CSeqDataFile: A sequence file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CSeqDataFile <./CCP4ModelData.html#CSeqDataFile>`__

.. _CCP4ModelData-contents-of-class-28:

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

.. _CCP4ModelData-qualifiers-for-class-28:

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
| guiLabel         | *str*  | **Sequence**     |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model_         |                  |
|                  |        | data#sequences** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['seq', 'pir', | A list of        |
|                  |        | 'fasta']**       | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applic         |                  |
|                  |        | ation/CCP4-seq** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Sequence       |                  |
| eTypeDescription |        | file**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **sequence**     | Label for file   |
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
| file             | *str*  | **CSequence**    | A string         |
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

CSeqDataFileList
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4ModelData.CSeqDataFileList <./CCP4ModelData.html#CSeqDataFileList>`__

.. _CCP4ModelData-contents-of-class-29:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CSeqDataFile <./CCP4ModelData.html#CSeqDataFile>`__

.. _CCP4ModelData-qualifiers-for-class-29:

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

CSequence: A string of sequence one-letter codes Need to be able to parse common seq file formats Do we need to support alternative residues What about nucleic/polysach?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CSequence <./CCP4ModelData.html#CSequence>`__

.. _CCP4ModelData-contents-of-class-30:

Contents of class:
^^^^^^^^^^^^^^^^^^

============ =====================================
moleculeType `CString <./CCP4Data.html#CString>`__
referenceDb  `CString <./CCP4Data.html#CString>`__
identifier   `CString <./CCP4Data.html#CString>`__
reference    `CString <./CCP4Data.html#CString>`__
sequence     `CString <./CCP4Data.html#CString>`__
============ =====================================

.. _CCP4ModelData-qualifiers-for-class-30:

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

CSequenceAlignment: An alignment of two or more sequences. Each sequence is obviously related to class CSequence, but will also contain gaps relevant to the alignment. We could implement the contents as a list of CSequence objects? The alignment is typically formatted in a file as consecutive or interleaved sequences.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CSequenceAlignment <./CCP4ModelData.html#CSequenceAlignment>`__

.. _CCP4ModelData-contents-of-class-31:

Contents of class:
^^^^^^^^^^^^^^^^^^

============ =====================================
moleculeType `CString <./CCP4Data.html#CString>`__
identifier   `CString <./CCP4Data.html#CString>`__
============ =====================================

.. _CCP4ModelData-qualifiers-for-class-31:

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

CSequenceMeta
~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4ModelData.CSequenceMeta <./CCP4ModelData.html#CSequenceMeta>`__

.. _CCP4ModelData-contents-of-class-32:

Contents of class:
^^^^^^^^^^^^^^^^^^

================ =====================================
expressionSystem `CString <./CCP4Data.html#CString>`__
uniprotId        `CString <./CCP4Data.html#CString>`__
organism         `CString <./CCP4Data.html#CString>`__
================ =====================================

.. _CCP4ModelData-qualifiers-for-class-32:

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

CTLSDataFile: A refmac TLS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4ModelData.CTLSDataFile <./CCP4ModelData.html#CTLSDataFile>`__

.. _CCP4ModelData-contents-of-class-33:

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

.. _CCP4ModelData-qualifiers-for-class-33:

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
| toolTip          | *str*  | **Definition of  |                  |
|                  |        | model domains    |                  |
|                  |        | for TLS          |                  |
|                  |        | refinement**     |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **TLS            |                  |
|                  |        | coefficients**   |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **model          |                  |
|                  |        | _data#tls_file** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['tls']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/refmac-TLS** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Refmac TLS     |                  |
| eTypeDescription |        | file**           |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **tls**          | Label for file   |
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
