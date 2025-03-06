CAltSpaceGroup
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__ ->
`CCP4XtalData.CAltSpaceGroup <./CCP4XtalData.html#CAltSpaceGroup>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | *                |                  |
|                  |        | *Hermann-Mauguin |                  |
|                  |        | space group      |                  |
|                  |        | name**           |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **crystal_da     |                  |
|                  |        | ta#space_group** |                  |
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

CAltSpaceGroupList
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CAltSpaceGroupList <./CCP4XtalData.html#CAltSpaceGroupList>`__

.. _CCP4XtalData-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CAltSpaceGroup <./CCP4XtalData.html#CAltSpaceGroup>`__

.. _CCP4XtalData-qualifiers-for-class-1:

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

CAnomalousColumnGroup: Selection of F/I and AnomF/I columns from MTZ. Expected to be part of ab initio phasing dataset ( CDataset)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CAnomalousColumnGroup <./CCP4XtalData.html#CAnomalousColumnGroup>`__

.. _CCP4XtalData-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-2:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
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
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **['The real     | Tooltips for     |
|                |        | part of the      | columns in group |
|                |        | experimental     |                  |
|                |        | structure        |                  |
|                |        | factors', 'The   |                  |
|                |        | anomalous part   |                  |
|                |        | of the           |                  |
|                |        | experimental     |                  |
|                |        | structure        |                  |
|                |        | factors']**      |                  |
+----------------+--------+------------------+------------------+

CAnomalousIntensityColumnGroup: Selection of I and AnomI columns from MTZ. Expected to be part of ab initio phasing dataset ( CDataset)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CAnomalousIntensityColumnGroup <./CCP4XtalData.html#CAnomalousIntensityColumnGroup>`__

.. _CCP4XtalData-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-3:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Intensity and  |                  |
|                |        | anomalous        |                  |
|                |        | intensity**      |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **['The real     | Tooltips for     |
|                |        | part of the      | columns in group |
|                |        | experimental     |                  |
|                |        | intensity', 'The |                  |
|                |        | anomalous part   |                  |
|                |        | of the           |                  |
|                |        | experimental     |                  |
|                |        | intensity']**    |                  |
+----------------+--------+------------------+------------------+

CAnomalousScatteringElement: Definition of a anomalous scattering element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.COneWord <./CCP4Data.html#COneWord>`__ ->
`CCP4ModelData.CElement <./CCP4ModelData.html#CElement>`__ ->
`CCP4XtalData.CAnomalousScatteringElement <./CCP4XtalData.html#CAnomalousScatteringElement>`__

.. _CCP4XtalData-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-4:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **Se**           |                  |
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
| charWidth        | *int*  | **4**            |                  |
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
| enumerators      | *list* | **['Br', 'Fe',   | A list of        |
|                  |        | 'Pt', 'Se']**    | allowed or       |
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

CAsuComponent: A component of the asymmetric unit. This is for use in MR, defining what we are searching for. There are similarities to CCrystalComponents and it should maybe be merged.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CAsuComponent <./CCP4XtalData.html#CAsuComponent>`__

.. _CCP4XtalData-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

============== ====================================================
moleculeType   `CString <./CCP4Data.html#CString>`__
numberOfCopies `CInt <./CCP4Data.html#CInt>`__
seqFile        `CSeqDataFile <./CCP4ModelData.html#CSeqDataFile>`__
============== ====================================================

.. _CCP4XtalData-qualifiers-for-class-5:

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

CAsuComponentList
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CAsuComponentList <./CCP4XtalData.html#CAsuComponentList>`__

.. _CCP4XtalData-contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CAsuComponent <./CCP4XtalData.html#CAsuComponent>`__

.. _CCP4XtalData-qualifiers-for-class-6:

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
| guiLabel       | *str*  | **Contents of    |                  |
|                |        | asymmetric       |                  |
|                |        | unit**           |                  |
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

CAtomicFormFactors: Table of form factors for element v wavelength
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CTable <./CCP4Data.html#CTable>`__ ->
`CCP4XtalData.CAtomicFormFactors <./CCP4XtalData.html#CAtomicFormFactors>`__

.. _CCP4XtalData-contents-of-class-7:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-7:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *dict* | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Atomic form    |                  |
|                  |        | factors for      |                  |
|                  |        | given element    |                  |
|                  |        | and wavelength** |                  |
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
| allo             | *bool* | **False**        |                  |
| wUndefinedHeader |        |                  |                  |
+------------------+--------+------------------+------------------+
| uniqueHeader     | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| rowLabel         | *str*  | **atomic         |                  |
|                  |        | element**        |                  |
+------------------+--------+------------------+------------------+
| columnLabel      | *str*  | **wavelength**   |                  |
+------------------+--------+------------------+------------------+

CCell: A unit cell
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CCell <./CCP4XtalData.html#CCell>`__

.. _CCP4XtalData-contents-of-class-8:

Contents of class:
^^^^^^^^^^^^^^^^^^

===== =================================================
a     `CCellLength <./CCP4XtalData.html#CCellLength>`__
c     `CCellLength <./CCP4XtalData.html#CCellLength>`__
b     `CCellLength <./CCP4XtalData.html#CCellLength>`__
beta  `CCellAngle <./CCP4XtalData.html#CCellAngle>`__
alpha `CCellAngle <./CCP4XtalData.html#CCellAngle>`__
gamma `CCellAngle <./CCP4XtalData.html#CCellAngle>`__
===== =================================================

.. _CCP4XtalData-qualifiers-for-class-8:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Cell lengths   |                  |
|                |        | and angles**     |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **cry            |                  |
|                |        | stal_data#cell** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+

CCellAngle: A cell angle
~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CFloat <./CCP4Data.html#CFloat>`__ ->
`CCP4XtalData.CCellAngle <./CCP4XtalData.html#CCellAngle>`__

.. _CCP4XtalData-contents-of-class-9:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-9:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+---------+------------------+------------------+
| allowUndefined  | *bool*  | **True**         |                  |
+-----------------+---------+------------------+------------------+
| default         | *float* | **None**         |                  |
+-----------------+---------+------------------+------------------+
| toolTip         | *str*   | **Cell angle in  |                  |
|                 |         | degrees**        |                  |
+-----------------+---------+------------------+------------------+
| guiLabel        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| guiDefinition   | *dict*  | **{}**           |                  |
+-----------------+---------+------------------+------------------+
| helpFile        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| saveToDb        | *bool*  | **False**        | Save this data   |
|                 |         |                  | in the database  |
+-----------------+---------+------------------+------------------+
| charWidth       | *int*   | **10**           |                  |
+-----------------+---------+------------------+------------------+
| min             | *float* | **0.0**          | The inclusive    |
|                 |         |                  | minimum value    |
+-----------------+---------+------------------+------------------+
| max             | *float* | **180.0**        | The inclusive    |
|                 |         |                  | maximum value    |
+-----------------+---------+------------------+------------------+
| onlyEnumerators | *bool*  | **False**        | If this is true  |
|                 |         |                  | then the         |
|                 |         |                  | enumerators are  |
|                 |         |                  | obligatory -     |
|                 |         |                  | otherwise they   |
|                 |         |                  | are treated as   |
|                 |         |                  | recommended      |
|                 |         |                  | values           |
+-----------------+---------+------------------+------------------+
| enumerators     | *list*  | **[]**           | A Python list of |
|                 |         |                  | allowed or       |
|                 |         |                  | recommended      |
|                 |         |                  | values - see     |
|                 |         |                  | onlyEnumerators  |
+-----------------+---------+------------------+------------------+
| menuText        | *list*  | **[]**           | A Python list of |
|                 |         |                  | strings,         |
|                 |         |                  | matching items   |
|                 |         |                  | in enumerators   |
|                 |         |                  | list, to appear  |
|                 |         |                  | on GUI menu      |
+-----------------+---------+------------------+------------------+

CCellLength: A cell length
~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CFloat <./CCP4Data.html#CFloat>`__ ->
`CCP4XtalData.CCellLength <./CCP4XtalData.html#CCellLength>`__

.. _CCP4XtalData-contents-of-class-10:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-10:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+---------+------------------+------------------+
| allowUndefined  | *bool*  | **False**        |                  |
+-----------------+---------+------------------+------------------+
| default         | *float* | **None**         |                  |
+-----------------+---------+------------------+------------------+
| toolTip         | *str*   | **Cell length in |                  |
|                 |         | A**              |                  |
+-----------------+---------+------------------+------------------+
| guiLabel        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| guiDefinition   | *dict*  | **{}**           |                  |
+-----------------+---------+------------------+------------------+
| helpFile        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| saveToDb        | *bool*  | **False**        | Save this data   |
|                 |         |                  | in the database  |
+-----------------+---------+------------------+------------------+
| charWidth       | *int*   | **10**           |                  |
+-----------------+---------+------------------+------------------+
| min             | *float* | **0.0**          | The inclusive    |
|                 |         |                  | minimum value    |
+-----------------+---------+------------------+------------------+
| max             | *float* | **None**         | The inclusive    |
|                 |         |                  | maximum value    |
+-----------------+---------+------------------+------------------+
| onlyEnumerators | *bool*  | **False**        | If this is true  |
|                 |         |                  | then the         |
|                 |         |                  | enumerators are  |
|                 |         |                  | obligatory -     |
|                 |         |                  | otherwise they   |
|                 |         |                  | are treated as   |
|                 |         |                  | recommended      |
|                 |         |                  | values           |
+-----------------+---------+------------------+------------------+
| enumerators     | *list*  | **[]**           | A Python list of |
|                 |         |                  | allowed or       |
|                 |         |                  | recommended      |
|                 |         |                  | values - see     |
|                 |         |                  | onlyEnumerators  |
+-----------------+---------+------------------+------------------+
| menuText        | *list*  | **[]**           | A Python list of |
|                 |         |                  | strings,         |
|                 |         |                  | matching items   |
|                 |         |                  | in enumerators   |
|                 |         |                  | list, to appear  |
|                 |         |                  | on GUI menu      |
+-----------------+---------+------------------+------------------+

CColumnGroup: Groups of columns in MTZ - probably from analysis by hklfile
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CColumnGroup <./CCP4XtalData.html#CColumnGroup>`__

.. _CCP4XtalData-contents-of-class-11:

Contents of class:
^^^^^^^^^^^^^^^^^^

=============== =======================================
columnGroupType `COneWord <./CCP4Data.html#COneWord>`__
columnList      `CList <./CCP4Data.html#CList>`__
selected        `CBoolean <./CCP4Data.html#CBoolean>`__
contentFlag     `CInt <./CCP4Data.html#CInt>`__
dataset         `CString <./CCP4Data.html#CString>`__
=============== =======================================

.. _CCP4XtalData-qualifiers-for-class-11:

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

CColumnGroupItem: Definition of set of columns that form a 'group'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CColumnGroupItem <./CCP4XtalData.html#CColumnGroupItem>`__

.. _CCP4XtalData-contents-of-class-12:

Contents of class:
^^^^^^^^^^^^^^^^^^

============= =========================================================
defaultList   `CString <./CCP4Data.html#CString>`__
columnName    `COneWord <./CCP4Data.html#COneWord>`__
partnerTo     `COneWord <./CCP4Data.html#COneWord>`__
columnType    `CColumnTypeList <./CCP4XtalData.html#CColumnTypeList>`__
partnerOffset `CInt <./CCP4Data.html#CInt>`__
============= =========================================================

.. _CCP4XtalData-qualifiers-for-class-12:

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

CColumnGroupList
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CColumnGroupList <./CCP4XtalData.html#CColumnGroupList>`__

.. _CCP4XtalData-contents-of-class-13:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CColumnGroup <./CCP4XtalData.html#CColumnGroup>`__

.. _CCP4XtalData-qualifiers-for-class-13:

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

CColumnType: A list of recognised MTZ column types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CColumnType <./CCP4XtalData.html#CColumnType>`__

.. _CCP4XtalData-contents-of-class-14:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-14:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **F**            |                  |
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
| enumerators      | *list* | **['H', 'J',     | A list of        |
|                  |        | 'F', 'D', 'Q',   | allowed or       |
|                  |        | 'G', 'L', 'K',   | recommended      |
|                  |        | 'M', 'E', 'P',   | values for       |
|                  |        | 'W', 'A', 'B',   | string           |
|                  |        | 'Y', 'I', 'R']** |                  |
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

CColumnTypeList: A list of acceptable MTZ column types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CColumnTypeList <./CCP4XtalData.html#CColumnTypeList>`__

.. _CCP4XtalData-contents-of-class-15:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CColumnType <./CCP4XtalData.html#CColumnType>`__

.. _CCP4XtalData-qualifiers-for-class-15:

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

CCrystalComponents: A list of sequences, monomers and anomalous scatterers expected in a crystal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CCrystalComponents <./CCP4XtalData.html#CCrystalComponents>`__

.. _CCP4XtalData-contents-of-class-16:

Contents of class:
^^^^^^^^^^^^^^^^^^

============ =================================
monomers     `CList <./CCP4Data.html#CList>`__
phasingGroup `CList <./CCP4Data.html#CList>`__
sequences    `CList <./CCP4Data.html#CList>`__
============ =================================

.. _CCP4XtalData-qualifiers-for-class-16:

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

CCrystalCompositionLabel: Serves as column header for CCrystalComposition - is a name for a composition model The composition model can be for one of three units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CCrystalCompositionLabel <./CCP4XtalData.html#CCrystalCompositionLabel>`__

.. _CCP4XtalData-contents-of-class-17:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =====================================
identifier `CString <./CCP4Data.html#CString>`__
unit       `CString <./CCP4Data.html#CString>`__
========== =====================================

.. _CCP4XtalData-qualifiers-for-class-17:

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

CCrystalCompositionTable: A table of crystal components v. composition models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4XtalData.CCrystalCompositionTable <./CCP4XtalData.html#CCrystalCompositionTable>`__

.. _CCP4XtalData-contents-of-class-18:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

CCrystalName
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CCrystalName <./CCP4XtalData.html#CCrystalName>`__

.. _CCP4XtalData-contents-of-class-19:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-18:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **False**        |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Unique         |                  |
|                  |        | identifier for   |                  |
|                  |        | crystal (one     |                  |
|                  |        | word)**          |                  |
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
| minLength        | *int*  | **1**            | Minimum length   |
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

CDataset: The experimental data model for ab initio phasing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CDataset <./CCP4XtalData.html#CDataset>`__

.. _CCP4XtalData-contents-of-class-20:

Contents of class:
^^^^^^^^^^^^^^^^^^

================ ===================================================
datasetName      `CDatasetName <./CCP4XtalData.html#CDatasetName>`__
formFactorSource `CString <./CCP4Data.html#CString>`__
crystalName      `CCrystalName <./CCP4XtalData.html#CCrystalName>`__
selected         `CBoolean <./CCP4Data.html#CBoolean>`__
obsDataFile      `CObsDataFile <./CCP4XtalData.html#CObsDataFile>`__
formFactors      `CFormFactor <./CCP4XtalData.html#CFormFactor>`__
================ ===================================================

.. _CCP4XtalData-qualifiers-for-class-19:

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

CDatasetList
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CDatasetList <./CCP4XtalData.html#CDatasetList>`__

.. _CCP4XtalData-contents-of-class-21:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CDataset <./CCP4XtalData.html#CDataset>`__

.. _CCP4XtalData-qualifiers-for-class-20:

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

CDatasetName
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CDatasetName <./CCP4XtalData.html#CDatasetName>`__

.. _CCP4XtalData-contents-of-class-22:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-21:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **False**        |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Unique         |                  |
|                  |        | identifier for   |                  |
|                  |        | dataset (one     |                  |
|                  |        | word)**          |                  |
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
| minLength        | *int*  | **1**            | Minimum length   |
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

CExperimentalDataType: Experimental data type e.g. native or peak
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CExperimentalDataType <./CCP4XtalData.html#CExperimentalDataType>`__

.. _CCP4XtalData-contents-of-class-23:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-22:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **SAD**          |                  |
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
| enumerators      | *list* | **['native',     | A list of        |
|                  |        | 'derivative',    | allowed or       |
|                  |        | 'SAD', 'peak',   | recommended      |
|                  |        | 'inflection',    | values for       |
|                  |        | 'high_remote',   | string           |
|                  |        | 'low_remote',    |                  |
|                  |        | '']**            |                  |
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

CFPairColumnGroup
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CFPairColumnGroup <./CCP4XtalData.html#CFPairColumnGroup>`__

.. _CCP4XtalData-contents-of-class-24:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-23:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Anomalous      |                  |
|                |        | structure        |                  |
|                |        | factors and      |                  |
|                |        | sigma**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CFSigFColumnGroup
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CFSigFColumnGroup <./CCP4XtalData.html#CFSigFColumnGroup>`__

.. _CCP4XtalData-contents-of-class-25:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-24:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Structure      |                  |
|                |        | factor and       |                  |
|                |        | sigma**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CFormFactor: The for factor (Fp and Fpp) for a giving element and wavelength
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CFormFactor <./CCP4XtalData.html#CFormFactor>`__

.. _CCP4XtalData-contents-of-class-26:

Contents of class:
^^^^^^^^^^^^^^^^^^

=== ===================================
Fp  `CFloat <./CCP4Data.html#CFloat>`__
Fpp `CFloat <./CCP4Data.html#CFloat>`__
=== ===================================

.. _CCP4XtalData-qualifiers-for-class-25:

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

CFreeRColumnGroup
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CFreeRColumnGroup <./CCP4XtalData.html#CFreeRColumnGroup>`__

.. _CCP4XtalData-contents-of-class-27:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-26:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Set of FreeR   |                  |
|                |        | flags**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CFreeRDataFile
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
-> `CCP4XtalData.CFreeRDataFile <./CCP4XtalData.html#CFreeRDataFile>`__

.. _CCP4XtalData-contents-of-class-28:

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

.. _CCP4XtalData-qualifiers-for-class-27:

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
| toolTip          | *str*  | **Set of         |                  |
|                  |        | reflections used |                  |
|                  |        | for FreeR        |                  |
|                  |        | calculation**    |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Free R set**   |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **da             |                  |
|                  |        | ta_files#FreeR** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **               |                  |
|                  |        | application/CCP4 |                  |
|                  |        | -mtz-freerflag** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **FreeR flag**   |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **freeRflag**    | Label for file   |
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
| file             | *str*  | **CMtzData**     | A string         |
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
| correctColumns   | *list* | **['I']**        | A list of        |
|                  |        |                  | coloumn data     |
|                  |        |                  | types expected   |
|                  |        |                  | in the file      |
+------------------+--------+------------------+------------------+
| colu             | *str*  | **[]**           |                  |
| mnGroupClassList |        |                  |                  |
+------------------+--------+------------------+------------------+
| sameCrystalAs    | *str*  | **               | Name of          |
|                  |        | NotImplemented** | CMtzDataFile     |
|                  |        |                  | object that      |
|                  |        |                  | crystal          |
|                  |        |                  | parameters       |
|                  |        |                  | should match -   |
|                  |        |                  | probably the     |
|                  |        |                  | observed data    |
+------------------+--------+------------------+------------------+

CGenericReflDataFile
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CGenericReflDataFile <./CCP4XtalData.html#CGenericReflDataFile>`__

.. _CCP4XtalData-contents-of-class-29:

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

.. _CCP4XtalData-qualifiers-for-class-28:

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
| toolTip          | *str*  | **A reflection   |                  |
|                  |        | data file in MTZ |                  |
|                  |        | or a non-CCP4    |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Reflection     |                  |
|                  |        | data**           |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **import_merge   |                  |
|                  |        | d#file_formats** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'hkl', | A list of        |
|                  |        | 'HKL', 'sca',    | strings          |
|                  |        | 'SCA', 'mmcif',  | containing       |
|                  |        | 'cif', 'ent']**  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applic         |                  |
|                  |        | ation/CCP4-gener |                  |
|                  |        | ic-reflections** |                  |
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
| file             | *str*  | **CUnmer         | A string         |
| ContentClassName |        | gedDataContent** | containing the   |
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

CHLColumnGroup
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
-> `CCP4XtalData.CHLColumnGroup <./CCP4XtalData.html#CHLColumnGroup>`__

.. _CCP4XtalData-contents-of-class-30:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-29:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Hend           |                  |
|                |        | rickson-Lattmann |                  |
|                |        | coefficients**   |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CIPairColumnGroup
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CIPairColumnGroup <./CCP4XtalData.html#CIPairColumnGroup>`__

.. _CCP4XtalData-contents-of-class-31:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-30:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Anomalous      |                  |
|                |        | intensities and  |                  |
|                |        | sigma**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CISigIColumnGroup
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CISigIColumnGroup <./CCP4XtalData.html#CISigIColumnGroup>`__

.. _CCP4XtalData-contents-of-class-32:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-31:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Intensity and  |                  |
|                |        | sigma**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CImageFile
~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CImageFile <./CCP4XtalData.html#CImageFile>`__

.. _CCP4XtalData-contents-of-class-33:

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

.. _CCP4XtalData-qualifiers-for-class-32:

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
| toolTip          | *str*  | **First image    |                  |
|                  |        | file in a        |                  |
|                  |        | directory**      |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Image file**   |                  |
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
| fileExtensions   | *list* | **['img']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/CCP4-image** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Image file**   |                  |
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

CImosflmXmlDataFile: An iMosflm data file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CImosflmXmlDataFile <./CCP4XtalData.html#CImosflmXmlDataFile>`__

.. _CCP4XtalData-contents-of-class-34:

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

.. _CCP4XtalData-qualifiers-for-class-33:

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
| guiLabel         | *str*  | **iMosflm data** |                  |
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
| fileExtensions   | *list* | **[              | A list of        |
|                  |        | 'imosflm.xml']** | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicati      |                  |
|                  |        | on/iMosflm-xml** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **iMosflm data** |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **imosflm**      | Label for file   |
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

CImportUnmerged
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CImportUnmerged <./CCP4XtalData.html#CImportUnmerged>`__

.. _CCP4XtalData-contents-of-class-35:

Contents of class:
^^^^^^^^^^^^^^^^^^

+------------------+--------------------------------------------------+
| crystalName      | `CString <./CCP4Data.html#CString>`__            |
+------------------+--------------------------------------------------+
| dataset          | `CString <./CCP4Data.html#CString>`__            |
+------------------+--------------------------------------------------+
| cell             | `CCell <./CCP4XtalData.html#CCell>`__            |
+------------------+--------------------------------------------------+
| excludeSelection | `CRan                                            |
|                  | geSelection <./CCP4Data.html#CRangeSelection>`__ |
+------------------+--------------------------------------------------+
| file             | `CUnmergedDat                                    |
|                  | aFile <./CCP4XtalData.html#CUnmergedDataFile>`__ |
+------------------+--------------------------------------------------+
| wavelength       | `                                                |
|                  | CWavelength <./CCP4XtalData.html#CWavelength>`__ |
+------------------+--------------------------------------------------+

.. _CCP4XtalData-qualifiers-for-class-34:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Imported data  |                  |
|                |        | file, cell       |                  |
|                |        | parameters and   |                  |
|                |        | crystal/dataset  |                  |
|                |        | identifiers**    |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **import_merge   |                  |
|                |        | d#file_formats** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+

CImportUnmergedList
~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CImportUnmergedList <./CCP4XtalData.html#CImportUnmergedList>`__

.. _CCP4XtalData-contents-of-class-36:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CImportUnmerged <./CCP4XtalData.html#CImportUnmerged>`__

.. _CCP4XtalData-qualifiers-for-class-35:

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

CMapCoeffsDataFile
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
->
`CCP4XtalData.CMapCoeffsDataFile <./CCP4XtalData.html#CMapCoeffsDataFile>`__

.. _CCP4XtalData-contents-of-class-37:

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

.. _CCP4XtalData-qualifiers-for-class-36:

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
| toolTip          | *str*  | **Electron       |                  |
|                  |        | density map      |                  |
|                  |        | coefficients:    |                  |
|                  |        | F,Phi**          |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Map            |                  |
|                  |        | coefficients**   |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **data_f         |                  |
|                  |        | iles#MapCoeffs** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicatio     |                  |
|                  |        | n/CCP4-mtz-map** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ F-phi**    |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **ma             | Label for file   |
|                  |        | p_coefficients** |                  |
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
| file             | *str*  | **CMtzData**     | A string         |
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
| correctColumns   | *list* | **['FP',         | A list of        |
|                  |        | 'FQP']**         | coloumn data     |
|                  |        |                  | types expected   |
|                  |        |                  | in the file      |
+------------------+--------+------------------+------------------+
| colu             | *str*  | **[]**           |                  |
| mnGroupClassList |        |                  |                  |
+------------------+--------+------------------+------------------+
| sameCrystalAs    | *str*  | **               | Name of          |
|                  |        | NotImplemented** | CMtzDataFile     |
|                  |        |                  | object that      |
|                  |        |                  | crystal          |
|                  |        |                  | parameters       |
|                  |        |                  | should match -   |
|                  |        |                  | probably the     |
|                  |        |                  | observed data    |
+------------------+--------+------------------+------------------+

CMapColumnGroup
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CMapColumnGroup <./CCP4XtalData.html#CMapColumnGroup>`__

.. _CCP4XtalData-contents-of-class-38:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-37:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Structure      |                  |
|                |        | factor and phase |                  |
|                |        | to define a      |                  |
|                |        | map**            |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CMapDataFile: A CCP4 Map file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMapDataFile <./CCP4XtalData.html#CMapDataFile>`__

.. _CCP4XtalData-contents-of-class-39:

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

.. _CCP4XtalData-qualifiers-for-class-38:

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
| toolTip          | *str*  | **A map in CCP4  |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Electron       |                  |
|                  |        | Density Map**    |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **data_f         |                  |
|                  |        | iles#map_files** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['map']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applic         |                  |
|                  |        | ation/CCP4-map** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **CCP4 Electron  |                  |
| eTypeDescription |        | density map**    |                  |
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

CMergeMiniMtz
~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CMergeMiniMtz <./CCP4XtalData.html#CMergeMiniMtz>`__

.. _CCP4XtalData-contents-of-class-40:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===========================================================
columnNames `CString <./CCP4Data.html#CString>`__
columnTag   `CString <./CCP4Data.html#CString>`__
fileName    `CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
=========== ===========================================================

.. _CCP4XtalData-qualifiers-for-class-39:

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

CMergeMiniMtzList
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CMergeMiniMtzList <./CCP4XtalData.html#CMergeMiniMtzList>`__

.. _CCP4XtalData-contents-of-class-41:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CMergeMiniMtz <./CCP4XtalData.html#CMergeMiniMtz>`__

.. _CCP4XtalData-qualifiers-for-class-40:

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
| saveToDb       | *bool* | **True**         | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| listMinLength  | *int*  | **2**            | Inclusive        |
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

CMiniMtzDataFile
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__

.. _CCP4XtalData-contents-of-class-42:

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

.. _CCP4XtalData-qualifiers-for-class-41:

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
| toolTip          | *str*  | **Mini-MTZ file  |                  |
|                  |        | containing       |                  |
|                  |        | reflect          |                  |
|                  |        | ion,phases,FreeR |                  |
|                  |        | set or map       |                  |
|                  |        | coefficients**   |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Experimental   |                  |
|                  |        | data**           |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | data_files#MTZ** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **application    |                  |
|                  |        | /CCP4-mtz-mini** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ            |                  |
| eTypeDescription |        | experimental     |                  |
|                  |        | data**           |                  |
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
| file             | *str*  | **CMtzData**     | A string         |
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
| correctColumns   | *list* | **['FQ', 'JQ',   | A list of        |
|                  |        | 'GLGL', 'KMKM',  | coloumn data     |
|                  |        | 'AAAA', 'PW',    | types expected   |
|                  |        | 'FP', 'I']**     | in the file      |
+------------------+--------+------------------+------------------+
| colu             | *str*  | **               |                  |
| mnGroupClassList |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| sameCrystalAs    | *str*  | **               | Name of          |
|                  |        | NotImplemented** | CMtzDataFile     |
|                  |        |                  | object that      |
|                  |        |                  | crystal          |
|                  |        |                  | parameters       |
|                  |        |                  | should match -   |
|                  |        |                  | probably the     |
|                  |        |                  | observed data    |
+------------------+--------+------------------+------------------+

CMiniMtzDataFileList
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CMiniMtzDataFileList <./CCP4XtalData.html#CMiniMtzDataFileList>`__

.. _CCP4XtalData-contents-of-class-43:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__

.. _CCP4XtalData-qualifiers-for-class-42:

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

CMmcifReflData: Reflection data in mmCIF format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4File.CMmcifData <./CCP4File.html#CMmcifData>`__ ->
`CCP4XtalData.CMmcifReflData <./CCP4XtalData.html#CMmcifReflData>`__

.. _CCP4XtalData-contents-of-class-44:

Contents of class:
^^^^^^^^^^^^^^^^^^

================ =================================================
cell             `CCell <./CCP4XtalData.html#CCell>`__
haveIobsColumn   `CBoolean <./CCP4Data.html#CBoolean>`__
haveFobsColumn   `CBoolean <./CCP4Data.html#CBoolean>`__
wavelength       `CWavelength <./CCP4XtalData.html#CWavelength>`__
spaceGroup       `CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__
haveFpmObsColumn `CBoolean <./CCP4Data.html#CBoolean>`__
haveIpmObsColumn `CBoolean <./CCP4Data.html#CBoolean>`__
haveFreeRColumn  `CBoolean <./CCP4Data.html#CBoolean>`__
================ =================================================

.. _CCP4XtalData-qualifiers-for-class-43:

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

CMmcifReflDataFile: A reflection file in mmCIF format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4File.CMmcifDataFile <./CCP4File.html#CMmcifDataFile>`__ ->
`CCP4XtalData.CMmcifReflDataFile <./CCP4XtalData.html#CMmcifReflDataFile>`__

.. _CCP4XtalData-contents-of-class-45:

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

.. _CCP4XtalData-qualifiers-for-class-44:

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
| toolTip          | *str*  | **A reflection   |                  |
|                  |        | file in mmCIF    |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **mmCIF          |                  |
|                  |        | reflection       |                  |
|                  |        | data**           |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **da             |                  |
|                  |        | ta_files#mmCIF** |                  |
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
| mimeTypeName     | *str*  | **               |                  |
|                  |        | chemical/x-cif** |                  |
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
| file             | *str*  | **               | A string         |
| ContentClassName |        | CMmcifReflData** | containing the   |
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

CMtzColumn: An MTZ column with column label and column type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CMtzColumn <./CCP4XtalData.html#CMtzColumn>`__

.. _CCP4XtalData-contents-of-class-46:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== =================================================
groupIndex  `CInt <./CCP4Data.html#CInt>`__
columnLabel `COneWord <./CCP4Data.html#COneWord>`__
columnType  `CColumnType <./CCP4XtalData.html#CColumnType>`__
dataset     `COneWord <./CCP4Data.html#COneWord>`__
=========== =================================================

.. _CCP4XtalData-qualifiers-for-class-45:

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

CMtzColumnGroup
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CMtzColumnGroup <./CCP4XtalData.html#CMtzColumnGroup>`__

.. _CCP4XtalData-contents-of-class-47:

Contents of class:
^^^^^^^^^^^^^^^^^^

+-----------+---------------------------------------------------------+
| groupType | `CMtzColum                                              |
|           | nGroupType <./CCP4XtalData.html#CMtzColumnGroupType>`__ |
+-----------+---------------------------------------------------------+
| columns   | `CList <./CCP4Data.html#CList>`__                       |
+-----------+---------------------------------------------------------+

.. _CCP4XtalData-qualifiers-for-class-46:

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

CMtzColumnGroupType
~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CColumnType <./CCP4XtalData.html#CColumnType>`__ ->
`CCP4XtalData.CMtzColumnGroupType <./CCP4XtalData.html#CMtzColumnGroupType>`__

.. _CCP4XtalData-contents-of-class-48:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-47:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **F**            |                  |
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
| enumerators      | *list* | **['H', 'J',     | A list of        |
|                  |        | 'F', 'D', 'Q',   | allowed or       |
|                  |        | 'G', 'L', 'K',   | recommended      |
|                  |        | 'M', 'E', 'P',   | values for       |
|                  |        | 'W', 'A', 'B',   | string           |
|                  |        | 'Y', 'I', 'R']** |                  |
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

CMtzData: Some of the data contents of an MTZ file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4XtalData.CMtzData <./CCP4XtalData.html#CMtzData>`__

.. _CCP4XtalData-contents-of-class-49:

Contents of class:
^^^^^^^^^^^^^^^^^^

+-----------------+---------------------------------------------------+
| cell            | `CCell <./CCP4XtalData.html#CCell>`__             |
+-----------------+---------------------------------------------------+
| wavelengths     | `CList <./CCP4Data.html#CList>`__                 |
+-----------------+---------------------------------------------------+
| datasets        | `CList <./CCP4Data.html#CList>`__                 |
+-----------------+---------------------------------------------------+
| crystalNames    | `CList <./CCP4Data.html#CList>`__                 |
+-----------------+---------------------------------------------------+
| spaceGroup      | `CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__ |
+-----------------+---------------------------------------------------+
| resolutionRange | `CResoluti                                        |
|                 | onRange <./CCP4XtalData.html#CResolutionRange>`__ |
+-----------------+---------------------------------------------------+
| datasetCells    | `CList <./CCP4Data.html#CList>`__                 |
+-----------------+---------------------------------------------------+
| listOfColumns   | `CList <./CCP4Data.html#CList>`__                 |
+-----------------+---------------------------------------------------+
| merged          | `CBoolean <./CCP4Data.html#CBoolean>`__           |
+-----------------+---------------------------------------------------+

.. _CCP4XtalData-qualifiers-for-class-48:

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

CMtzDataFile: An MTZ experimental data file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__

.. _CCP4XtalData-contents-of-class-50:

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

.. _CCP4XtalData-qualifiers-for-class-49:

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
| toolTip          | *str*  | **Experimental   |                  |
|                  |        | data in CCP4's   |                  |
|                  |        | MTZ format**     |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Experimental   |                  |
|                  |        | data**           |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | data_files#MTZ** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applic         |                  |
|                  |        | ation/CCP4-mtz** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ            |                  |
| eTypeDescription |        | experimental     |                  |
|                  |        | data**           |                  |
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
| file             | *str*  | **CMtzData**     | A string         |
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

CMtzDataset
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CMtzDataset <./CCP4XtalData.html#CMtzDataset>`__

.. _CCP4XtalData-contents-of-class-51:

Contents of class:
^^^^^^^^^^^^^^^^^^

============ =====================================
columnGroups `CList <./CCP4Data.html#CList>`__
name         `CString <./CCP4Data.html#CString>`__
============ =====================================

.. _CCP4XtalData-qualifiers-for-class-50:

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

CObsDataFile
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
-> `CCP4XtalData.CObsDataFile <./CCP4XtalData.html#CObsDataFile>`__

.. _CCP4XtalData-contents-of-class-52:

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

.. _CCP4XtalData-qualifiers-for-class-51:

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
| toolTip          | *str*  | **Observed       |                  |
|                  |        | structure        |                  |
|                  |        | factors or       |                  |
|                  |        | intensities**    |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Reflections**  |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | data_files#Obs** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | *                |                  |
|                  |        | *application/CCP |                  |
|                  |        | 4-mtz-observed** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ observed** |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | *                | Label for file   |
|                  |        | *observed_data** |                  |
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
| file             | *str*  | **CMtzData**     | A string         |
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
| correctColumns   | *list* | **['KMKM',       | A list of        |
|                  |        | 'GLGL', 'JQ',    | coloumn data     |
|                  |        | 'FQ']**          | types expected   |
|                  |        |                  | in the file      |
+------------------+--------+------------------+------------------+
| colu             | *str*  | **[, , , ]**     |                  |
| mnGroupClassList |        |                  |                  |
+------------------+--------+------------------+------------------+
| sameCrystalAs    | *str*  | **               | Name of          |
|                  |        | NotImplemented** | CMtzDataFile     |
|                  |        |                  | object that      |
|                  |        |                  | crystal          |
|                  |        |                  | parameters       |
|                  |        |                  | should match -   |
|                  |        |                  | probably the     |
|                  |        |                  | observed data    |
+------------------+--------+------------------+------------------+

CPhaserSolDataFile
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CPhaserSolDataFile <./CCP4XtalData.html#CPhaserSolDataFile>`__

.. _CCP4XtalData-contents-of-class-53:

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

.. _CCP4XtalData-qualifiers-for-class-52:

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
| toolTip          | *str*  | **Possible       |                  |
|                  |        | solutions passed |                  |
|                  |        | between runs of  |                  |
|                  |        | the Phaser       |                  |
|                  |        | program**        |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Phaser         |                  |
|                  |        | solutions**      |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **data_f         |                  |
|                  |        | iles#phasersol** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['ph           | A list of        |
|                  |        | aser_sol.pkl']** | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/phaser-sol** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Phaser         |                  |
| eTypeDescription |        | solution file**  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **phaser_sol**   | Label for file   |
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

CPhasingGroup
~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CPhasingGroup <./CCP4XtalData.html#CPhasingGroup>`__

.. _CCP4XtalData-contents-of-class-54:

Contents of class:
^^^^^^^^^^^^^^^^^^

======= =====================================
formula `CString <./CCP4Data.html#CString>`__
type    `CString <./CCP4Data.html#CString>`__
======= =====================================

.. _CCP4XtalData-qualifiers-for-class-53:

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

CPhiFomColumnGroup
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__
->
`CCP4XtalData.CPhiFomColumnGroup <./CCP4XtalData.html#CPhiFomColumnGroup>`__

.. _CCP4XtalData-contents-of-class-55:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-54:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **Phase and      |                  |
|                |        | figure of        |                  |
|                |        | merit**          |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CPhsDataFile
~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
-> `CCP4XtalData.CPhsDataFile <./CCP4XtalData.html#CPhsDataFile>`__

.. _CCP4XtalData-contents-of-class-56:

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

.. _CCP4XtalData-qualifiers-for-class-55:

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
| toolTip          | *str*  | **Phases in      |                  |
|                  |        | Hend             |                  |
|                  |        | rickson-Lattmann |                  |
|                  |        | or Phi/FOM       |                  |
|                  |        | form**           |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Phases**       |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | data_files#Phs** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'cif', | A list of        |
|                  |        | 'ent']**         | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **application/C  |                  |
|                  |        | CP4-mtz-phases** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ phases**   |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **phases**       | Label for file   |
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
| file             | *str*  | **CMtzData**     | A string         |
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
| correctColumns   | *list* | **['AAAA',       | A list of        |
|                  |        | 'PW']**          | coloumn data     |
|                  |        |                  | types expected   |
|                  |        |                  | in the file      |
+------------------+--------+------------------+------------------+
| colu             | *str*  | **[, ]**         |                  |
| mnGroupClassList |        |                  |                  |
+------------------+--------+------------------+------------------+
| sameCrystalAs    | *str*  | **               | Name of          |
|                  |        | NotImplemented** | CMtzDataFile     |
|                  |        |                  | object that      |
|                  |        |                  | crystal          |
|                  |        |                  | parameters       |
|                  |        |                  | should match -   |
|                  |        |                  | probably the     |
|                  |        |                  | observed data    |
+------------------+--------+------------------+------------------+

CProgramColumnGroup: A group of MTZ columns required for program input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__

.. _CCP4XtalData-contents-of-class-57:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-56:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *list* | **[]**           | Preferred values |
|                |        |                  | for column names |
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
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+
| toolTipList    | *list* | **[]**           | Tooltips for     |
|                |        |                  | columns in group |
+----------------+--------+------------------+------------------+

CProgramColumnGroup0
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CProgramColumnGroup0 <./CCP4XtalData.html#CProgramColumnGroup0>`__

.. _CCP4XtalData-contents-of-class-58:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== =========================================================
datasetName `CString <./CCP4Data.html#CString>`__
columnGroup `CMtzColumnGroup <./CCP4XtalData.html#CMtzColumnGroup>`__
=========== =========================================================

.. _CCP4XtalData-qualifiers-for-class-57:

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
| groupTypes     | *list* | **[]**           | Type of          |
|                |        |                  | columnGroup      |
|                |        |                  | required by      |
|                |        |                  | program          |
+----------------+--------+------------------+------------------+
| mtzFileKey     | *str*  |                  | The key for a    |
|                |        |                  | CMtxDataFile in  |
|                |        |                  | the same         |
|                |        |                  | CContainer       |
+----------------+--------+------------------+------------------+
| mustExist      | *bool* | **False**        | Flag if the      |
|                |        |                  | parameter must   |
|                |        |                  | be set at run    |
|                |        |                  | time             |
+----------------+--------+------------------+------------------+

CReindexOperator
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CReindexOperator <./CCP4XtalData.html#CReindexOperator>`__

.. _CCP4XtalData-contents-of-class-59:

Contents of class:
^^^^^^^^^^^^^^^^^^

= =====================================
h `CString <./CCP4Data.html#CString>`__
k `CString <./CCP4Data.html#CString>`__
l `CString <./CCP4Data.html#CString>`__
= =====================================

.. _CCP4XtalData-qualifiers-for-class-58:

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

CResolutionRange
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CResolutionRange <./CCP4XtalData.html#CResolutionRange>`__

.. _CCP4XtalData-contents-of-class-60:

Contents of class:
^^^^^^^^^^^^^^^^^^

==== ===================================
high `CFloat <./CCP4Data.html#CFloat>`__
low  `CFloat <./CCP4Data.html#CFloat>`__
==== ===================================

.. _CCP4XtalData-qualifiers-for-class-59:

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

CRunBatchRange
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CRunBatchRange <./CCP4XtalData.html#CRunBatchRange>`__

.. _CCP4XtalData-contents-of-class-61:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ===============================
runNumber   `CInt <./CCP4Data.html#CInt>`__
fileNumber  `CInt <./CCP4Data.html#CInt>`__
batchRange0 `CInt <./CCP4Data.html#CInt>`__
batchRange1 `CInt <./CCP4Data.html#CInt>`__
=========== ===============================

.. _CCP4XtalData-qualifiers-for-class-60:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Specify range  |                  |
|                |        | of reflections   |                  |
|                |        | to treat as one  |                  |
|                |        | run**            |                  |
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

CRunBatchRangeList
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CRunBatchRangeList <./CCP4XtalData.html#CRunBatchRangeList>`__

.. _CCP4XtalData-contents-of-class-62:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CRunBatchRange <./CCP4XtalData.html#CRunBatchRange>`__

.. _CCP4XtalData-qualifiers-for-class-61:

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

CShelxFADataFile
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CShelxFADataFile <./CCP4XtalData.html#CShelxFADataFile>`__

.. _CCP4XtalData-contents-of-class-63:

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

.. _CCP4XtalData-qualifiers-for-class-62:

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
| toolTip          | *str*  | **Data used by   |                  |
|                  |        | Shelx programs** |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Shelx FA**     |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **data           |                  |
|                  |        | _files#shelxfa** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['hkl']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **application    |                  |
|                  |        | /CCP4-shelx-FA** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Shelx FA**     |                  |
| eTypeDescription |        |                  |                  |
+------------------+--------+------------------+------------------+
| fileLabel        | *str*  | **shelx_FA**     | Label for file   |
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

CShelxLabel
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CShelxLabel <./CCP4XtalData.html#CShelxLabel>`__

.. _CCP4XtalData-contents-of-class-64:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-63:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **UNDEFINED**    |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Hint to Shelx  |                  |
|                  |        | for the use of   |                  |
|                  |        | the dataset**    |                  |
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
| enumerators      | *list* | **['UNDEFINED',  | A list of        |
|                  |        | 'HREM', 'LREM',  | allowed or       |
|                  |        | 'PEAK', 'INFL',  | recommended      |
|                  |        | 'NAT', 'DERI']** | values for       |
|                  |        |                  | string           |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **['undefined',  | A list of        |
|                  |        | 'high remote',   | strings          |
|                  |        | 'low remote',    | equivalent to    |
|                  |        | 'peak',          | the enumerators  |
|                  |        | 'inflection',    | that will appear |
|                  |        | 'native',        | in the GUI       |
|                  |        | 'derivative']**  |                  |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+

CSpaceGroup: A string holding the space group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4XtalData.CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__

.. _CCP4XtalData-contents-of-class-65:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-64:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | *                |                  |
|                  |        | *Hermann-Mauguin |                  |
|                  |        | space group      |                  |
|                  |        | name**           |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **crystal_da     |                  |
|                  |        | ta#space_group** |                  |
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

CSpaceGroupCell: Cell space group and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4XtalData.CSpaceGroupCell <./CCP4XtalData.html#CSpaceGroupCell>`__

.. _CCP4XtalData-contents-of-class-66:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =================================================
cell       `CCell <./CCP4XtalData.html#CCell>`__
spaceGroup `CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__
========== =================================================

.. _CCP4XtalData-qualifiers-for-class-65:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Space group    |                  |
|                |        | and cell length  |                  |
|                |        | and angles**     |                  |
+----------------+--------+------------------+------------------+
| guiLabel       | *str*  | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| guiDefinition  | *dict* | **{}**           |                  |
+----------------+--------+------------------+------------------+
| helpFile       | *str*  | *                |                  |
|                |        | *crystal_data#ce |                  |
|                |        | ll_space_group** |                  |
+----------------+--------+------------------+------------------+
| saveToDb       | *bool* | **False**        | Save this data   |
|                |        |                  | in the database  |
+----------------+--------+------------------+------------------+

CUnmergedDataContent
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFileContent <./CCP4File.html#CDataFileContent>`__ ->
`CCP4XtalData.CUnmergedDataContent <./CCP4XtalData.html#CUnmergedDataContent>`__

.. _CCP4XtalData-contents-of-class-67:

Contents of class:
^^^^^^^^^^^^^^^^^^

=============== ===================================================
knownwavelength `CBoolean <./CCP4Data.html#CBoolean>`__
format          `CString <./CCP4Data.html#CString>`__
batchs          `CString <./CCP4Data.html#CString>`__
wavelength      `CWavelength <./CCP4XtalData.html#CWavelength>`__
spaceGroup      `CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__
knowncell       `CBoolean <./CCP4Data.html#CBoolean>`__
datasetName     `CDatasetName <./CCP4XtalData.html#CDatasetName>`__
numberLattices  `CInt <./CCP4Data.html#CInt>`__
crystalName     `CCrystalName <./CCP4XtalData.html#CCrystalName>`__
lowRes          `CFloat <./CCP4Data.html#CFloat>`__
cell            `CCell <./CCP4XtalData.html#CCell>`__
highRes         `CFloat <./CCP4Data.html#CFloat>`__
merged          `CString <./CCP4Data.html#CString>`__
=============== ===================================================

.. _CCP4XtalData-qualifiers-for-class-66:

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

CUnmergedDataFile: Handle MTZ, XDS and scalepack files. Allow wildcard filename
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CUnmergedDataFile <./CCP4XtalData.html#CUnmergedDataFile>`__

.. _CCP4XtalData-contents-of-class-68:

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

.. _CCP4XtalData-qualifiers-for-class-67:

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
| toolTip          | *str*  | **Unmerged       |                  |
|                  |        | experimental     |                  |
|                  |        | data in any      |                  |
|                  |        | format**         |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Unmerged       |                  |
|                  |        | reflections**    |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **data_files     |                  |
|                  |        | #unmerged_data** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz', 'hkl', | A list of        |
|                  |        | 'HKL', 'sca',    | strings          |
|                  |        | 'SCA']**         | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | **applicat       |                  |
|                  |        | ion/CCP4-unmerge |                  |
|                  |        | d-experimental** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **Unmerged       |                  |
| eTypeDescription |        | experimental     |                  |
|                  |        | data**           |                  |
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
| file             | *str*  | **CUnmer         | A string         |
| ContentClassName |        | gedDataContent** | containing the   |
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

CUnmergedDataFileList
~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4XtalData.CUnmergedDataFileList <./CCP4XtalData.html#CUnmergedDataFileList>`__

.. _CCP4XtalData-contents-of-class-69:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CUnmergedDataFile <./CCP4XtalData.html#CUnmergedDataFile>`__

.. _CCP4XtalData-qualifiers-for-class-68:

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

CUnmergedMtzDataFile
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4File.CDataFile <./CCP4File.html#CDataFile>`__ ->
`CCP4XtalData.CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ ->
`CCP4XtalData.CUnmergedMtzDataFile <./CCP4XtalData.html#CUnmergedMtzDataFile>`__

.. _CCP4XtalData-contents-of-class-70:

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

.. _CCP4XtalData-qualifiers-for-class-69:

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
| toolTip          | *str*  | **Unmerged       |                  |
|                  |        | experimental     |                  |
|                  |        | data in CCP4's   |                  |
|                  |        | MTZ format**     |                  |
+------------------+--------+------------------+------------------+
| guiLabel         | *str*  | **Unmerged MTZ   |                  |
|                  |        | reflections**    |                  |
+------------------+--------+------------------+------------------+
| guiDefinition    | *dict* | **{}**           |                  |
+------------------+--------+------------------+------------------+
| helpFile         | *str*  | **               |                  |
|                  |        | data_files#MTZ** |                  |
+------------------+--------+------------------+------------------+
| saveToDb         | *bool* | **True**         | Save the name of |
|                  |        |                  | this file in the |
|                  |        |                  | database         |
+------------------+--------+------------------+------------------+
| fileExtensions   | *list* | **['mtz']**      | A list of        |
|                  |        |                  | strings          |
|                  |        |                  | containing       |
|                  |        |                  | allowed file     |
|                  |        |                  | extensions (no   |
|                  |        |                  | dot)             |
+------------------+--------+------------------+------------------+
| mimeTypeName     | *str*  | *                |                  |
|                  |        | *application/CCP |                  |
|                  |        | 4-mtz-unmerged** |                  |
+------------------+--------+------------------+------------------+
| mim              | *str*  | **MTZ unmerged   |                  |
| eTypeDescription |        | experimental     |                  |
|                  |        | data**           |                  |
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

CWavelength: Wavelength in Angstrom
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CFloat <./CCP4Data.html#CFloat>`__ ->
`CCP4XtalData.CWavelength <./CCP4XtalData.html#CWavelength>`__

.. _CCP4XtalData-contents-of-class-71:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4XtalData-qualifiers-for-class-70:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+---------+------------------+------------------+
| allowUndefined  | *bool*  | **True**         |                  |
+-----------------+---------+------------------+------------------+
| default         | *float* | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| toolTip         | *str*   | **Data           |                  |
|                 |         | collection       |                  |
|                 |         | wavelength in    |                  |
|                 |         | Angstrom**       |                  |
+-----------------+---------+------------------+------------------+
| guiLabel        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| guiDefinition   | *dict*  | **{}**           |                  |
+-----------------+---------+------------------+------------------+
| helpFile        | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| saveToDb        | *bool*  | **False**        | Save this data   |
|                 |         |                  | in the database  |
+-----------------+---------+------------------+------------------+
| charWidth       | *int*   | **10**           |                  |
+-----------------+---------+------------------+------------------+
| min             | *float* | **0.0**          | The inclusive    |
|                 |         |                  | minimum value    |
+-----------------+---------+------------------+------------------+
| max             | *float* | **None**         | The inclusive    |
|                 |         |                  | maximum value    |
+-----------------+---------+------------------+------------------+
| onlyEnumerators | *bool*  | **False**        | If this is true  |
|                 |         |                  | then the         |
|                 |         |                  | enumerators are  |
|                 |         |                  | obligatory -     |
|                 |         |                  | otherwise they   |
|                 |         |                  | are treated as   |
|                 |         |                  | recommended      |
|                 |         |                  | values           |
+-----------------+---------+------------------+------------------+
| enumerators     | *list*  | **[]**           | A Python list of |
|                 |         |                  | allowed or       |
|                 |         |                  | recommended      |
|                 |         |                  | values - see     |
|                 |         |                  | onlyEnumerators  |
+-----------------+---------+------------------+------------------+
| menuText        | *list*  | **[]**           | A Python list of |
|                 |         |                  | strings,         |
|                 |         |                  | matching items   |
|                 |         |                  | in enumerators   |
|                 |         |                  | list, to appear  |
|                 |         |                  | on GUI menu      |
+-----------------+---------+------------------+------------------+
