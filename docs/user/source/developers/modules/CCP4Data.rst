CBoolean: A Boolean
~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CBoolean <./CCP4Data.html#CBoolean>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *bool* **NotImplemented** 
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
charWidth      *int*  **10**             
============== ====== ================== ==============================

CFloat: A float
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CFloat <./CCP4Data.html#CFloat>`__

.. _CCP4Data-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-1:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+---------+------------------+------------------+
| allowUndefined  | *bool*  | **True**         |                  |
+-----------------+---------+------------------+------------------+
| default         | *float* | **               |                  |
|                 |         | NotImplemented** |                  |
+-----------------+---------+------------------+------------------+
| toolTip         | *str*   | **               |                  |
|                 |         | NotImplemented** |                  |
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
| min             | *float* | **None**         | The inclusive    |
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

CFloatRange: Two floats defining start and end of range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CRange <./CCP4Data.html#CRange>`__ ->
`CCP4Data.CFloatRange <./CCP4Data.html#CFloatRange>`__

.. _CCP4Data-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

===== ===================================
start `CFloat <./CCP4Data.html#CFloat>`__
end   `CFloat <./CCP4Data.html#CFloat>`__
===== ===================================

.. _CCP4Data-qualifiers-for-class-2:

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
| compare        | *int*  | **               | If value is 1/-1 |
|                |        | NotImplemented** | the end value    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than start.      |
+----------------+--------+------------------+------------------+

CFollowFromJob
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CUUID <./CCP4Data.html#CUUID>`__ ->
`CCP4Data.CFollowFromJob <./CCP4Data.html#CFollowFromJob>`__

.. _CCP4Data-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-3:

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

CI2DataType
~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CI2DataType <./CCP4Data.html#CI2DataType>`__

.. _CCP4Data-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-4:

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
| onlyEnumerators  | *bool* | **False**        | If this is true  |
|                  |        |                  | then the         |
|                  |        |                  | enumerators are  |
|                  |        |                  | obligatory -     |
|                  |        |                  | otherwise they   |
|                  |        |                  | are treated as   |
|                  |        |                  | recommended      |
|                  |        |                  | values           |
+------------------+--------+------------------+------------------+
| enumerators      | *list* | **               | A list of        |
|                  |        | ['CPdbDataFile', | allowed or       |
|                  |        | 'CSeqDataFile',  | recommended      |
|                  |        | 'CObsDataFile',  | values for       |
|                  |        | 'CPhsDataFile',  | string           |
|                  |        | 'CMap            |                  |
|                  |        | CoeffsDataFile', |                  |
|                  |        | '                |                  |
|                  |        | CFreeRDataFile', |                  |
|                  |        | 'CMtzDataFile',  |                  |
|                  |        | 'CDictDataFile', |                  |
|                  |        | 'CDataFile',     |                  |
|                  |        | 'CInt',          |                  |
|                  |        | 'CFloat',        |                  |
|                  |        | 'CString']**     |                  |
+------------------+--------+------------------+------------------+
| menuText         | *list* | **['Atomic model | A list of        |
|                  |        | - CPdbDataFile', | strings          |
|                  |        | 'Sequence -      | equivalent to    |
|                  |        | CSeqDataFile',   | the enumerators  |
|                  |        | 'Reflections -   | that will appear |
|                  |        | CObsDataFile',   | in the GUI       |
|                  |        | 'Phases -        |                  |
|                  |        | CPhsDataFile',   |                  |
|                  |        | 'Map             |                  |
|                  |        | coefficients -   |                  |
|                  |        | CMap             |                  |
|                  |        | CoeffsDataFile', |                  |
|                  |        | 'Free R set -    |                  |
|                  |        | CFreeRDataFile', |                  |
|                  |        | 'Experimental    |                  |
|                  |        | data -           |                  |
|                  |        | CMtzDataFile',   |                  |
|                  |        | 'Geometry        |                  |
|                  |        | dictionary -     |                  |
|                  |        | CDictDataFile',  |                  |
|                  |        | 'CDataFile',     |                  |
|                  |        | 'CInt',          |                  |
|                  |        | 'CFloat',        |                  |
|                  |        | 'CString']**     |                  |
+------------------+--------+------------------+------------------+
| allowedCharsCode | *int*  | **0**            | Flag if the text |
|                  |        |                  | is limited to    |
|                  |        |                  | set of allowed   |
|                  |        |                  | characters       |
+------------------+--------+------------------+------------------+

CInt: An integer
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CInt <./CCP4Data.html#CInt>`__

.. _CCP4Data-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-5:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+--------+------------------+------------------+
| allowUndefined  | *bool* | **True**         |                  |
+-----------------+--------+------------------+------------------+
| default         | *int*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| toolTip         | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| guiLabel        | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| guiDefinition   | *dict* | **{}**           |                  |
+-----------------+--------+------------------+------------------+
| helpFile        | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| saveToDb        | *bool* | **False**        | Save this data   |
|                 |        |                  | in the database  |
+-----------------+--------+------------------+------------------+
| charWidth       | *int*  | **10**           |                  |
+-----------------+--------+------------------+------------------+
| min             | *int*  | **None**         | The inclusive    |
|                 |        |                  | maximum allowed  |
|                 |        |                  | value            |
+-----------------+--------+------------------+------------------+
| max             | *int*  | **None**         | The inclusive    |
|                 |        |                  | minimum allowed  |
|                 |        |                  | value            |
+-----------------+--------+------------------+------------------+
| onlyEnumerators | *bool* | **False**        | If this is true  |
|                 |        |                  | then the         |
|                 |        |                  | enumerators are  |
|                 |        |                  | obligatory -     |
|                 |        |                  | otherwise they   |
|                 |        |                  | are treated as   |
|                 |        |                  | recommended      |
|                 |        |                  | values           |
+-----------------+--------+------------------+------------------+
| enumerators     | *list* | **[]**           | A Python list of |
|                 |        |                  | allowed or       |
|                 |        |                  | recommended      |
|                 |        |                  | values - see     |
|                 |        |                  | onlyEnumerators  |
+-----------------+--------+------------------+------------------+
| menuText        | *list* | **[]**           | A Python list of |
|                 |        |                  | strings,         |
|                 |        |                  | matching items   |
|                 |        |                  | in enumerators   |
|                 |        |                  | list, to appear  |
|                 |        |                  | on GUI menu      |
+-----------------+--------+------------------+------------------+

CIntRange: Two integers defining start and end of range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CRange <./CCP4Data.html#CRange>`__ ->
`CCP4Data.CIntRange <./CCP4Data.html#CIntRange>`__

.. _CCP4Data-contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

===== ===============================
start `CInt <./CCP4Data.html#CInt>`__
end   `CInt <./CCP4Data.html#CInt>`__
===== ===============================

.. _CCP4Data-qualifiers-for-class-6:

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
| compare        | *int*  | **               | If value is 1/-1 |
|                |        | NotImplemented** | the end value    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than start.      |
+----------------+--------+------------------+------------------+

CJobStatus
~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CInt <./CCP4Data.html#CInt>`__ ->
`CCP4Data.CJobStatus <./CCP4Data.html#CJobStatus>`__

.. _CCP4Data-contents-of-class-7:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-7:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+--------+------------------+------------------+
| allowUndefined  | *bool* | **True**         |                  |
+-----------------+--------+------------------+------------------+
| default         | *int*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| toolTip         | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| guiLabel        | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| guiDefinition   | *dict* | **{}**           |                  |
+-----------------+--------+------------------+------------------+
| helpFile        | *str*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| saveToDb        | *bool* | **False**        | Save this data   |
|                 |        |                  | in the database  |
+-----------------+--------+------------------+------------------+
| charWidth       | *int*  | **10**           |                  |
+-----------------+--------+------------------+------------------+
| min             | *int*  | **None**         | The inclusive    |
|                 |        |                  | maximum allowed  |
|                 |        |                  | value            |
+-----------------+--------+------------------+------------------+
| max             | *int*  | **None**         | The inclusive    |
|                 |        |                  | minimum allowed  |
|                 |        |                  | value            |
+-----------------+--------+------------------+------------------+
| onlyEnumerators | *bool* | **False**        | If this is true  |
|                 |        |                  | then the         |
|                 |        |                  | enumerators are  |
|                 |        |                  | obligatory -     |
|                 |        |                  | otherwise they   |
|                 |        |                  | are treated as   |
|                 |        |                  | recommended      |
|                 |        |                  | values           |
+-----------------+--------+------------------+------------------+
| enumerators     | *list* | **[]**           | A Python list of |
|                 |        |                  | allowed or       |
|                 |        |                  | recommended      |
|                 |        |                  | values - see     |
|                 |        |                  | onlyEnumerators  |
+-----------------+--------+------------------+------------------+
| menuText        | *list* | **[]**           | A Python list of |
|                 |        |                  | strings,         |
|                 |        |                  | matching items   |
|                 |        |                  | in enumerators   |
|                 |        |                  | list, to appear  |
|                 |        |                  | on GUI menu      |
+-----------------+--------+------------------+------------------+

CJobTitle
~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CJobTitle <./CCP4Data.html#CJobTitle>`__

.. _CCP4Data-contents-of-class-8:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-8:

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

COneWord: A single word string - no white space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.COneWord <./CCP4Data.html#COneWord>`__

.. _CCP4Data-contents-of-class-9:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-9:

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

CPatchSelection
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CPatchSelection <./CCP4Data.html#CPatchSelection>`__

.. _CCP4Data-contents-of-class-10:

Contents of class:
^^^^^^^^^^^^^^^^^^

======== =====================================
taskName `CString <./CCP4Data.html#CString>`__
patch    `CString <./CCP4Data.html#CString>`__
======== =====================================

.. _CCP4Data-qualifiers-for-class-10:

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

CRange: Base class for CIntRange and CFloatRange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CRange <./CCP4Data.html#CRange>`__

.. _CCP4Data-contents-of-class-11:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-11:

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
| compare        | *int*  | **               | If value is 1/-1 |
|                |        | NotImplemented** | the end value    |
|                |        |                  | must be          |
|                |        |                  | greater/less     |
|                |        |                  | than start.      |
+----------------+--------+------------------+------------------+

CRangeSelection
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CRangeSelection <./CCP4Data.html#CRangeSelection>`__

.. _CCP4Data-contents-of-class-12:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-12:

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

CString: A string
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__

.. _CCP4Data-contents-of-class-13:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-13:

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

CTable
~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CTable <./CCP4Data.html#CTable>`__

.. _CCP4Data-contents-of-class-14:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-14:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
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
| saveToDb         | *bool* | **False**        | Save this data   |
|                  |        |                  | in the database  |
+------------------+--------+------------------+------------------+
| allo             | *bool* | **False**        |                  |
| wUndefinedHeader |        |                  |                  |
+------------------+--------+------------------+------------------+
| uniqueHeader     | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| rowLabel         | *str*  | **row**          |                  |
+------------------+--------+------------------+------------------+
| columnLabel      | *str*  | **column**       |                  |
+------------------+--------+------------------+------------------+

CUUID
~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Data.CUUID <./CCP4Data.html#CUUID>`__

.. _CCP4Data-contents-of-class-15:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Data-qualifiers-for-class-15:

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
