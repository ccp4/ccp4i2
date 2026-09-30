CCustomComFile
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4CustomTaskManager.CCustomComFile <./CCP4CustomTaskManager.html#CCustomComFile>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

==== =====================================
text `CString <./CCP4Data.html#CString>`__
name `CString <./CCP4Data.html#CString>`__
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

CCustomComFileList
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4CustomTaskManager.CCustomComFileList <./CCP4CustomTaskManager.html#CCustomComFileList>`__

.. _CCP4CustomTaskManager-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CCustomComFile <./CCP4CustomTaskManager.html#CCustomComFile>`__

.. _CCP4CustomTaskManager-qualifiers-for-class-1:

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

CCustomTaskDefinition
~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Container.CContainer <./CCP4Container.html#CContainer>`__ ->
`CCP4CustomTaskManager.CCustomTaskDefinition <./CCP4CustomTaskManager.html#CCustomTaskDefinition>`__

.. _CCP4CustomTaskManager-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4CustomTaskManager-qualifiers-for-class-2:

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

CCustomTaskFileFunction
~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4CustomTaskManager.CCustomTaskFileFunction <./CCP4CustomTaskManager.html#CCustomTaskFileFunction>`__

.. _CCP4CustomTaskManager-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4CustomTaskManager-qualifiers-for-class-3:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+--------+------------------+------------------+
| allowUndefined  | *bool* | **True**         |                  |
+-----------------+--------+------------------+------------------+
| default         | *str*  | **               |                  |
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
| charWidth       | *int*  | **-1**           |                  |
+-----------------+--------+------------------+------------------+
| minLength       | *int*  | **None**         | Minimum length   |
|                 |        |                  | of string        |
+-----------------+--------+------------------+------------------+
| maxLength       | *int*  | **None**         | Maximum length   |
|                 |        |                  | of string        |
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
| enumerators     | *list* | **['unknown',    | A list of        |
|                 |        | 'input',         | allowed or       |
|                 |        | 'output',        | recommended      |
|                 |        | 'control         | values for       |
|                 |        | parameter',      | string           |
|                 |        | 'log']**         |                  |
+-----------------+--------+------------------+------------------+
| menuText        | *list* | **[]**           | A list of        |
|                 |        |                  | strings          |
|                 |        |                  | equivalent to    |
|                 |        |                  | the enumerators  |
|                 |        |                  | that will appear |
|                 |        |                  | in the GUI       |
+-----------------+--------+------------------+------------------+

CCustomTaskParam
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4CustomTaskManager.CCustomTaskParam <./CCP4CustomTaskManager.html#CCustomTaskParam>`__

.. _CCP4CustomTaskManager-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

+---------------------+-----------------------------------------------+
| obligatory          | `CBoolean <./CCP4Data.html#CBoolean>`__       |
+---------------------+-----------------------------------------------+
| function            | `CCustomTaskFileFunction <./CCP4Custo         |
|                     | mTaskManager.html#CCustomTaskFileFunction>`__ |
+---------------------+-----------------------------------------------+
| name                | `COneWord <./CCP4Data.html#COneWord>`__       |
+---------------------+-----------------------------------------------+
| saveDataToDb        | `CBoolean <./CCP4Data.html#CBoolean>`__       |
+---------------------+-----------------------------------------------+
| requiredContentType | `CList <./CCP4Data.html#CList>`__             |
+---------------------+-----------------------------------------------+
| dataType            | `CI2DataType <./CCP4Data.html#CI2DataType>`__ |
+---------------------+-----------------------------------------------+
| mergeTo             | `CString <./CCP4Data.html#CString>`__         |
+---------------------+-----------------------------------------------+
| outputFilePath      | `CString <./CCP4Data.html#CString>`__         |
+---------------------+-----------------------------------------------+
| splitColumns        | `CString <./CCP4Data.html#CString>`__         |
+---------------------+-----------------------------------------------+
| label               | `CString <./CCP4Data.html#CString>`__         |
+---------------------+-----------------------------------------------+

.. _CCP4CustomTaskManager-qualifiers-for-class-4:

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

CCustomTaskParamList
~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4CustomTaskManager.CCustomTaskParamList <./CCP4CustomTaskManager.html#CCustomTaskParamList>`__

.. _CCP4CustomTaskManager-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of
`CCustomTaskParam <./CCP4CustomTaskManager.html#CCustomTaskParam>`__

.. _CCP4CustomTaskManager-qualifiers-for-class-5:

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
