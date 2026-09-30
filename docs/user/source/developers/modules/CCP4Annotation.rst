CAnnotation: Annotation text with user id and time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CAnnotation <./CCP4Annotation.html#CAnnotation>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

====== ===========================================
text   `CString <./CCP4Data.html#CString>`__
author `CUserId <./CCP4Annotation.html#CUserId>`__
time   `CTime <./CCP4Annotation.html#CTime>`__
====== ===========================================

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **Enter your     |                  |
|                |        | comments**       |                  |
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

CAnnotationList: A list of annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4Annotation.CAnnotationList <./CCP4Annotation.html#CAnnotationList>`__

.. _CCP4Annotation-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CAnnotation <./CCP4Annotation.html#CAnnotation>`__

.. _CCP4Annotation-qualifiers-for-class-1:

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

CAuthor: Placeholder for bibliographic author
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Annotation.CAuthor <./CCP4Annotation.html#CAuthor>`__

.. _CCP4Annotation-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Annotation-qualifiers-for-class-2:

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

CBibReference: Bibliographic reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CBibReference <./CCP4Annotation.html#CBibReference>`__

.. _CCP4Annotation-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =======================================
authorList `CList <./CCP4Data.html#CList>`__
title      `CString <./CCP4Data.html#CString>`__
url        `CString <./CCP4Data.html#CString>`__
selected   `CBoolean <./CCP4Data.html#CBoolean>`__
source     `CString <./CCP4Data.html#CString>`__
pmid       `CInt <./CCP4Data.html#CInt>`__
========== =======================================

.. _CCP4Annotation-qualifiers-for-class-3:

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

CBibReferenceGroup: Set of bibliographic references for a task
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CBibReferenceGroup <./CCP4Annotation.html#CBibReferenceGroup>`__

.. _CCP4Annotation-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =====================================
taskName   `CString <./CCP4Data.html#CString>`__
version    `CString <./CCP4Data.html#CString>`__
references `CList <./CCP4Data.html#CList>`__
title      `CString <./CCP4Data.html#CString>`__
========== =====================================

.. _CCP4Annotation-qualifiers-for-class-4:

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

CFont: Simplified Qt font options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CFont <./CCP4Annotation.html#CFont>`__

.. _CCP4Annotation-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

========= =====================================
style     `CInt <./CCP4Data.html#CInt>`__
pointSize `CInt <./CCP4Data.html#CInt>`__
weight    `CInt <./CCP4Data.html#CInt>`__
family    `CString <./CCP4Data.html#CString>`__
========= =====================================

.. _CCP4Annotation-qualifiers-for-class-5:

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

CHostName: Computer name
~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Annotation.CHostName <./CCP4Annotation.html#CHostName>`__

.. _CCP4Annotation-contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Annotation-qualifiers-for-class-6:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Hostname as    |                  |
|                  |        | mycomputer       |                  |
|                  |        | .myplace.ac.uk** |                  |
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

CHostname
~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Annotation.CHostName <./CCP4Annotation.html#CHostName>`__ ->
`CCP4Annotation.CHostname <./CCP4Annotation.html#CHostname>`__

.. _CCP4Annotation-contents-of-class-7:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Annotation-qualifiers-for-class-7:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **Hostname as    |                  |
|                  |        | mycomputer       |                  |
|                  |        | .myplace.ac.uk** |                  |
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

CMetaDataTag: This class will extend list of enumerators if new value for string is entered
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CMetaDataTag <./CCP4Annotation.html#CMetaDataTag>`__

.. _CCP4Annotation-contents-of-class-8:

Contents of class:
^^^^^^^^^^^^^^^^^^

=== =====================================
tag `CString <./CCP4Data.html#CString>`__
=== =====================================

.. _CCP4Annotation-qualifiers-for-class-8:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+----------------+----------------+----------------+
| allowUndefined | *bool*         | **True**       |                |
+----------------+----------------+----------------+----------------+
| default        | *dict*         | **No           |                |
|                |                | tImplemented** |                |
+----------------+----------------+----------------+----------------+
| toolTip        | *str*          | **No           |                |
|                |                | tImplemented** |                |
+----------------+----------------+----------------+----------------+
| guiLabel       | *str*          | **No           |                |
|                |                | tImplemented** |                |
+----------------+----------------+----------------+----------------+
| guiDefinition  | *dict*         | **{}**         |                |
+----------------+----------------+----------------+----------------+
| helpFile       | *str*          | **No           |                |
|                |                | tImplemented** |                |
+----------------+----------------+----------------+----------------+
| saveToDb       | *bool*         | **False**      | Save this data |
|                |                |                | in the         |
|                |                |                | database       |
+----------------+----------------+----------------+----------------+
| enume          | *i             | **None**       |                |
| ratorsFunction | nstancemethod* |                |                |
+----------------+----------------+----------------+----------------+
| addEnum        | *i             | **None**       |                |
| eratorFunction | nstancemethod* |                |                |
+----------------+----------------+----------------+----------------+

CMetaDataTagList
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4Annotation.CMetaDataTagList <./CCP4Annotation.html#CMetaDataTagList>`__

.. _CCP4Annotation-contents-of-class-9:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of `CMetaDataTag <./CCP4Annotation.html#CMetaDataTag>`__

.. _CCP4Annotation-qualifiers-for-class-9:

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

CServerGroup: One or more compute servers used in "remote" running
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CServerGroup <./CCP4Annotation.html#CServerGroup>`__

.. _CCP4Annotation-contents-of-class-10:

Contents of class:
^^^^^^^^^^^^^^^^^^

================ =========================================
customCodeFile   `CDataFile <./CCP4File.html#CDataFile>`__
keyFilename      `CString <./CCP4Data.html#CString>`__
ccp4Dir          `CString <./CCP4Data.html#CString>`__
tempDir          `CString <./CCP4Data.html#CString>`__
validate         `CString <./CCP4Data.html#CString>`__
maxTries         `CInt <./CCP4Data.html#CInt>`__
sge_root         `CString <./CCP4Data.html#CString>`__
name             `CString <./CCP4Data.html#CString>`__
mechanism        `CString <./CCP4Data.html#CString>`__
timeout          `CFloat <./CCP4Data.html#CFloat>`__
queueOptionsFile `CDataFile <./CCP4File.html#CDataFile>`__
serverList       `CList <./CCP4Data.html#CList>`__
userExtensible   `CBoolean <./CCP4Data.html#CBoolean>`__
================ =========================================

.. _CCP4Annotation-qualifiers-for-class-10:

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

CTime: The time. Uses Python time module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CInt <./CCP4Data.html#CInt>`__ ->
`CCP4Annotation.CTime <./CCP4Annotation.html#CTime>`__

.. _CCP4Annotation-contents-of-class-11:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Annotation-qualifiers-for-class-11:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+-----------------+--------+------------------+------------------+
| allowUndefined  | *bool* | **True**         |                  |
+-----------------+--------+------------------+------------------+
| default         | *int*  | **               |                  |
|                 |        | NotImplemented** |                  |
+-----------------+--------+------------------+------------------+
| toolTip         | *str*  | **Time and date  |                  |
|                 |        | as hh:mm         |                  |
|                 |        | dd/mm/yyyy**     |                  |
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
| min             | *int*  | **0**            | The inclusive    |
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
| format          | *str*  | **%H:%M          | Argument to      |
|                 |        | %d/%b/%y**       | Python           |
|                 |        |                  | time.strftime to |
|                 |        |                  | display time in  |
|                 |        |                  | human readable   |
|                 |        |                  | format           |
+-----------------+--------+------------------+------------------+

CUserAddress: User id and platform node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Annotation.CUserAddress <./CCP4Annotation.html#CUserAddress>`__

.. _CCP4Annotation-contents-of-class-12:

Contents of class:
^^^^^^^^^^^^^^^^^^

============ ===========================================
userId       `CUserId <./CCP4Annotation.html#CUserId>`__
platformNode `CString <./CCP4Data.html#CString>`__
============ ===========================================

.. _CCP4Annotation-qualifiers-for-class-12:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+----------------+--------+------------------+------------------+
| allowUndefined | *bool* | **True**         |                  |
+----------------+--------+------------------+------------------+
| default        | *dict* | **               |                  |
|                |        | NotImplemented** |                  |
+----------------+--------+------------------+------------------+
| toolTip        | *str*  | **User id as     |                  |
|                |        | me@myplace.ac.uk |                  |
|                |        | and machine      |                  |
|                |        | name**           |                  |
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

CUserId: A user ID
~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CString <./CCP4Data.html#CString>`__ ->
`CCP4Annotation.CUserId <./CCP4Annotation.html#CUserId>`__

.. _CCP4Annotation-contents-of-class-13:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4Annotation-qualifiers-for-class-13:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

+------------------+--------+------------------+------------------+
| allowUndefined   | *bool* | **True**         |                  |
+------------------+--------+------------------+------------------+
| default          | *str*  | **               |                  |
|                  |        | NotImplemented** |                  |
+------------------+--------+------------------+------------------+
| toolTip          | *str*  | **User id as     |                  |
|                  |        | me               |                  |
|                  |        | @myplace.ac.uk** |                  |
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
