CAngle: An angle
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CBaseData <./CCP4Data.html#CBaseData>`__ ->
`CCP4Data.CFloat <./CCP4Data.html#CFloat>`__ ->
`CCP4MathsData.CAngle <./CCP4MathsData.html#CAngle>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

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

CEulerRotation
~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4MathsData.CEulerRotation <./CCP4MathsData.html#CEulerRotation>`__

.. _CCP4MathsData-contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

===== ========================================
alpha `CAngle <./CCP4MathsData.html#CAngle>`__
beta  `CAngle <./CCP4MathsData.html#CAngle>`__
gamma `CAngle <./CCP4MathsData.html#CAngle>`__
===== ========================================

.. _CCP4MathsData-qualifiers-for-class-1:

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

CMatrix33
~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4MathsData.CMatrix33 <./CCP4MathsData.html#CMatrix33>`__

.. _CCP4MathsData-contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _CCP4MathsData-qualifiers-for-class-2:

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

CTransformation
~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4MathsData.CTransformation <./CCP4MathsData.html#CTransformation>`__

.. _CCP4MathsData-contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

=========== ========================================================
translation `CXyz <./CCP4MathsData.html#CXyz>`__
rotation    `CEulerRotation <./CCP4MathsData.html#CEulerRotation>`__
=========== ========================================================

.. _CCP4MathsData-qualifiers-for-class-3:

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

CXyz
~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4MathsData.CXyz <./CCP4MathsData.html#CXyz>`__

.. _CCP4MathsData-contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

= ===================================
y `CFloat <./CCP4Data.html#CFloat>`__
x `CFloat <./CCP4Data.html#CFloat>`__
z `CFloat <./CCP4Data.html#CFloat>`__
= ===================================

.. _CCP4MathsData-qualifiers-for-class-4:

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

CXyzBox
~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4MathsData.CXyzBox <./CCP4MathsData.html#CXyzBox>`__

.. _CCP4MathsData-contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

==== ===================================
zMax `CFloat <./CCP4Data.html#CFloat>`__
yMax `CFloat <./CCP4Data.html#CFloat>`__
zMin `CFloat <./CCP4Data.html#CFloat>`__
xMax `CFloat <./CCP4Data.html#CFloat>`__
xMin `CFloat <./CCP4Data.html#CFloat>`__
yMin `CFloat <./CCP4Data.html#CFloat>`__
==== ===================================

.. _CCP4MathsData-qualifiers-for-class-5:

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
