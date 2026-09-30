CWorkflowContainerList
~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CDict <./CCP4Data.html#CDict>`__ ->
`CCP4WorkflowManager.CWorkflowContainerList <./CCP4WorkflowManager.html#CWorkflowContainerList>`__

Contents of class:
^^^^^^^^^^^^^^^^^^

CDict of `CContainer <./CCP4Container.html#CContainer>`__

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **{}**             
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================

CWorkflowDataFlow
~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4WorkflowManager.CWorkflowDataFlow <./CCP4WorkflowManager.html#CWorkflowDataFlow>`__

.. _contents-of-class-1:

Contents of class:
^^^^^^^^^^^^^^^^^^

========== =====================================
fromKey    `CString <./CCP4Data.html#CString>`__
toKey      `CString <./CCP4Data.html#CString>`__
fromJob    `CString <./CCP4Data.html#CString>`__
annotation `CString <./CCP4Data.html#CString>`__
========== =====================================

.. _qualifiers-for-class-1:

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

CWorkflowDataFlowList
~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CList <./CCP4Data.html#CList>`__ ->
`CCP4WorkflowManager.CWorkflowDataFlowList <./CCP4WorkflowManager.html#CWorkflowDataFlowList>`__

.. _contents-of-class-2:

Contents of class:
^^^^^^^^^^^^^^^^^^

List of
`CWorkflowDataFlow <./CCP4WorkflowManager.html#CWorkflowDataFlow>`__

.. _qualifiers-for-class-2:

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

CWorkflowDefinition
~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Container.CContainer <./CCP4Container.html#CContainer>`__ ->
`CCP4WorkflowManager.CWorkflowDefinition <./CCP4WorkflowManager.html#CWorkflowDefinition>`__

.. _contents-of-class-3:

Contents of class:
^^^^^^^^^^^^^^^^^^

This is a simple class without contents

.. _qualifiers-for-class-3:

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

CWorkflowFileOut
~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4WorkflowManager.CWorkflowFileOut <./CCP4WorkflowManager.html#CWorkflowFileOut>`__

.. _contents-of-class-4:

Contents of class:
^^^^^^^^^^^^^^^^^^

========= =====================================
className `CString <./CCP4Data.html#CString>`__
key       `CString <./CCP4Data.html#CString>`__
========= =====================================

.. _qualifiers-for-class-4:

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

CWorkflowJobDefinition
~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4WorkflowManager.CWorkflowJobDefinition <./CCP4WorkflowManager.html#CWorkflowJobDefinition>`__

.. _contents-of-class-5:

Contents of class:
^^^^^^^^^^^^^^^^^^

+----------------+----------------------------------------------------+
| input          | `CWorkflowDataFlowList <./                         |
|                | CCP4WorkflowManager.html#CWorkflowDataFlowList>`__ |
+----------------+----------------------------------------------------+
| taskName       | `CString <./CCP4Data.html#CString>`__              |
+----------------+----------------------------------------------------+
| allOutputFiles | `CList <./CCP4Data.html#CList>`__                  |
+----------------+----------------------------------------------------+
| output         | `CWorkflowDataFlowList <./                         |
|                | CCP4WorkflowManager.html#CWorkflowDataFlowList>`__ |
+----------------+----------------------------------------------------+

.. _qualifiers-for-class-5:

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

CWorkflowJobDefinitionDict
~~~~~~~~~~~~~~~~~~~~~~~~~~

`CCP4Data.CData <./CCP4Data.html#CData>`__ ->
`CCP4Data.CCollection <./CCP4Data.html#CCollection>`__ ->
`CCP4Data.CDict <./CCP4Data.html#CDict>`__ ->
`CCP4WorkflowManager.CWorkflowJobDefinitionDict <./CCP4WorkflowManager.html#CWorkflowJobDefinitionDict>`__

.. _contents-of-class-6:

Contents of class:
^^^^^^^^^^^^^^^^^^

CDict of
`CWorkflowJobDefinition <./CCP4WorkflowManager.html#CWorkflowJobDefinition>`__

.. _qualifiers-for-class-6:

Qualifiers for class:
^^^^^^^^^^^^^^^^^^^^^

============== ====== ================== ==============================
allowUndefined *bool* **True**           
default        *dict* **{}**             
toolTip        *str*  **NotImplemented** 
guiLabel       *str*  **NotImplemented** 
guiDefinition  *dict* **{}**             
helpFile       *str*  **NotImplemented** 
saveToDb       *bool* **False**          Save this data in the database
============== ====== ================== ==============================
