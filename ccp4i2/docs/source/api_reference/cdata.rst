CData
=====

.. automodule:: core.CCP4Data
   :members:
   :undoc-members:
   :show-inheritance:

CData is the base class for all data objects in the CCP4i2 container hierarchy.

Key Classes
-----------

CData
~~~~~

The fundamental data object providing:

- Hierarchical parent/child relationships
- Value get/set operations
- Validation through ``validity()`` method
- XML serialization via ``getEtree()``

CList
~~~~~

Container for lists of CData objects:

- ``append(item)`` - Add item to list
- ``makeItem()`` - Create new item of appropriate type
- ``listMinLength``, ``listMaxLength`` - Length constraints

Common Methods
--------------

Value Operations
~~~~~~~~~~~~~~~~

- ``set(value)`` - Set value from dict, CData, or primitive
- ``get()`` - Get value as dict
- ``isSet(field_name)`` - Check if field has value
- ``unSet(field_name)`` - Clear field value
- ``setToDefault()`` - Reset to default value

Validation
~~~~~~~~~~

- ``validity()`` - Return CErrorReport with validation errors
- ``getQualifier(name)`` - Get qualifier value
- ``setQualifier(name, value)`` - Set qualifier value

Serialization
~~~~~~~~~~~~~

- ``getEtree()`` - Convert to XML ElementTree
- ``setEtree(element)`` - Load from XML ElementTree
