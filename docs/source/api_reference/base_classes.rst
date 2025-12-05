Base Classes
============

.. automodule:: core.base_object.base_classes
   :members:
   :undoc-members:
   :show-inheritance:

The base class hierarchy provides the foundation for all CCP4i2 data objects.

Class Hierarchy
---------------

::

    HierarchicalObject
        └── CData
            ├── CContainer
            │   ├── inputData
            │   ├── outputData
            │   └── controlParameters
            └── CDataFile
                ├── CPdbDataFile
                ├── CMtzDataFile
                └── ...

HierarchicalObject
~~~~~~~~~~~~~~~~~~

Provides:

- Parent/child relationships
- Signal/slot mechanism
- Object naming and paths
- Event system integration

CData
~~~~~

Adds:

- Value storage and retrieval
- Qualifier system (allowUndefined, mustExist, etc.)
- Validation via ``validity()``
- XML serialization
