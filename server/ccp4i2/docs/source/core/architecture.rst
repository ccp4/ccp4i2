Architecture
============

CCP4i2 uses a hierarchical data model for parameters and files.

Class Hierarchy
---------------

::

    HierarchicalObject (base_classes.py)
        │
        └── CData (base_classes.py)
            │
            ├── CContainer (base_classes.py)
            │   ├── inputData
            │   ├── outputData
            │   └── controlParameters
            │
            ├── Fundamental Types (fundamental_types.py)
            │   ├── CInt
            │   ├── CFloat
            │   ├── CString
            │   └── CBoolean
            │
            └── CDataFile (CCP4File.py)
                ├── CPdbDataFile
                ├── CMtzDataFile
                └── ...

CPluginScript
-------------

``CPluginScript`` extends ``CData`` with plugin infrastructure:

- Process lifecycle (startProcess, processOutputFiles)
- Sub-plugin creation (makePluginObject)
- Signal/slot connections (connectSignal)
- Error reporting (appendErrorReport)
- Database integration
