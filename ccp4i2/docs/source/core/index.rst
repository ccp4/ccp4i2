Core Modules
============

The core module hierarchy provides the foundation for CCP4i2.

.. toctree::
   :maxdepth: 2

   architecture
   containers

Architecture Overview
---------------------

::

    core/
    ├── base_object/          # Base class infrastructure
    │   ├── base_classes.py   # CData, CContainer
    │   ├── fundamental_types.py  # CInt, CFloat, CString, CBoolean
    │   ├── error_reporting.py    # CErrorReport
    │   └── signal_system.py      # Signal/slot mechanism
    ├── CCP4PluginScript.py   # Pipeline/wrapper base class
    ├── CCP4Container.py      # Container classes
    ├── CCP4Data.py           # Data type definitions
    └── CCP4File.py           # File data types
