CDataFile
=========

.. automodule:: core.CCP4File
   :members:
   :undoc-members:
   :show-inheritance:

CDataFile and related classes handle file-based data in CCP4i2.

Key Classes
-----------

CDataFile
~~~~~~~~~

Base class for all file data types. Provides:

- ``fullPath`` - Full filesystem path
- ``annotation`` - Human-readable description
- ``mustExist`` - Whether file must exist for validation

Specialized File Types
~~~~~~~~~~~~~~~~~~~~~~

- ``CPdbDataFile`` - PDB/mmCIF coordinate files
- ``CMtzDataFile`` - MTZ reflection files
- ``CMapDataFile`` - CCP4 map files
- ``CSeqDataFile`` - Sequence files

Common Methods
--------------

.. code-block:: python

    # Set file path
    file_obj.setFullPath("/path/to/file.pdb")

    # Get file path
    path = file_obj.getFullPath()

    # Check if file exists
    if file_obj.exists():
        ...

    # Get file content info
    content = file_obj.getFileContent()
