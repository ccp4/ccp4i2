CDataFile
=========

.. automodule:: ccp4i2.core.CCP4File
   :members:
   :undoc-members:
   :show-inheritance:

CDataFile and related classes handle file-based data in CCP4i2.

Key Classes
-----------

CDataFile
~~~~~~~~~

Base class for all file data types. Its fields, read from the class itself:

.. cdata-fields:: ccp4i2.core.base_object.cdata_file.CDataFile

The ``fullPath`` property composes ``project``, ``relPath`` and ``baseName``
into a filesystem path; ``mustExist`` is a qualifier rather than a field, and
governs whether validation requires the file to be present.

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

Error codes
-----------

.. cdata-errors:: ccp4i2.core.base_object.cdata_file.CDataFile
