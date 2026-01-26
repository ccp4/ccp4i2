CPluginScript
=============

.. automodule:: core.CCP4PluginScript
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

CPluginScript is the base class for all CCP4i2 wrappers and pipelines.

Key Methods
-----------

Lifecycle Methods
~~~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   CPluginScript.process
   CPluginScript.startProcess
   CPluginScript.processOutputFiles
   CPluginScript.postProcess
   CPluginScript.reportStatus

Validation Methods
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   CPluginScript.validity
   CPluginScript.checkInputData
   CPluginScript.checkOutputData

Plugin Creation
~~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   CPluginScript.makePluginObject
   CPluginScript.connectSignal

Error Handling
~~~~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   CPluginScript.appendErrorReport

Class Attributes
----------------

Status Codes
~~~~~~~~~~~~

- ``SUCCEEDED = 0`` - Job completed successfully
- ``FAILED = 1`` - Job failed
- ``RUNNING = 2`` - Job is running
- ``UNSATISFACTORY = 3`` - Job completed with warnings

Plugin Metadata
~~~~~~~~~~~~~~~

- ``TASKNAME`` - Unique task identifier
- ``TASKTITLE`` - Display title for GUI
- ``TASKCOMMAND`` - Executable name
- ``ASYNCHRONOUS`` - Whether to run asynchronously
