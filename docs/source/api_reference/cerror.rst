Error Handling
==============

.. automodule:: core.base_object.error_reporting
   :members:
   :undoc-members:
   :show-inheritance:

CErrorReport
------------

``CErrorReport`` is the container for tracking validation errors and warnings.

Severity Levels
~~~~~~~~~~~~~~~

- ``SEVERITY_OK = 0`` - No error
- ``SEVERITY_UNDEFINED = 1`` - Value not set
- ``SEVERITY_WARNING = 2`` - Warning (non-blocking)
- ``SEVERITY_UNDEFINED_ERROR = 3`` - Required value missing
- ``SEVERITY_ERROR = 4`` - Fatal error

Key Methods
~~~~~~~~~~~

.. autosummary::
   :nosignatures:

   CErrorReport.append
   CErrorReport.extend
   CErrorReport.maxSeverity
   CErrorReport.report
   CErrorReport.getErrors
   CErrorReport.getEtree

Usage Example
~~~~~~~~~~~~~

.. code-block:: python

    from core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR

    # Create error report
    error = CErrorReport()

    # Add an error
    error.append(
        klass="MyPlugin",
        code=101,
        details="File not found",
        name="XYZIN",
        severity=SEVERITY_ERROR
    )

    # Check severity
    if error.maxSeverity() >= SEVERITY_ERROR:
        print(error.report())

CException
----------

``CException`` combines ``CErrorReport`` with Python's ``Exception`` class,
allowing errors to be both raised and tracked.

.. code-block:: python

    from core.base_object.error_reporting import CException, SEVERITY_ERROR

    raise CException("MyPlugin", 101, "File not found", "XYZIN", SEVERITY_ERROR)
