Error codes
===========

When a job fails, CCP4i2 reports an error code alongside a message. If you have
only the code — from a log, a report, or a colleague — this page is where to
look it up.

Codes are numbered **per task**, not globally: only about forty distinct
numbers are in use across nearly eighty tasks, so code 301 means one thing in
``molrep_mr`` and something else in ``servalcat_pipe``. Find the row for your
code *and* the task that raised it.

Task error codes
----------------

.. error-code-index::

Data object error codes
-----------------------

These are raised by the data objects themselves rather than by a task, and are
shared by every parameter in the system. They generally mean a value is
missing, is of the wrong type, or could not be read from or written to a file.

.. cdata-errors:: ccp4i2.core.base_object.cdata.CData
