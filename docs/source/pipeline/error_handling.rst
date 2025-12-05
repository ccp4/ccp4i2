Error Handling
==============

CCP4i2 provides comprehensive error handling through ``CErrorReport`` and the ``appendErrorReport()`` method.

Key Concepts
------------

ERROR_CODES Dictionary
~~~~~~~~~~~~~~~~~~~~~~

Define error codes in your pipeline class:

.. code-block:: python

    class MyPipeline(CPluginScript):
        ERROR_CODES = {
            101: {'description': 'Failed to read input file'},
            102: {'description': 'Processing error'},
        }

try/except Pattern
~~~~~~~~~~~~~~~~~~

Wrap pipeline phases in try/except:

.. code-block:: python

    def startProcess(self, processId):
        try:
            self.phase1()
        except Exception as e:
            self.appendErrorReport(101, f'Phase 1 failed: {e}')
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

See ``mddocs/pipeline/ERROR_HANDLING_PATTERNS.md`` for comprehensive documentation.
