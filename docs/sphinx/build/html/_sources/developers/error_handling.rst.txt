################################
CCP4I2 Developers Error Handling
################################

Error handling
~~~~~~~~~~~~~~

The *CCP4ErrorHandling* module has two classes *CErrorReport* which
holds a list of errors and *CException* which just subclasses
*CErrorReport* and the Python *Exception*. I use *CCP4ErrorReport* where
a function/method returns a list of errors (e.g. *CData.validity()*)
and\ *CException* when raising an exception but I am not sure sure that
this distinction is necessary. The error report is a list because many
of the context likely to raise errors are looping over many items (e.g.
reading a file, drawing widgets) and we usually do not want errors to
stop the overall process but we do want a report of all errors.

Each library class that might raise an exception defines the error codes
that it uses - for example

::

     class CRange(CData):
         ERROR_CODES = { 101 : { 'description' : 'End of range less than start' },
                         102 : { 'description' : 'End of range greater than start' } }

The description should ideally be helpful to a thoughtful user. The code
which catches the exception can also add more context information to the
error message. There is a slight complication in the numbering of
ERROR_CODES in that errors may also come from a base class so for
derived classes I start the code numbering at 101, 201 etc.. There are a
*lot* of these ERROR_CODES definitions - this is an attempt to give
helpful error messages.

The errors have a severity rating defined in the module:

+--------------------------+---+-------------------------------+
| SEVERITY_OK              | 0 | No error                      |
+--------------------------+---+-------------------------------+
| SEVERITY_UNDEFINED       | 1 | Data undefined for an object  |
|                          |   | for which allowUndefined is   |
|                          |   | True                          |
+--------------------------+---+-------------------------------+
| SEVERITY_WARNING         | 2 | A warning                     |
+--------------------------+---+-------------------------------+
| SEVERITY_UNDEFINED_ERROR | 3 | Data undefined for an object  |
|                          |   | for which allowUndefined is   |
|                          |   | False                         |
+--------------------------+---+-------------------------------+
| SEVERITY_ERROR           | 4 | An error                      |
+--------------------------+---+-------------------------------+

By default the ERROR_CODES are assumed to be SEVERITY_ERROR but can be
specified otherwise:

::

      ERROR_CODES = {
           ....
                          7 : { 'severity' : SEVERITY_WARNING,
                                'description' : 'Unrecognised qualifier in data input' },
           ....
                     }

The *CCP4ErrorHandling.CErrorReport* methods are:

+---------------------+-------------+--------------------------------+
| CErrorReport method | Argument    | Description                    |
+=====================+=============+================================+
| \__init_\_          |             |                                |
+---------------------+-------------+--------------------------------+
| append              |             | Append an error to the report  |
+---------------------+-------------+--------------------------------+
|                     | cls         | The class reporting the error  |
+---------------------+-------------+--------------------------------+
|                     | code        | An error code number unique    |
|                     |             | for the class                  |
+---------------------+-------------+--------------------------------+
|                     | details     | A Python string containing     |
|                     |             | further information            |
+---------------------+-------------+--------------------------------+
|                     | name        | A Python string with an        |
|                     |             | identifier for the object.     |
|                     |             | *CData.objectPath()* gives     |
|                     |             | suitable value for this.       |
+---------------------+-------------+--------------------------------+
| extend              |             | Append errors from another     |
|                     |             | *CErrorReport* object to this  |
|                     |             | object                         |
+---------------------+-------------+--------------------------------+
|                     | other       | A *CErrorReport* or            |
|                     |             | *CException*                   |
+---------------------+-------------+--------------------------------+
| count               |             | Return the number of errors in |
|                     |             | list with the optional         |
|                     |             | specified class or severity    |
+---------------------+-------------+--------------------------------+
|                     | cls         | Optional. A class              |
+---------------------+-------------+--------------------------------+
|                     | severity    | Optional. A minimum severity   |
+---------------------+-------------+--------------------------------+
| maxSeverity         |             | Return the maximum severity of |
|                     |             | the errors in the report       |
+---------------------+-------------+--------------------------------+
| report              |             | Return a text string summary   |
|                     |             | or the errors.                 |
+---------------------+-------------+--------------------------------+
| \__str_\_           |             | As *report()*                  |
+---------------------+-------------+--------------------------------+
| warningMessage      |             | In graphical context presents  |
|                     |             | a dialog box with warning      |
|                     |             | message                        |
+---------------------+-------------+--------------------------------+
|                     | windowTitle | Title for dialog box window    |
+---------------------+-------------+--------------------------------+
|                     | message     | Additional message - should    |
|                     |             | explain context of error       |
+---------------------+-------------+--------------------------------+

*CCP4ErrorHandling.CException* has the same methods but can,optionally,
be initiated with the information for one error:

+-------------------+----------+-------------------------------------+
| CException method | Argument | Description                         |
+===================+==========+=====================================+
| \__init_\_        |          |                                     |
+-------------------+----------+-------------------------------------+
|                   | cls      | The class reporting the error       |
+-------------------+----------+-------------------------------------+
|                   | code     | An error code number unique for the |
|                   |          | class                               |
+-------------------+----------+-------------------------------------+
|                   | details  | A Python string containing further  |
|                   |          | information                         |
+-------------------+----------+-------------------------------------+
|                   | name     | A Python string with an identifier  |
|                   |          | for the object.                     |
|                   |          | *CData.objectPath()* gives suitable |
|                   |          | value for this.                     |
+-------------------+----------+-------------------------------------+

The individual items in the error list can be accessed by the usual list
mechanism (square backets) and each item has the attributes: class,
code, details, name.

Error handling in plugin scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CPluginScript maintains a list of errors and warnings that are reported
in the *diagnostic.xml* file that is presented in the job report if the
job fails. There is a convenience function to append errors to the list:

+----------------------+----------+----------------------------+
| CPluginScript method | Argument | Description                |
+======================+==========+============================+
| appendErrorReport    |          | Append an error to the     |
|                      |          | report                     |
+----------------------+----------+----------------------------+
|                      | code     | An error code number       |
|                      |          | unique for the class       |
+----------------------+----------+----------------------------+
|                      | details  | A Python string containing |
|                      |          | further information        |
+----------------------+----------+----------------------------+
|                      | name     | A Python string with an    |
|                      |          | identifier for the object. |
|                      |          | *CData.objectPath()* gives |
|                      |          | suitable value for this.   |
+----------------------+----------+----------------------------+

For a subclass can define its own error codes - for example for
buccaneer-refmac:

::

   class buccaneer_build_refine(CPluginScript):

       ...
       ERROR_CODES = {  200 : { 'description' : 'Buccaneer task failed' },
                        201 : { 'description' : 'Refmac task failed' },
                        202 : { 'description' : 'Buccaneer output coordinate file not found' },
                        203 : { 'description' : 'Refmac output coordinate file not found' }
                        }

       def process(self):
         ...
           # buccaneer run
           btask = self.makePluginObject('buccaneer')
           ...
           # run the buccaneer sub-task
           jobStatus = btask.process()
           # check sub-task has worked and produced output file
           if jobStatus == CPluginScript.FAILED:
             self.appendErrorReport(200)
             self.reportStatus(jobStatus)
             return
           if not btask.container.outputData.XYZOUT.exists():
             self.appendErrorReport(202,'Expected file: '+str(btask.container.outputData.XYZOUT))
             self.reportStatus(CPluginScript.FAILED)
             return
         

Note that the number codes start at 200 to avoid clashes with the base
class error codes.

The three *CPluginScript* methods that a subclass is likely to
reimplement (*processInputFiles()*, *processOutputFiles()* and
*makeCommandAndScript()*) should return *CPluginScript.FAILED* if there
is an error that will prevent the task continuing (and otherwise return
*CPluginScript.SUCCEEDED*). Before returning a failed status the method
should add an error report explaining the problem.

--------------

Last modified: Thu May 8 15:01:31 BST 2014
