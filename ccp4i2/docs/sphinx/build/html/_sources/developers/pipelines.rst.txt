CCP4I2 Developers: Wrappers and Pipelines
=========================================

Introduction
~~~~~~~~~~~~

To make a program or script available in CCP4i2 it must be 'wrapped' in
a Python script that is a sub-class of the *CPluginScript* class.
*CPluginScript* provides methods to handle parameters, run programs,
maintain an error report and record progress to a database so writing a
simple script requires very minimal lines of code. We usually call a
script that runs just one program a 'wrapper'; and a 'pipeline' is a
script that calls two or more wrappers in order to run multiple
programs. We also refer to either wrappers or pipelines as 'plugins'. We
hope that developers will contribute wrappers to the *ccp4i2/wrappers*
directory to be available to all pipeline developers and the GUI. This
is the place to look for examples to follow. There is a separate
*ccp4i2/pipelines* directory for contributed pipelines.

For each script we need to clearly define the interface (i.e. the
parameters that are passed into the script); this is vital both for
communication with the GUI and for enabling the database to track the
provenance of crystallographic data. The parameters of the script
interface must follow the CCP4i2 data model that has classes for all the
common crystallographic data (well we are working on it!). These
parameters are listed in a *DEF* file; this is an XML file (see
`appendix <#xml_files>`__) and there is one for each *CPluginScript*
script file. The *DEF* can be created by hand-editing but can more
easily be created by the developer using the *defEd* GUI that provides a
list and documentation for all the classes of the CCP4i2 data model.
Note that the *DEF* file is a definition of the parameters; it does not
contain the actual values for parameters.

Most *CPluginScript* scripts are made accessible in the CCP4i2 GUI. The
GUI uses the same data model as the scripts and has widgets to represent
each data class so a GUI for the script can be auto-generated from the
*DEF* file. The auto-generated GUI is useful for testing purposes; a
better structured and more attractive GUI can be written based on the
*DEF* file. When a job is run the GUI writes an *input_params* file that
contains the values of parameters that are passed to the script.

The biggest job in creating a wrapper script is usually generating an
input command file to the program that is correctly formatted and
interprets the parameters passed to the script. The easiest way to do
this is to write a template for the program command file and then use
the *CComTemplate* class to interpret the template and the script
parameters to create a program command file.

In the CCP4i2 GUI the scripts are referred to as 'tasks' and every time
the user runs a task it is called a 'job'. If a *CPluginScript* calls
another *CPluginScript* then this is called a 'sub-job'. The CCP4I2
database records all jobs and any sub-jobs and the input and output data
for each job. This input and output data is mostly in the form of
references to PDB and MTZ files but other data (usually referred to as
'extra' data) is also recorded. The database is notified when a job is
started and when it finishes. When a job starts from the GUI the
database is notified automatically. When a *CPluginScript* script calls
other scripts it should use the *makePluginObject()* method that
instantiates the new *CPluginScript*, passes necessary 'database'
information to the new object and notifies the database that a sub-job
is started. The completion of jobs is notified by the *reportStatus()*
method that is called by the base class *postProcess()* method. The
*reportStatus()* method also saves the script's parameters to a *PARAMS*
file that is read automatically to load the appropriate data to the
database. The *PARAMS* file has the same format and mostly the same
content as the *input_params* file but may contain additional 'output'
information.

*CPluginScript* maintains a list of all errors or warnings generated
within the class or by functions called from the class. This serves as a
complete log - note that it is possible to add 'error' reports that are
not an error (severity is SEVERITY_OK).

Input and Output Data
~~~~~~~~~~~~~~~~~~~~~

Firstly to sort out a naming convention. There are two sorts of xml file
that are relvant here: the 'def' file is a definition of the parameters
and the 'params' file contain actual data for a particular job. I used
to call the 'params' file a 'data' file but that is confusing with the
real data files such as PDB/MTX.

The contents of a def file are split into three or four sub-containers.
It is very important that data is put in the right sub-container. The
*inputData* and *outputData* containers are for the input and output
experimental data and model data. This is usually the PDB and MTZ files
but is also the 'extra' data such as NCS and domains which are an
essential part of the description of the model. The *controlParameters*
subcontainer should contain just about everything else. Parameters
relevant for implementation of the gui go in the *guiAdmin*
sub-container.

The contents of the *inputData* and *outputData* define the *signature*
of a script and is vital for the correct management of the data in the
overall i2 system. Many programs can have alternative input and output
dependent on the exact function. As a general principle we could say
that the different signatures ought to be treated as different tasks in
CCP4i2 and so should have different wrappers. There are probably good
reasons for not totally enforcing that principle and writing multiple
wrappers is inefficient. See how it goes.

The GUI automatically populates the input data based on the previous
jobs but the user can override this and select files or data. The task
interfaces do not allow the user to set output data - the user can
export (i.e. copy) files after the job is finished. The GUI will write a
params file, *whatever.input_params.def.xml* that is read by the wrapper
(this is done automatically as the wrapper is instatiated). The names of
output files are set automatically by *CPluginScript.checkOutputData()*
method that is called as part of the base class *process()* method and
if you reimplement this then you must ensure that *checkOutputData()* is
called before *makeCommandAndScript()*. The *makeCommandAndScript()*
will pass-on the names of ouput files to the program arguments or input
file. The ouput data that does not go into a file, i.e. the 'extra' data
must be handled by the script. What needs doing will depend on the data
and the program but typically files might need to be copied to the
'correct' path specified by the *outputDta* params and other data might
be read from the program log file and loaded into the appropriate
*outputData* parameters. The updated version of the parameters are
written to a file *whatever.params.xml* by the
*CPluginScript.saveParams()* method that is called as part of
*CPluginScript.postProcess()*. The main GUI process extracts the ouput
data from this file and saves to the database.

It is obviously important to identify the 'extra' data that is passed
between programs and to use a standard data model.

A Simple Wrapper Example - *pdbset*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This wrapper changes the cell parameter in a PDB file using the *pdbset*
program. (A better way to do this now would use Python interface to
*mmdb* library). The input to the script is a PDB file name and cell
values. The output is another PDB file name. The DEF file to define the
script parameters is
`$CCP4I2/wrappers/pdbset/script/pdbset.def.xml <../../wrappers/pdbset/script/pdbset.def.xml>`__
which contains (excluding the meta-data header):

::

       <container id="inputData">
         <content id="XYZIN">
           <className>CPdbDataFile</className>
           <qualifiers>
              <mustExist>True</mustExist>
           </qualifiers>
         </content>
         <content id="CELL">
           <className>CCell</className>
         </content>
       </container>
       <container id="outputData">
         <content id="XYZOUT">
           <className>CPdbDataFile</className>
         </content>
       </container>
       <container id="controlParameters">
       </container>

The two important types of element in this XML are <container> and
<content>. The *container* is used to group together the parameters into
three groups: *inputData*, *controlParameters* and *outputData* (this
example does not have any *controlParameters*). The *content* elements
specify the name (as the id attribute), data class (as *className*
element), and *qualifiers* for each parameter. So, for example, the
*XYZIN* parameter is class *CPdbDataFile* and has a qualifier
*mustExist* that is *True* (which means the file must exist for the the
*CPdbDataFile* to be valid). A file like this can be created using the
*defEd* GUI.

When a script is run from CCP4i2 a *input_params* parameter file is
created by the GUI and for *pdbset* has the form (excluding the
meta-data header):

::

       <inputData>
         <XYZIN>
           <project>my_project</project>
           <baseName>broken.pdb</baseName>
           <relPath/>
         </XYZIN>
         <CELL>
           <a>64.897</a>
           <b>78.323</b>
           <c>38.792</c>
           <alpha>90.0</alpha>
           <beta>90.0</beta>
           <gamma>90.0</gamma>
         </CELL>
       </inputData>
       <outputData>
         <XYZOUT>
           <project></project>
           <baseName></baseName>
           <relPath></relPath>
         </XYZOUT>
       </outputData>
       <controlParameters/>

This file follows the structure specified in the DEF file with
containers *<inputData>*\ and *<outputData>*, and parameters XYZIN, CELL
and XYZOUT. Note that no value has been given for XYZOUT.

The actual script to run *pdbset* is in the file
`$CCP4I2/wrappers/pdbset/script/pdbset.py <../../wrappers/pdbset/script/pdbset.py>`__
and is:

::

   from core.CCP4PluginScript import CPluginScript
        
   class pdbset(CPluginScript):

       TASKMODULE = 'utility'
       TASKTITLE = 'PDBSet'
       TASKNAME = 'pdbset'
       TASKCOMMAND = 'pdbset'
       TASKVERSION= 0.0
       COMLINETEMPLATE = '''1 XYZIN $XYZIN XYZOUT $XYZOUT'''
       COMTEMPLATE = '''1 CELL $CELL.a $CELL.b $CELL.c $CELL.alpha $CELL.beta $CELL.gamma
   1 END'''

The *pdbset* class is a sub-class of *CPluginScript* and requires only
that some class attributes are set for the base class methods to do the
all the work. The class attribute TASKMODULE and TASKTITLE specify how
the script appears in the GUI and TASKVERSION is a version for the
script. Information to run the program is: TASKCOMMAND - the name of the
executable; COMLINETEMPLATE - a template for the command line;
COMTEMPLATE - a template for the program command file. The templates are
described in more detail `here <./command_templates.html>`__ but two
important point:

-  The first word on a line is evaluated (it could be a parameter or a
   short Python script) and if it is true then the line is written to
   the command file.
-  Any word beginning with $ is substituted by the value of the
   parameter.

The main *CPluginScript* methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *CPluginScript.__init__()* method creates a container (a
`CContainer <./data_classes.html#containers>`__ object) and reads the
script's DEF file to specify the allowed contents of the container; note
that this does not specify the actual parameter values. The parameter
values can be passed to the container programmatically (if the script is
been run from another pipeline script) or are read from a file. When a
script is run from CCP4i2 a PARAMS file (*input_params.xml*) is created
by the GUI and loaded by the *CPluginScript.__init__()* method.

The *CPluginScript.process()* method does most of the work. It calls the
following methods:

-  *checkInputData()* - checks that data in the *inputData* container is
   valid and will cause the script to fail if it is not so. Particularly
   this checks that the input files exist and return a list of invalid
   file parameters.
-  *checkOutputData()* - checks that output files have valid names.
   Usually the GUI does not provide output file names and this method
   will set appropriate file names. This method is expected to fix
   things and not cause the script to fail.
-  *processInputFiles()* - The base class method does nothing but many
   scripts will need to reimplement this method - see below.
-  *makeCommandAndScript()* - The base class method uses the
   COMLINETEMPLATE and COMTEMPLATE templates to generate the program
   command line and command file. As an alternative to providing these
   templates a script developer can reimplement the
   *makeCommandAndScript()* method.
-  *startProcess()* - uses the PROCESSMANAGER module to start the
   program executable specified by TASKCOMMAND in another process. The
   command line and file create by *makeCommandAndScript()* are
   used.This is reimplemented for classes that do not run an external
   process.

When the program run is finished the *postProcess()* method is called.
This calls three other methods:

-  *postProcessCheck()* - retrieves information from the PROCESSMANAGER
   to find out if the program process ran successfully. It returns
   CPluginScript.SUCCEEDED or CPluginScript.FAILED.
-  *processOutputFiles()* - The base class method does nothing but many
   scripts will need to reimplement this method - see below.
-  *reportStatus()* - this does three things:

   -  Saves the parameters from the container to a PARAMS file. These
      parameters will be similar to the *input_params* file, if it
      exists, but could include output file names added by
      *checkOutputData()* or any changes to parameters made by more
      sophisticated scripts.
   -  Reports the jobs completion to the database.
   -  Emits a signal to say it has finished that may be received by a
      controlling pipeline script.

Normally a wrapper may need to reimplement:

-  *processInputFiles()* to perform any manipulations on input data or
   files before calling the main program
-  *processOutputFiles()* to perform any manipulations on output data or
   files after the program has run
-  *makeCommandAndScript()* to create the command line and script that
   will run the program

All of these methods must return either *CPluginScript.SUCCEEDED* or
*CPluginScript.FAILED* which will determine whether the script will
continue running or will abort. There is a third possible return status
for *processOutputFiles()*: *CPluginScript.UNSATISFACTORY* for when the
job has run without error but not produced a useful result. This return
flag leads to an 'Unsatisifactory' database status for the job and the
job is presented differently in the GUI: the icon on the Job List
indicates it is unsatisfactory but the normal job report is presented
(rather than the generic failed job report.). The job report creation
code is initialised with the status flag 'Unsatisfactory' and should
provide an appropriate report.

For pipelines it is necessary to reimplement the *process()* method to
call the sub-tasks but it should not be necessary to re-implement
*process()* or *postProcess()* methods for wrappers. If you have a
wrapper issue that can not be addressed with the three methods above
please talk to Liz.

The CPluginScript process handling will handle the case of an externally
run program clearly failing (i.e. the Process Manager has a 'failed'
error code - see `the Process Manager appendix <#process_manager>`__ but
does not check the exit status of the external process as the output
from differnt programs is unpredictable. The task implementer can do
this in the *processOutputFiles()* method and should make an error
report before returning a CPluginScipt.FAILED:

::

   from core.CCP4Modules import PROCESSMANAGER

   class myPluginScript(CPluginScript):
     ERROR_CODES = { 201 : {'description' : 'My program failed with error code' } }

     ...

     def processOutputFiles(self):
       exitStatus = PROCESSMANAGER().getJobData(pid=self.getProcessId(),attribute='exitStatus')
       if exitStatus != 0:     #Beware what is expect from given program and os?
         self.appendErrorReport(201,'Exit status: '+str(exitStatus))
         return CPluginScipt.FAILED

       ....

See also the `error handling in
plugins <./error_handling.html#plugins>`__.

A Simple Pipeline Example - *demo_copycell*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a demo of the usual *mtzdump* and *pdbset* pipeline that uses
*mtzdump* to extract the cell parameters from an MTZ file and *pdbset*
to set these cell parameters in a PDB file. This would not be the
recommended way to achieve this result - there are much simpler Python
tools! To run this script from the gui tick the preferences check-box:
``Preferences -> Advanced -> Show testing/development tasks``
so the *Demo copy
cell* task is displayed in the *Test code for developers only* module at
the bottom of the *Task menu*. Note this script does not have a good
report generator so there is no report.

The code is in
`$CCP4I2/wrappers2/demo_copycell/script/demo_copycell.py <../../wrappers2/demo_copycell/script/demo_copycell.py>`__
and looks like this:

::

   class demo_copycell(CPluginScript):

       TASKMODULE = 'utility'
       TASKTITLE = 'Demo copy cell'
       TASKNAME = 'demo_copycell'
       TASKVERSION= 0.0
       ASYNCHRONOUS=True

       def process(self):
         
         # Check all input files exist
         nonExFiles = self.checkInputData()
         if len(nonExFiles)>0:
           self.reportStatus(CPluginScript.FAILED)
           return

         # Provide default output file names if necessary
         self.checkOutputData()

         # Create instance of mtzdump class (and get it registered with the database)
         self.mtzdump = self.makePluginObject('mtzdump')
         self.mtzdump.container.inputData.HKLIN = self.container.inputData.HKLIN
         # Run mtzdump to get the cell parameters
         self.connectSignal(self.mtzdump,'finished',self.process_1)
         self.mtzdump.process()

The TASKMODULE and TASKTITLE attributes are provided so the script can
be placed in the GUI but there are no TASKCOMMAND, COMTEMPLATE or
COMLINETEMPLATE attributes defined. The *process()* method calls
*checkInputData()* and *checkOutputData()*. *checkInputData()* returns a
list of names of non-existant input files. If this list has length
greater than zero then *reportStatus()* is called with a failed status
and the method returns. An instance of *mtzdump* is created using the
*makePluginObject()* method. The only data *mtzdump* needs is the input
MTZ file, HKLIN; this parameter value is copied from the *demo_copycell*
container to the *mtzdump* container. Before running *mtzdump.process()*
the *self.connectSignal()* method is called to ensure that when
*mtzdump* finishes (and emits a finished signal) the next method in the
pipeline (*process_1()*) is called.

When mtzdump finishes (by calling *postProcess()*) it emits a signal
with an additional status parameter that is either
*CPluginScript.SUCCEEDED* or *CPluginScript.FAILED*. The
*demo_copycell.process_1* has an argument for the status:

::

       def self.process_1(self,status):
         if status == CPluginScript.FAILED:
           self.reportStatus(status)
           return

         # Create an instance of pdbset and copy parameters
         self.pdbset = self.makePluginObject('pdbset')
         self.pdbset.container.inputData.XYZIN = self.container.inputData.XYZIN
         self.pdbset.container.inputData.CELL = self.mtzdump.container.outputData.CELL
         self.pdbset.container.outputData.XYZOUT = self.container.outputData.XYZOUT

         # run pdbset
         self.connectSignal(self.pdbset,'finished',self.postProcessWrapper)
         self.pdbset.process()

*process_1()* first checks the wrapper exit status and finishes if the
wrapper failed. It then creates an instance of the *pdbset* wrapper and
copies the appropriate parameter values to the pdbset container. It then
runs the *pdbset.process()* method after connecting the finish signal to
the *CPluginScript.postProcessWrapper()* method. *postProcessWrapper()*
accepts the exit status from the wrapper and finishes the pipeline
appropriately.

The above script runs asynchronously - after the *mtzdump* and pdbset
sub-tasks are started by calling *process()* it does not wait for then
to complete. This is discussed further `below <#running%20processes>`__.
The alternative approach of would lead to simpler code with only one
*process()* method:

::

       ASYNCHRONOUS=False

       def process(self):
         
         # Check all input files exist
         nonExFiles = self.checkInputData()
         if len(nonExFiles)>0:
           self.reportStatus(CPluginScript.FAILED)
           return

         # Provide default output file names if necessary
         self.checkOutputData()

         # Create instance of mtzdump class (and get it registered with the database)
         self.mtzdump = self.makePluginObject('mtzdump')
         self.mtzdump.container.inputData.HKLIN = self.container.inputData.HKLIN
         # Run mtzdump to get the cell parameters
         status = self.mtzdump.process()
         if status == CPluginScript.FAILED:
           self.reportStatus(status)
           return

         # Create an instance of pdbset and copy parameters
         self.pdbset = self.makePluginObject('pdbset')
         self.pdbset.container.inputData.XYZIN = self.container.inputData.XYZIN
         self.pdbset.container.inputData.CELL = self.mtzdump.container.outputData.CELL
         self.pdbset.container.outputData.XYZOUT = self.container.outputData.XYZOUT

         # run pdbset
         status = self.pdbset.process()
         if status == CPluginScript.FAILED:
           self.reportStatus(status)
           return

         return CPluginScript.SUCEEDED

Here the *ASYNCHRONOUS* attribute is set *False* to inform the core
system how to run this and there is not need for calls to
*connectSignal()*. The sub-task *process()* method returns a status
flag, either *FAILED* or *SUCCEEDED* and if this is *FAILED* then the
pipeline should either handle the problem or terminate itself by calling
*self.reportStatus(status)*. The *process()* method itself should return
*SUCCEEDED* or *FAILED*.

Pipeline and Wrapper File Organisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CCP4I2 directory has a *wrapper* sub-directory to which can be added
wrapper scripts. Each script has its own directory with *script* and
*test_data* sub-directories. The script directory should contain a
*wrapper*.py file containing the *CPluginScript* script and a
*wrapper*.def.xml file. There is a *test_data* directory for any data
needed to test the script. So the part of the *wrapper* directory for
the two demo wrappers:

::

      ccp4i2 - wrappers - mtzdump - script ---- mtzdump.py
                       |         |        |
                       |         |         ---- mtzdump.def.xml
                       |         |
                       |          - test_data - gere_nat.mtz
                       |                     |
                       |                      - mtzdump.xml
                       |
                        -- pdbset - script ---- pdbset.py
                                 |        |
                                 |         ---- pdbset.def.xml
                                 |
                                   - test_data - 1df7.pdb

There is a slightly more complicated directory structure for pipelines
which should allow for developers' more varied requirements when
developing a pipeline. The *defEd* utility will create the directory
structure for a new wrapper or pipeline (look under the *Tools* menu).

There is no technical difference between the *wrappers* and *pipelines*
directories but the understanding that developers are encouraged to use
code in *wrappers* and extend its functionality if necessary. There is
also a *wrappers2* directory for some unconventional wrappers used by
the core system.

Unittests for Wrappers and Pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unittests are useful for testing the non-graphical scripts particularly
while in development but `project based
testing <./project_based_testing.html>`__ is now preferred.

The CCP4i2 core code uses the standard Python *unittest* module for
testing (see the `Python
documentation <http://docs.python.org/release/2.6.6/library/unittest.html>`__).
There are currently a couple of issues to running a *CPluginScript*\ s
in *unittest*:

-  *CPluginScript* requires that Qt *QApplication* is created
-  The default mode for the PROCESSMANAGER running a process is
   asynchronously - this does not play nicely with *unittest* so it must
   be set to 'wait for finish' mode.

These issues can be addressed by using the *unittest.TestCase.setUp()*
and *unittest.TestCase.tearDown()* methods thus:

::

      def setUp(self):
       # Create a Qt application
       from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
       self.app = QTAPPLICATION()
       # make all background jobs wait for completion
       PROCESSMANAGER().setWaitForFinished(10000)

      def tearDown(self):
       from core.CCP4Modules import PROCESSMANAGER
       PROCESSMANAGER().setWaitForFinished(-1)

To ensure portability of the tests you need to avoid any fixed filepaths
in your tests. Please use a 'standard' temporary work directory:

::

      from core import CCP4Utils
      workDirectory = CCP4Utils.getTestTmpDir()

This will set *workDirectory* to your $HOME/ccp4i2_test. If this
directory does not exist then it will be created. If you do not like
this choice of directory then edit your local copy of CCP4Utils.

Typically the test code at the bottom of a pipeline module will look
like this:

::

   import unittest,os

   class test_my_pipeline(unittest.TestCase):

      def setUp(self):
       from core import CCP4Modules
       self.app = CCP4Modules.QTAPPLICATION()
       # make all background jobs wait for completion
       # this is essential for unittest to work
       CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

      def tearDown(self):
       from core import CCP4Modules
       CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

      def test_1(self):
        from core import CCP4Modules, CCP4Utils
        import os

        workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test1')
        if not os.path.exists(workDirectory): os.mkdir(workDirectory)

        self.wrapper = my_pipeline(parent=CCP4Modules.QTAPPLICATION(),name='test1',workDirectory=workDirectory)
        self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'pipelines','my_pipeline','test_data','test1.data.xml'))

        self.wrapper.setWaitForFinished(1000000)
        pid = self.wrapper.process()
        self.wrapper.setWaitForFinished(-1)
        if len(self.wrapper.errorReport)>0: print self.wrapper.errorReport.report()
        #test exit status correct
        exitStatus = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitStatus')
        self.assertEqual(exitStatus,0,'Process exit status non-zero: '+str(exitStatus))
        #test if output file created
        self.assertEqual(os.path.exists( str(self.wrapper.container.outputData.HKLOUT) ),1,'Failed to create copied pdb file '+ str(self.wrapper.container.outputData.HKLOUT))

   def TESTSUITE():
     suite = unittest.TestLoader().loadTestsFromTestCase(test_my_pipeline)
     return suite

   def testModule():
     suite = TESTSUITE()
     unittest.TextTestRunner(verbosity=2).run(suite)

The *unittest* sub-class *test_my_pipeline* has one or more methods
called *test_\** which each perform one test. They do this by
instantiating a plugin, loading suitable test data, calling the
*process()* method to run the plugin and then testing that the plugin
has run sucessfully. In this example the input data for *test_1* come
from a 'PARAMS' xml file (called ..my_pipeline/test_data/test1.data.xml
in the example code) that is loaded into the wrapper container by the
*loadDataFromXml()* method. Typically such a file could be set up by
using the GUI to run *my_pipeline* and then looking in the job directory
for the *input_params.xml* file and copying it (along with any input
data files) to the wrapper/pipeline *test_data* directory. The 'PARAMS'
xml file will probably need to be editted to set the appropriate paths
for input and output files. An input file will be in the CCP4I2 source
code tree and should have *project* set to *CCP4I2* (which will find the
the CCP4I2 root directory) and *relpath* set to the path to the
*test_data* directory:

::

       <HKLIN>
           <project>CCP4I2</project>
           <relPath>pipelines/my_pipeline/test_data</relPath>
           <baseName>whatsit.mtz</baseName>
       <HKLIN>

And for an output files:

::

       <HKLOUT>
           <project>CCP4I2_TEST</project>
           <baseName>whatsit.mtz</baseName>
       </HKLOUT>

CCP4I2_TEST is the 'temporary' directory returned by
CCP4Utils.getTestTmpDir().

The *unittest* documentation describes the test methods such as
*assertEqual()* and typically a plugin test could check the return codes
from the process and check that output files exist.

The sub-class of *unittest* (*test_my_pipeline* in the example above)
does the real work but in the code above there are two extra functions:
**TESTSUITE()** which should return the suite of tests and is used by
the CCP4i2 overall test mechanism to find the test suite for that
module. You could possibly have additional tests that are not returned
by *TESTSUITE()* but those that are returned are part of CCP4i2 release
tests. **testModule()** is a non-essential utility that can make running
the module test from the Python prompt very quick:

::

   > import my_pipeline
   > my_pipeline.testModule()

The script to run all the CCP4i2 module tests is *$CCP4I2/bin/runtest*
which will search all of the *wrappers* and *pipeline* directories for
*TESTSUITE()* functions and run them. The *runtest* script ignores
scripts which specify their *TASKMODULE* as *demo* or *test* but reports
any other module that does not have a *TESTSUITE()* function.

Testing and Debugging Scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a job fails when run from the gui the key diagnostic should be
presented in a report page; otherwise check in the *Project directory*
view on the gui where you can see the job directory and should look at
the *stout.txt*, *stderr.txt* and *diagnostic.xml* files. The
*diagnostic.xml* file lists all Python exceptions from the script.

It is possible to run scripts and external programs from the console.
After you have attempted to run a script from the gui there will be a
job directory in the project directory (*project_path/CCP4_JOBS/job_nnn*
where nnn is the job number listed in the gui) and this contains an
input params file (*input_params.xml*). To rerun the script with these
input parameters:

::

   > $CCP4I2_TOP/bin/runTask project_path/CCP4_JOBS/job_nnn/input_params.xml

Note that as this script is run it will update the database and progress
will be shown if the gui is open. Beware that rerunning jobs may have
the potential for confusing the database so avoid doing this in
important projects.

It is also possible to start the external process (i.e. running the
program) from a Python console. The stdout from the script will list
calls to the PROCESSMANAGER().startProcess() method that runs the
wrapped program as an external process. The stdout output is formatted
so that you can just cut-n-paste to rerun the external process from the
ccp4i2 Python prompt, *ccp4i2/bin/pyi2* (see `Console
mode <./introduction.html#console>`__). The command will look something
like:

::

   >>> PROCESSMANAGER().startProcess("cbuccaneer", ["-stdin"], "/Users/lizp/Desktop/test_projects/t3/CCP4_JOBS/job_127/com.txt", "/Users/lizp/Desktop/test_projects/t3/CCP4_JOBS/job_127/log.txt")

The arguments to startProcess() are:

**executable** - the executable to run

**args** - the command line arguments as a Python list of single words

**inputFile** - a command file

**logFile** - the output log file

Buccaneer - a realistic example wrapper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| See the `buccaneer
  script <../../wrappers/buccaneer/script/buccaneer.py>`__. This
  re-implements three methods (in the order that they are called):
| *processInputFiles()* which convert the input data files into a form
  suitable for the wrapped program particularly using *makeHklin()* to
  convert mini-MTZs to a single input MTZ
| *makeCommandAndScript()* creates the command line and command file for
  the program
| *processOutputFiles()* converts the files output by the program to the
  standard forms of CCP4i2. It should also add appropriate default
  *CDataFile.annotation* (i.e. the label for a file that is seen in the
  gui) and appropriate *contentFlag* and *subType* attributes to the
  output mini-MTZ data objects.

Note that it should be unnecessary to re-implement the more generic and
complex *process()* and *postProcess()* methods for wrappers. If you
have a wrapper issue that can not be addressed with these three methods
please talk to Liz. But pipelines do need the *process()* method
re-implementing.

Also note that all of these methods must return either
*CPluginScript.SUCCEEDED* or *CPluginScript.FAILED* which will determine
whether the script will continue running or will abort.

Handling Mini-MTZs
~~~~~~~~~~~~~~~~~~

Please read the user's introduction to mini-MTZs
`here <../general/miniMTZs.html>`__. Each mini-MTZ contains 1 to 4 data
columns to provide one consistent group of data and should never contain
anything else. The columns in the mini-MTZ always have the same
conventional names and types.

There is a CCP4XtalData.CMiniMtzDataFile sub-class of CCP4File.CDataFile
to handle mini-MTZs. This is further subclassed to *CObsDataFile*
(observed SFs or intensities), *CPhsDataFile* (phases represented as HL
or phi/fom), *CMapCoeffsDataFile* (F/phi map coeficients) and
*CFreeRDataFile*.

Two additional parameters have been added to the *CDataFile* class to
keep track of the content of miniMTZs:

-  *contentFlag* flags the representation and currently only applies to
   *CObsDataFile* and *CPhsDataFile* and the allowed values are listed
   in those classes (see below).
-  *subType* (not yet fully implemented) flags the provenance of the
   data (e.g. observed or calculated) and will be used to guide
   appropriate use of the data.

Note that other, non-MTZ files also have these parameters which might be
put to use some time in the future. There are also two additional
columns in the *Files* table of the database holding the same
information.

| These CMiniMtzDataFile classes also have several attributes that
  wrapper developers may need to access:
| *CONTENT_FLAG_xxxx* An allowed value for CMiniMtzDataFile.contentFlag
  parameter indicating a particular representation in the file
| *CONTENT_SIGNATURE_LIST* a list of the conventional column names
| *SUBTYPE_xxxx*\ An allowed value for CMiniMtzDataFile.subType
  parameter which has slightly different meaning for the different data
  types but generally indicates the provenance of the data.

CCP4 intends to move to using mini-MTZs but as few programs currently
support this way of working it is necessary for program wrappers to
convert multiple mini-MTZs to and from a single MTZ that is the usual
input to or output from a program. This should be done in the wrappers
re-implementation of *CPluginScript.processInputFiles()* and
*CPluginScript.processOutputFiles()*. There are two methods to help with
this conversion:

+----------------------+--------------------+-----------------------+
| CPluginScript method | Argument           | Description           |
+======================+====================+=======================+
| makeHklin            |                    | Merge two or more     |
|                      |                    | mini-MTZs to make an  |
|                      |                    | 'HKLIN' file          |
+----------------------+--------------------+-----------------------+
|                      | miniMtzsIn         | A list of either: the |
|                      |                    | names of              |
|                      |                    | *CMiniMtzDataFile*    |
|                      |                    | objects in the        |
|                      |                    | *c                    |
|                      |                    | ontainer.inputParams* |
|                      |                    | or a sub-list of the  |
|                      |                    | *CMiniMtzDataFile*    |
|                      |                    | name and the required |
|                      |                    | representation type   |
|                      |                    | (i.e. a               |
|                      |                    | *CONTENT_FLAG_xxxx*   |
|                      |                    | value).               |
+----------------------+--------------------+-----------------------+
|                      | hklin              | The basename of the   |
|                      |                    | output mtz - default  |
|                      |                    | 'hklin'               |
+----------------------+--------------------+-----------------------+
| splitHklout          |                    | Split an 'HKLOUT'     |
|                      |                    | file to multiple      |
|                      |                    | mini-MTZs             |
+----------------------+--------------------+-----------------------+
|                      | miniMtzsOut        | A list of the names   |
|                      |                    | of *CMiniMtzDataFile* |
|                      |                    | objects in the        |
|                      |                    | *co                   |
|                      |                    | ntainer.outputParams* |
+----------------------+--------------------+-----------------------+
|                      | programColumnNames | A list of lists of    |
|                      |                    | column names in the   |
|                      |                    | ''HKLOUT' file        |
|                      |                    | specifying the        |
|                      |                    | columns to copy to    |
|                      |                    | the files in          |
|                      |                    | *miniMtzsOut* list    |
+----------------------+--------------------+-----------------------+
|                      | inFile             | Name of the 'HKLOUT'  |
|                      |                    | file - defaults to    |
|                      |                    | *working              |
|                      |                    | Directory/hklout.mtz* |
+----------------------+--------------------+-----------------------+
|                      | logFile            | Name of the log file  |
|                      |                    | for program runs -    |
|                      |                    | defaults to           |
|                      |                    | *workingDi            |
|                      |                    | rectory/splitmtz.log* |
+----------------------+--------------------+-----------------------+
| splitHkloutList      |                    | Split a list of       |
|                      |                    | 'HKLOUT' files to     |
|                      |                    | mini-MTZs (see        |
|                      |                    | phaser_mr.py for      |
|                      |                    | example)              |
+----------------------+--------------------+-----------------------+
|                      | miniMtzsOut        | A list of the names   |
|                      |                    | of *CMiniMtzDataFile* |
|                      |                    | objects in the        |
|                      |                    | *co                   |
|                      |                    | ntainer.outputParams* |
+----------------------+--------------------+-----------------------+
|                      | programColumnNames | A list of lists of    |
|                      |                    | column names in the   |
|                      |                    | ''HKLOUT' file        |
|                      |                    | specifying the        |
|                      |                    | columns to copy to    |
|                      |                    | the files in          |
|                      |                    | *miniMtzsOut* list    |
+----------------------+--------------------+-----------------------+
|                      | outputBaseName     | base name for output  |
|                      |                    | files to which a      |
|                      |                    | numerical index and   |
|                      |                    | filetype              |
|                      |                    | specification will be |
|                      |                    | appended              |
+----------------------+--------------------+-----------------------+
|                      | outputContentFlags | a list mapping to     |
|                      |                    | splitHkloutList with  |
|                      |                    | the content flag      |
|                      |                    | parameters to be set  |
|                      |                    | in the output         |
|                      |                    | mini-MTZs             |
+----------------------+--------------------+-----------------------+
|                      | inFileList         | A list of the         |
|                      |                    | 'HKLOUT' files to be  |
|                      |                    | processed to miniMTZs |
+----------------------+--------------------+-----------------------+
|                      | logFile            | Name of the log file  |
|                      |                    | for program runs -    |
|                      |                    | defaults to           |
|                      |                    | *workingDi            |
|                      |                    | rectory/splitmtz.log* |
+----------------------+--------------------+-----------------------+

*CPluginScript* Class Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following class-wide attributes can be set in a *CPluginScript*
sub-class. Examples of their use are given above. These parameters will
usually be overruled by the same paramter if it appears in the gui
definition *CTaskWidget* sub-class. Note that it is also possible to
have the *TASKTITLE* defined in the *def* file (actually called
*PLUGINTITLE* as part of the file header) but this should not normally
be used except possibly if there is more than one *def* file associated
with the script and you need to distinguish them.

+----------------------------------+----------------------------------+
| **CPluginScript Class            | **Description**                  |
| Attributes**                     |                                  |
+----------------------------------+----------------------------------+
| TASKNAME                         | The obligatory name of the       |
|                                  | script essential for             |
|                                  | cross-refererence to other files |
+----------------------------------+----------------------------------+
| TASKVERSION                      | A version number for the script. |
|                                  | Please give one!                 |
+----------------------------------+----------------------------------+
| TASKMODULE                       | Refers to where the task should  |
|                                  | appear in the task menu, can be  |
|                                  | list of multiple modules see     |
|                                  | `Task                            |
|                                  | Men                              |
|                                  | u <./task_guis.html#taskmenu>`__ |
+----------------------------------+----------------------------------+
| TASKTITLE                        | The name for the task to appear  |
|                                  | in the GUI see `Task             |
|                                  | Men                              |
|                                  | u <./task_guis.html#taskmenu>`__ |
+----------------------------------+----------------------------------+
| DESCRIPTION                      | More details to appear in the    |
|                                  | GUI see `Task                    |
|                                  | Men                              |
|                                  | u <./task_guis.html#taskmenu>`__ |
+----------------------------------+----------------------------------+
| MAINTAINER                       | Email address of first-call      |
|                                  | maintainer of this code          |
+----------------------------------+----------------------------------+
| TASKCOMMAND                      | The name of the executable that  |
|                                  | the script will run              |
+----------------------------------+----------------------------------+
| SUBTASKS                         | List of main sub-tasks that      |
|                                  | might be run by a pipeline       |
|                                  | script - currently used for      |
|                                  | `bibliographic                   |
|                                  | references <./extras.h           |
|                                  | tml#bibliographic_references>`__ |
+----------------------------------+----------------------------------+
| COMTEMPLATE                      | A template for the program       |
|                                  | command file                     |
+----------------------------------+----------------------------------+
| COMTEMPLATEFILE                  | The name of a file containing a  |
|                                  | program command template         |
+----------------------------------+----------------------------------+
| COMLINETEMPLATE                  | A template for the program       |
|                                  | command line                     |
+----------------------------------+----------------------------------+
| INTERRUPTABLE                    | true/false can the script be     |
|                                  | interrupted?                     |
+----------------------------------+----------------------------------+
| RESTARTABLE                      | true/false - can the script be   |
|                                  | restatred                        |
+----------------------------------+----------------------------------+
| INTERRUPTLABEL                   | Text label to appear on button   |
|                                  | for interruptable sub-job        |
+----------------------------------+----------------------------------+
| DBOUTPUTDATA                     | A list of the output files to be |
|                                  | saved to database - leave        |
|                                  | undefined for all to be saved    |
+----------------------------------+----------------------------------+
| ASYNCHRONOUS                     | Boolean flag if the script, or   |
|                                  | any of its children, run         |
|                                  | external processes               |
|                                  | asynchronously                   |
+----------------------------------+----------------------------------+
| PERFORMANCECLASS                 | Name of class (in                |
|                                  | core/CCP4PerformanceData.py)     |
|                                  | that is used as the performance  |
|                                  | indicator for this plugin        |
+----------------------------------+----------------------------------+
| RUNEXTERNALPROCESS               | Boolean flag (default True)      |
|                                  | indicating that an external      |
|                                  | program will be run              |
+----------------------------------+----------------------------------+
| CLONEABLE                        | Boolean flag (default True)      |
|                                  | indicating if the task can be    |
|                                  | cloned                           |
+----------------------------------+----------------------------------+

| 

+------------------------------+--------------------------------------+
| **CPluginScript Attributes** | **Description**                      |
+------------------------------+--------------------------------------+
| command                      | The name of 'executable' that the    |
|                              | script will run.                     |
+------------------------------+--------------------------------------+
| container                    | A *CContainer* that contains the     |
|                              | parameters for a script. it is       |
|                              | initialised by the DEF file in the   |
|                              | plugin script directory              |
+------------------------------+--------------------------------------+
| workDirectory                | The directory for files produced by  |
|                              | the script                           |
+------------------------------+--------------------------------------+
| errorReport                  | A *CErrorReport* containing a log of |
|                              | script run                           |
+------------------------------+--------------------------------------+
| commandLine                  | A Python list of the words in the    |
|                              | command line to the 'executable'.    |
|                              | Should normally only be accessed     |
|                              | using methods below.                 |
+------------------------------+--------------------------------------+
| commandScript                | A Python list of strings - one item  |
|                              | for each line of a command file.     |
|                              | Should normally only be accessed     |
|                              | using methods below.                 |
+------------------------------+--------------------------------------+
| async                        | Flag if external processes should be |
|                              | run asynchronously. Default False.   |
+------------------------------+--------------------------------------+
| timeout                      | Timeout period for an external       |
|                              | process run in blocking mode.        |
+------------------------------+--------------------------------------+
| \_db\*                       | Parameters for interface to database |
|                              | and default directory and file       |
|                              | names. Do NOT mess with these!       |
+------------------------------+--------------------------------------+
| runningProcessId             | An i2-specific id for the last       |
|                              | external process started by this     |
|                              | plugin. This can be used to access   |
|                              | info in the *PROCESSMANAGER*         |
+------------------------------+--------------------------------------+
|                              |                                      |
+------------------------------+--------------------------------------+

| 
+----------------------+------------------+------------------------+
| CPluginScript method | Argument         | Description            |
+======================+==================+========================+
| \__init_\_           |                  |                        |
+----------------------+------------------+------------------------+
|                      | parent           | A parent class. Use    |
|                      |                  | *CCP4M                 |
|                      |                  | odule.QTAPPLICATION()* |
|                      |                  | if nothing else is     |
|                      |                  | appropriate.           |
+----------------------+------------------+------------------------+
|                      | name             | A unique name for the  |
|                      |                  | instance of the class. |
|                      |                  | If                     |
|                      |                  | *makePluginObject()*   |
|                      |                  | is used to instantiate |
|                      |                  | object then it will    |
|                      |                  | create a name          |
|                      |                  | dependent on the       |
|                      |                  | project and job id.    |
+----------------------+------------------+------------------------+
|                      | workDirectory    | A directory for        |
|                      |                  | 'temporary' files.     |
+----------------------+------------------+------------------------+
| makePluginObject     |                  | Create an instance of  |
|                      |                  | a CPluginScript        |
|                      |                  | sub-class - should be  |
|                      |                  | used rather than       |
|                      |                  | instantiating directly |
+----------------------+------------------+------------------------+
|                      | pluginName       | Name of                |
|                      |                  | wrapper/pipeline       |
|                      |                  | CPluginScript          |
|                      |                  | sub-class              |
+----------------------+------------------+------------------------+
|                      | reportToDatabase | Report the new plugin  |
|                      |                  | instance (i.e. job) to |
|                      |                  | the database - default |
|                      |                  | *True*                 |
+----------------------+------------------+------------------------+
| loadContentsFromXml  |                  | Set the contents of    |
|                      |                  | the container from the |
|                      |                  | DEF file in the plugin |
|                      |                  | directory.             |
+----------------------+------------------+------------------------+
| checkInputData       |                  | The base class         |
|                      |                  | implementation checks  |
|                      |                  | that all files         |
|                      |                  | specified in the       |
|                      |                  | *container.inputData*  |
|                      |                  | actually exist and     |
|                      |                  | returns a list of      |
|                      |                  | non-existant file      |
|                      |                  | (this is wrong - it    |
|                      |                  | should return a        |
|                      |                  | CErrorReport)          |
+----------------------+------------------+------------------------+
| checkOutputData      |                  | Provide output file    |
|                      |                  | names if missing       |
+----------------------+------------------+------------------------+
|                      | container        | *CContainer* to check  |
|                      |                  | - default *None* will  |
|                      |                  | check *self.container* |
+----------------------+------------------+------------------------+
| appendErrorReport    |                  | Append an error        |
|                      |                  | description to the     |
|                      |                  | overall plugin error   |
|                      |                  | report.                |
+----------------------+------------------+------------------------+
|                      | code             | An error code number   |
|                      |                  | unique for the class - |
|                      |                  | must be key in the     |
|                      |                  | class ERROR_CODES      |
|                      |                  | dictionary             |
+----------------------+------------------+------------------------+
|                      | details          | A Python string        |
|                      |                  | containing further     |
|                      |                  | information (e.g. name |
|                      |                  | of file that you could |
|                      |                  | not read). Default     |
|                      |                  | *None*                 |
+----------------------+------------------+------------------------+
|                      | name             | A Python string with   |
|                      |                  | an identifier for the  |
|                      |                  | data object.           |
|                      |                  | *CData.objectPath()*   |
|                      |                  | gives suitable value   |
|                      |                  | for this.Default       |
|                      |                  | *None*                 |
+----------------------+------------------+------------------------+
|                      | cls              | The class reporting    |
|                      |                  | the error. Default     |
|                      |                  | *None* will be         |
|                      |                  | interpreted as current |
|                      |                  | class.                 |
+----------------------+------------------+------------------------+
| extendErrorReport    |                  | Append errors from     |
|                      |                  | another *CErrorReport* |
|                      |                  | object to this object  |
+----------------------+------------------+------------------------+
|                      | other            | A *CErrorReport* or    |
|                      |                  | *CException*           |
+----------------------+------------------+------------------------+
| appendCommandLine    |                  | Append words to the    |
|                      |                  | command line to run an |
|                      |                  | 'executable'           |
+----------------------+------------------+------------------------+
|                      | \*args           | Arguments can be       |
|                      |                  | Python strings,        |
|                      |                  | *CData* object, or a   |
|                      |                  | list (Python list or   |
|                      |                  | *CList*) of these.     |
+----------------------+------------------+------------------------+
| appendCommandScript  |                  | Add a line to the      |
|                      |                  | command script         |
+----------------------+------------------+------------------------+
|                      | text             | A Python string or     |
|                      |                  | object with \__str__() |
|                      |                  | method or a list       |
|                      |                  | (Python or CList) of   |
|                      |                  | these.                 |
+----------------------+------------------+------------------------+
| writeCommandFile     |                  | Write text from the    |
|                      |                  | comand script to a     |
|                      |                  | file (the name is      |
|                      |                  | determined             |
|                      |                  | automatically).        |
+----------------------+------------------+------------------------+
| makeFileName         |                  | Generate an            |
|                      |                  | appropriate file name  |
|                      |                  | based on the work      |
|                      |                  | directory and unique   |
|                      |                  | plugin name.           |
+----------------------+------------------+------------------------+
|                      | format           | The file format -      |
|                      |                  | either                 |
|                      |                  | 'COM','LOG','PARAMS'   |
|                      |                  | or 'PRORAMXML'         |
+----------------------+------------------+------------------------+
|                      | ext              | Optional file          |
|                      |                  | extension - overrides  |
|                      |                  | the extension implied  |
|                      |                  | by the format          |
+----------------------+------------------+------------------------+
|                      | baseName         | Optional base filename |
|                      |                  | that overrides the     |
|                      |                  | plugin name.           |
+----------------------+------------------+------------------------+
| process              |                  | Default method to run  |
|                      |                  | an executable. Runs    |
|                      |                  | *checkInputData*,      |
|                      |                  | *makeCommandAndScript* |
|                      |                  | and *startProcess*     |
+----------------------+------------------+------------------------+
| makeCommandAndScript |                  | Creates command line   |
|                      |                  | and script from        |
|                      |                  | COMTEMPLATE or         |
|                      |                  | COMTEMPLATEFILE and    |
|                      |                  | COMLINETEMPLATE. Or    |
|                      |                  | may be reimplemented   |
|                      |                  | in sub-classes         |
+----------------------+------------------+------------------------+
| startProcess         |                  | Call                   |
|                      |                  | *PROCESSMAN            |
|                      |                  | AGER().startProcess()* |
|                      |                  | to start external      |
|                      |                  | process. Returns job   |
|                      |                  | status                 |
|                      |                  | C                      |
|                      |                  | PluginScript.SUCCEEDED |
|                      |                  | or                     |
|                      |                  | CPluginScript.FAILED.  |
+----------------------+------------------+------------------------+
|                      | command          | The external           |
|                      |                  | 'executable'.          |
+----------------------+------------------+------------------------+
|                      | reportStatus     | Boolean, default True. |
|                      |                  | If True startProcess   |
|                      |                  | calls reportStatus     |
|                      |                  | when complete.         |
+----------------------+------------------+------------------------+
| postProcess          |                  | Called after an        |
|                      |                  | external process       |
|                      |                  | started by             |
|                      |                  | *startProcess*         |
|                      |                  | finishes.Can be        |
|                      |                  | reimplemented but must |
|                      |                  | call *reportStatus()*  |
+----------------------+------------------+------------------------+
|                      | processId        |                        |
+----------------------+------------------+------------------------+
|                      | data             |                        |
+----------------------+------------------+------------------------+
| postProcessWrapper   |                  | Can be called at end   |
|                      |                  | of pipeline script     |
|                      |                  | instead of             |
|                      |                  | *postProcess* to clean |
|                      |                  | up after last wrapper  |
|                      |                  | plugin                 |
+----------------------+------------------+------------------------+
|                      | finishStatus     | Finish status of last  |
|                      |                  | wrapper -              |
|                      |                  | *CP                    |
|                      |                  | luginScript.SUCCEEDED* |
|                      |                  | or                     |
|                      |                  | *                      |
|                      |                  | CPluginScript.FAILED*. |
+----------------------+------------------+------------------------+
| reportStatus         |                  | Reports completion to  |
|                      |                  | database, writes       |
|                      |                  | PARAMS file, emits a   |
|                      |                  | 'finished' signal with |
|                      |                  | an job status.         |
+----------------------+------------------+------------------------+
|                      | finishStatus     | Finish status of job - |
|                      |                  | *CP                    |
|                      |                  | luginScript.SUCCEEDED* |
|                      |                  | or                     |
|                      |                  | *                      |
|                      |                  | CPluginScript.FAILED*. |
+----------------------+------------------+------------------------+
| postProcessCheck     |                  | Query the process      |
|                      |                  | manager to find out if |
|                      |                  | the external job       |
|                      |                  | succeeded. Returns     |
|                      |                  | *CPluginScript.FAILED* |
|                      |                  | or                     |
|                      |                  | *CP                    |
|                      |                  | luginScript.SUCCEEDED* |
+----------------------+------------------+------------------------+
|                      | processId        | The process id - as    |
|                      |                  | returned by            |
|                      |                  | *startProcess()*       |
+----------------------+------------------+------------------------+
| logFileText          |                  | Return the text of a   |
|                      |                  | log file created by    |
|                      |                  | the external process.  |
+----------------------+------------------+------------------------+
| updateJobStatus      |                  | Report job status to   |
|                      |                  | database - provide     |
|                      |                  | *status* or            |
|                      |                  | *finishStatus*         |
+----------------------+------------------+------------------------+
|                      | status           | A CCP4DbApi            |
|                      |                  | JOBSTATUS_\*           |
+----------------------+------------------+------------------------+
|                      | finishStatus     | *CP                    |
|                      |                  | luginScript.SUCCEEDED* |
|                      |                  | or                     |
|                      |                  | *CPluginScript.FAILED* |
+----------------------+------------------+------------------------+
| saveParams           |                  | Save the content of    |
|                      |                  | self.container to a    |
|                      |                  | file                   |
|                      |                  | *baseName.params.xml*  |
+----------------------+------------------+------------------------+
| getProcessId         |                  | Return the processId   |
|                      |                  | which can be used to   |
|                      |                  | query                  |
|                      |                  | PROCES                 |
|                      |                  | SMANAGER()getJobData() |
+----------------------+------------------+------------------------+
| joinMtz              |                  | Join multiple MTZs     |
|                      |                  | into one MTZ           |
+----------------------+------------------+------------------------+
|                      | outFile          | Full path name of the  |
|                      |                  | output file            |
+----------------------+------------------+------------------------+
|                      | inFiles          | A list, each element   |
|                      |                  | is a list of two       |
|                      |                  | elements: the full     |
|                      |                  | path name of a input   |
|                      |                  | file and a text string |
|                      |                  | of comma-separated     |
|                      |                  | column names.          |
+----------------------+------------------+------------------------+
| splitMtz             |                  | Split one MTZ into     |
|                      |                  | multiple MTZs          |
+----------------------+------------------+------------------------+
|                      | inFile           | Full path name of the  |
|                      |                  | input file             |
+----------------------+------------------+------------------------+
|                      | outFiles         | A list, each element   |
|                      |                  | is a list of three     |
|                      |                  | elements: the full     |
|                      |                  | path name of a output  |
|                      |                  | file, the input        |
|                      |                  | columns and the        |
|                      |                  | corresponding output   |
|                      |                  | columns. Thr input and |
|                      |                  | output columns are a   |
|                      |                  | text string of         |
|                      |                  | comma-separated column |
|                      |                  | names.                 |
+----------------------+------------------+------------------------+

Moving Data Between Plugins
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data is held in a *CContainer* which is described
`here <./data_classes.html#containers>`__. Pipeline scripts usually need
to copy data from their own container to to wrapper containers or
between wrapper containers. Single items of data can be copied simply by
equivalencing:

::

         self.pdbset.container.inputData.XYZIN = self.container.inputData.XYZIN

| If you need to copy more data between plugins then use
  *CContainer.copyData()* which copies data into the container. There
  are two arguments:
|    *otherContainer* - container to copy from
|    *dataList* - list of items to copy. If no list is provided it will
  copy all items with the same name (use with care!).
| For example:

::

      self.pdbset.container.inputData.copyData(self.container.inputData,['XYZIN'])

Asynchronous Pipelines
~~~~~~~~~~~~~~~~~~~~~~

The *CPluginScript* method uses *PROCESSMAMANGER().startProcess()* to
start the programs in a separate process. By default external processes
are run in blocking mode i.e. the script does not continue until after
the external program has completed. This is fine for single wrappers or
pipelines performing a series of consecutive steps but for some
pipelines you may want to have several program runs start concurrently
and then process the result for each run as it completes. There is a
trivial example of how to do this in the pipeline
`demo_multi_mtzdump <../../pipelines/demo_multi_mtzdump/script/demo_multi_mtzdump.py/>`__
that runs multiple instances of the *mtzdump* program that extracts data
from an MTZ file. The `mtzdump
wrapper <../../wrappers/mtzdump/script/mtzdump.py/>`__ is a simple
wrapper requiring no specific features to support this functionality.
The *demo_multi_mtzdump.process()* method makes a list of mtz files and
then creates one instance of the *mtzdump* wrapper for each file using
the usual *makePluginObject()* method. In this example the call to
*makePluginObject()* has the *reportToDatabase* argument set False to
prevent the individual tasks being saved to the database but this is not
at all a requirement for plugins that are run asynchronously. The name
of the input mtz is set in *mtzdump* object and then the code that is
specific for asynchronous running:

::

           # Set to run asynchronously and set a callback
           mtzdump.async = True
           mtzdump.setFinishHandler(self.handleDone)
           #Start process
           rv = mtzdump.process()
           #The mtzdump instance must be saved to keep it in scope and everything else can be got from that.
           self.subProcessList.append(mtzdump)

The *mtzdump.async* parameter is set *True* and callback is set to
*self.handleDone* a method defined later in the script. The *mtzdump* is
started by calling *mtzdump.process()* and the mtzdump object is saved
to a list; this is essential to ensure it does not go out of scope when
the *demo_multi_mtzdump.process()* method returns. The *handleDone()*
method is called as each mtzdump process finishes and is passed the
*jobId* and a *processId* from the completed run. The *jobId* will be
None if the script is run in a context without the CCP4i2 database or if
(as here) the *reportToDatabase* is False. So it is best, as in this
example, to use the pid to find which mtzdump instances has finished:

::

         # Get the exit status and if successful get the CELL from the outputData
         if self.postProcessCheck(pid) == CPluginScript.SUCCEEDED:
           for p in self.subProcessList:
             if p.processId == pid:
               print p.container.inputData.HKLIN,p.container.outputData.CELL
               self.subProcessList.remove(p)

In this trivial example a result is printed out and the finished mtzdump
instance is removed from the list of sub-processes. It is also very
important that this callback method tests for when all sub-processes are
finished and calls *CPluginScript.reportStatus()* with a suitable status
flag:

::

         if len(self.subProcessList)==0:
           self.reportStatus(CPluginScript.SUCCEEDED)

One other essential feature is that the *demo_multi_mtzdump* class has
an attribute *ASYNCHRONOUS* set to True at the top of the class
definition:

::

     class demo_multi_mtzdump(CPluginScript):
       TASKTITLE = 'Demo multi MTZ dump'
       TASKNAME = 'demo_multi_mtzdump'
       ASYNCHRONOUS = True

This is a flag to the infrastructure running the task to ensure that the
necessary event loops and threads are set up before creating an instance
of the class.

Connecting a pipeline with connectSignal()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An alternative approach to running pipelines asynchronously is after
creating the first subtask but before calling its *process()* method,
use *connectSignal()* to connect the 'finished' signal from that first
subtask to a target function that will start the second subtask. This is
shown below in an extract from *aimless_pipe.py* where the first
*pointless* plugin is created and then a *connectSignal()* call ensures
that the 'finish' signal will invoke *process_aimless()* method. The
*process_aimless()* method first checks the finish status of the
*pointless* subtask and then runs *aimless* and calls *connectSignal()*
to get *process_cycle_ctruncate()* called when the *aimless* plugin
finishes .. and so on.

One issue with using *connectSignal()* is that the target function can
crash without any exception being passed up to a calling function so
there is no evidence of what happened and no call to *reportStatus()* to
cleanup and ensure that the job is recognised to have failed. So for the
user the job appears to be still running. The way to handle this is to
wrap all of the functionality in the target method in a *try..except*
statement with the *except* body calling *appendErrorReport()* and
*self.reportStatus(CPluginScript.FAILED)*. This gets the exception shown
in the job's diagnostic and ensures that the pipeline is correctly
reported as failed. The *CPluginScript* error code 39 is a generic
'Unknown error'. It would be better to have several *try..except*
statements and to append more informative *CException* (ie ccp4i2
exception) to the pipeline error log using *appendErrorReport()*.

::

     def process(self):
         print "### Starting aimless pipeline ###"
         ...
         self.process_pointless()

       def process_pointless(self):
         print "### Running pointless ###"
         ...
         self.pointless = self.makePluginObject('pointless')
         ...
         self.connectSignal(self.pointless,'finished',self.process_aimless)
         self.pointless.process()

       def process_aimless(self,status):
         if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
            return
         try:
           print "### Running aimless ###"
           crash #Test the fail mechanism

           # Run aimless to edit model
           self.aimless = self.makePluginObject('aimless')
           ...
           self.connectSignal(self.aimless,'finished',self.process_cycle_ctruncate)
           self.aimless.process()
         except Exception as e:
           print e
           self.appendErrorReport(CPluginScript,39,str(e))
           self.reportStatus(CPluginScript.FAILED)

A 'Wrapper' that does not run an External Process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible that a wrapper might do the processing itself without
running another program in another external process. In this case the
developer should reimplement the *startProcess()* method to do the
necessary processing and should set the class parameter
*RUNEXTERNALPROCESS* to *False*. An example of a script which does this
is `mergeMtz <../../wrappers/mergeMtz/script/mergeMtz.py>`__ which just
calls the utility *CPluginScript.joinMtz()* to perform its
functionality. Note that *startProcess()* should return
*CPluginScipt.SUCCEEDED* or *CPluginScipt.FAILED*.

::

   class mergeMtz(CPluginScript):

       TASKTITLE = 'Merge experimental data objects to MTZ'     # A short title for gui menu
       TASKNAME = 'mergeMtz'                                    # Task name - should be same as class name
       RUNEXTERNALPROCESS = False                               # There is no external process

       def startProcess(self,command,**kw):

         '''
             Some data preparation
         '''
             
         status = self.joinMtz(self.container.outputData.HKLOUT.fullPath.__str__(),inFiles)
         return status

Interrupting and restarting processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an example see
*pipelines/buccaneer_build_refine/scripts/buccaneer_build_refine.py*
which will stop at the end of a buccaneer-refmac cycle if it receives an
interrupt signal. The pipeline can test if there is an interrupt signal
at a suitable break point by calling *CPluginScript.testForInterrupt()*
which will return True if an interrupt request has been made. At this
point the script should:

do whatever to terminate the script cleanly with meaningful results

if the script supports restart it should save any necessary parameters -
e.g. buccaneer_build_refine has an extra *interruptStatus* container in
the def file to hold these parameters.

call *self.reportStatus(CPluginScript.INTERRUPTED)* to register with the
database that the job is interrupted

*return* from the script

If the script supports restart then the start of the script should check
the saved parameters (e.g. in the *interruptStatus* container and set
the starting parameters appropriately. To notify the system that a
script supports interrupt or restart the class attributes
*INTERRUPTABLE* and *RESTARTABLE* should be set to True.

There is an alternative mechanism to support stopping a sub-job in a
pipeline so that the pipeline can continue. This is best seen in the
CRANK2 pipeline. The interruptable task should have the task attribute
INTERRUPTLABEL set to a text string that is a label to appear on button
underneath the task report. If the user clicks this button then there is
a signal that can be checked for with the
*CPluginScript.testForInterrupt()* method as described above. The
pipeline should have a mechanism to continue to the next step when it
gets this signal. To make this tool usable the job report should be
showing the progress of the sub-job so the user can make a sensible
decision when to interrupt it and the subsequent report should probably
contain a note that the job was interrupted.

The sub-job interrupt was implemented for the benefit of CRANK2 and has
some obvious limitations - if a task has INTERRUPTLABEL attribute then
it applies in all pipelines that run the task - subclassing the task
class could be used to get around this. The system currently only
suuports interrupting the immediate children of a pipeline task.

Performance Indicators
~~~~~~~~~~~~~~~~~~~~~~

A performance indicator is one or more simple parameters that give an
indication of the effectiveness of a job - for example R-factors and
free R-factors can be performance indicators for refinement. CCP4i2's
use of performance indicators is currently in development but the basic
approach is that a plugin script saves the key parameters in the
*outputData* section of the *def* file and they will be picked up from
here by the CCP4i2 system and saved to the database. These parameters
are then displayed in the project window. In future the parameters
stored in the database might help to guide automated structure solution.

The file
`core/CCP4PerformanceData.py <../../core/CCP4PerformanceData.py>`__
contains the definitions of performance indicator classes that are
sub-classed from *CPerformanceIndicator*. It is expected that there
should be one performance indicator class for each stage in structure
solution so tasks that perform a similar function should use the same
performance indicator class (the classes will provide interconversion of
different representations and whatever else is necessary to make this
work). The currently defined classes hold several parameters that are
all either floats or strings and the mechanism to save data to the
database currently only handles floats (or integers converted to floats)
or strings (255 character limit) but this convention could be changed if
necessary. The parameters in each class are up for discussion.

To use a performance indicator a plugin class should set the class
variable PERFORMANCECLASS to the name of the performance class used
(this is there to help in efficiently generating the display in the
project window). In the *def* file *outputData* section a variable
(usually called *PERFORMANCEINDICATOR* but that is optional) is declared
with the appropriate class. At the end of the process, the plugin script
should set the parameters in the performnce indicator object. To see an
example of this look for 'PERFORMANCE' in
`../../pipelines/bucref/script/buccaneer_build_refine.def.xml <../../pipelines/bucref/script/buccaneer_build_refine.def.xml>`__
and
`../../pipelines/bucref/script/buccaneer_build_refine.py <../../pipelines/bucref/script/buccaneer_build_refine.py>`__

Project Defaults
~~~~~~~~~~~~~~~~

This is a mechanism to set some parameter(s) for a given task for every
time that task is run within a project. Typically it might be used to
control restraints in Refmac runs when these need customising for the
particular coordinate data in the project. A task that wants to set some
default parameters (for subsequent runs of the same task or for another
task) creates a fragment of the *params* file for the task containing
only the required parameters and saves this file to the
*CCP4_PROJECT_FILES* sub-directory of the project directory. When the
task plugin is created the appropriate *default params* file is
automatically loaded.

The implementation of this feature is under development but presently
there are two *CPluginScript* methods to enable a developer to create or
edit a *default params* file. For example to set the *WEIGHT* parameter
for a Refmac wrapper.

::

      defaultsContainer = self.getProjectDefaultParameters(taskName='refmac_martin',paramsList=['WEIGHT'])
      if defaultsContainer.controlParameters.WEIGHT != 2.0:
        defaultsContainer.controlParameters.WEIGHT = 2.0
        self.saveProjectDefaultParameters(defaultsContainer)

Here *CPluginScript.getProjectDefaultParameters()* has been called with
the name of a task and a list of the parameters that will be set. This
method returns a 'top-level' container with the usual sub-containers
(*inputData*,\ *controlParameters* and *outputData*) but only the listed
parameters. Except if there is already a *default params* file for this
task then the contents of that file plus the listed parameters are
returned. You should check that the required defaults have not already
being set (this is to avoid unnecessary extra versions of the defaults
file), set the required defaults and then save the edited defaults using
*CPluginScript.saveProjectDefaultParameters()*. Note that the defaults
file is saved with the header containing the jobId and jobNumber of the
job that created it and the *default params* files are versioned with
file names of the form *CCP4_PROJECT_FILES/refmac_martin_001.params.xml*

The defaults file is read and loaded by
*CPluginScript.loadProjectDefaults()* called from
*CPluginScript.process()* (currently commented out - Dec 2014).

There are several issues with this mechanism - please discuss with Liz
before using it.

Incorporating non-CCP4i2 pipelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to get a pipeline not written with CCP4I2 *CPluginScript*
to appear in the i2 database and GUI. This could avoid writing wrappers
for every program called by the pipeline but does require following some
I2 conventions:

The pipeline 'top' script must be a CPluginScript

The pipeline must place its output files (data file, log files) in the
specified directory ( CPluginScript.workDirectory)

Each sub-process that you want to be visible to i2 must have output and
log files placed in the i2 specified directory

For each sub-process a dummy *CPluginScript* must by created and its
data container loaded with info on input and output data.

The *dummy* pipeline (*dumMee* in the *test* folder of the GUI) has a
script *$CCP4I2_TOP/pipelines/dummy/script/dummy.py* that demonstates
how to do this.

The CRANK2 pipeline is an example of this.

.. _testing-and-debugging-scripts-1:

Testing and Debugging Scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to run scripts and external programs from the console.
After you have attempted to run a script from the gui there will be a
job directory in the project directory (*project_path/CCP4_JOBS/job_nnn*
where nnn is the job number listed in the gui) and this contains an
input params file (*input_params.xml*). To rerun the script with these
input parameters:

::

   > $CCP4I2_TOP/bin/Python $CCP4I2_TOP/bin/runTask.py >project_path/CCP4_JOBS/job_nnn/input_params.xml

Exporting a task as a compressed file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a mechanism to create a compressed file (currently a zip) to
make available to a limited number of users as an alternative to
checking code into the normal CCP4i2 repository for general
distribution. The code to export must be a sub-directory of the
*pipelines* directory and should be cleaned of *.pyc* and other
superfluous files. The *Export task* option is on the *Developer tools*
sub-menu of the *Utilities* menu. You will be asked to choose the
directory to export and the compressed file will be placed in the
*pipelines* directory. The user can import the task via the *Import
task* option on the *Utitlities* menu.

Exporting files from a job
~~~~~~~~~~~~~~~~~~~~~~~~~~

| In the project window job list the job context menu can include an
  option to export files from the job. This is likely to be particularly
  useful for exporting 'monster' MTZ files for use outside CCP4i2 but it
  could work for any file that either exists in the job directory or
  sub-directory or can be created from the data there. The task script
  module should provide two functions (note functions not methods of the
  CPluginScript class):
| **exportJobFileMenu( jobId )** - is called when the user opens the job
  context menu and should return details of items to appear on an
  'Export' menu. The function is passed the jobId of the job chosen by
  the user and should return a list where each item is a sub-list of
  three strings:
|     an identifier that will be *mode* input to *exportJobFile()*
|     the text string that will appear on the *Export* sub-menu - should
  inform the user exactly what will be exported
|     the mime type of the exported file (see CCP4CustomMimeTypes) -
  used by the file selection dialog for file extensions etc.
| **exportJobFile( jobId , mode)** - is called when the user selects one
  of the options from the *Export* sub-menu. The function is passed the
  jobId of the selected job and the identifier from the sub-menu item.
  The function should return the filename of the file to export - it may
  need to create a temporary file. There is an example of this in the
  aimless_pipe script. A code fragment that might be useful in this
  context:

::

       from core import CCP4Modules
       childJobs =CCP4Modules.PROJECTSMANAGER().db().getChildJobs(jobId=jobId,details=True,decendents=False)
       for jobNo,jobId,taskName  in childJobs:
         if taskName == 'ctruncate':
           f = os.path.join( CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False),'HKLOUT.mtz')
           if os.path.exists(f): return f

Here the *CDbApi.getChildJobs()* returns a list of the subJobs for the
*jobId* with the job number, jobId and taskname for each sub-job
returned. The *CProjectsMaanger.jobDirectory()* returns the job
directory for any job or sub-job.

| A likely requirement to create suitable export files is merging two or
  more MTZs. The usual CCP4i2 approach using *cmtzjoin* expects to work
  with column groups corresponding to i2 data objects and may not work
  for the 'old style' MTZs that are being exported. There is a simple
  interface to the CAD program via *CMtzDataFile.runCad()*. This has
  arguments:
|   *hklout* - output file name
|   *hklinList* - list of mtz to be merged with the file referenced by
  *CMtzDataFile*
|   *comLines* - list of strings to go into command file (END not
  necessary)
| If *comLines* is zero length then the compulsary LABIN input will be
  set to 'ALL' for all input files. If you require anything other than
  all the files merging you should set *comLines* to the necessary CAD
  command lines. This return hklout(str),err(CErrorReport).

The *hklout* should be set to something like 'export.mtz' in the job
directory and *runCad()* will create this file and *export_cad.com* and
*export_cad.log* in the job directory.

Example:

::

   def exportJobFile(jobId=None,mode=None):
       import os
       from core import CCP4Modules
       from core import CCP4XtalData

       jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId,create=False)
       exportFile = os.path.join(jobDir,'exportMtz.mtz')
       if os.path.exists(exportFile): return exportFile

       ...
       
       m = CCP4XtalData.CMtzDataFile(truncateOut)
       #print m.runCad.__doc__   #Print out docs for the function
       outfile,err = m.runCad(exportFile,[ freerflagOut ] )
       return   outfile

Temporary files and file cleanup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Considerations:
| Ideally wrappers and task guis should put intermediate files in the
  appropriate job directory so they can be found and (potentially)
  bundled in when a job or project is exported. Or found and deleted.
| There is a need to limit the size of job and project directories since
  users may want to export them and compressing and moving large files
  is slow.
| Many programs create large intermediate files. Many pipelines have
  little used intermediate data files.

In order to limit disk usage and user information overload some 'scratch
files' or little used files are deleted. At the end of a job run the
core system automatically applies a file cleanup and there are user
interface options to cleanup files from an individual job or whole
project. To support this the core system keeps a list of files that can
be deleted and individual tasks should have a list of files specific to
the task. The *CCP4ProjectsManager.CPurgeProject* class handles this and
the parameters shown below come from that class. The multi-step process
of assigning files to a catagory and then deleting certain catagories in
a given context allows for user customisation and the inevitable change
of requirements.

The deletable files are assigned to a catagory:

::

     PURGECODES = { 1 : 'Scratch file',
                    2 : 'Scratch file potentially useful diagnotic',
                    3 : 'Diagnostic file',
                    4 : 'File redundant after report created',
                    5 : 'Intermediate data file',
                    6 : 'Redundant on project completion',
                    7 : 'Large file that could be deleted',
                    0 : 'Retained file - used to override default' }

In any context the catagories of file that can be deleted are specified
(e.g. [1,2] on completion of a job). The *CPurgeProject* class has a
list of the contexts:

::

     CONTEXTLOOKUP = { 'script_finish' : [1, 2],
                       'script_finish_fail' : [1],
                       'script_finish_remote' : [1, 2],
                       'temporary' : [1, 2, 4],
                       'intermediate' : [1, 2, 4, 5],
                       'extended_intermediate' : [1, 2, 5, 7],
                       'project_complete' : [1, 2, 3, 4, 5, 6] }

The *CPurgeProject* class also has
a default list of files and purge catagories that specify files created
by many or all tasks. This currently is:

::

     SEARCHLIST = [ [ '*-coordinates_*.*', 1 ],
                    [ 'report.previous_*.html' , 1 ],
                    [ 'report_tmp.previous_*.html' , 1 ],
                    [ 'params.previous_*.xml' , 1 ],
                    [ 'diagnostic.previous_*.xml' , 1 ],
                    [ 'hklin.mtz' , 1 ],
                    [ 'hklout.mtz', 1 ],
                    [ '*mtzsplit.log', 2 ],
                    [ 'log_mtz*.txt', 2 ],
                    [ 'sftools', 2 ],
                    [ 'ctruncate', 2 ],
                    [ 'scratch', 2 ],
                    [ 'stderr.txt' , 3],
                    [ 'stdout.txt' , 3],
                    [ 'diagnostic.xml' , 3],
                    [ 'report_tmp.html' , 4 ],
                    [ 'program.xml' , 4 ],
                    [ 'XMLOUT.xml' , 4 ],
                    [ 'log.txt', 4 ],
                    [ '*.scene.xml' , 6 ],
                    [ '*observed_data_as*', 6 ],
                    [ '*phases_as*', 6 ],
                    [ '*%*/report.html', 6 ],
                    [ '*%*/tables_as_csv_files', 6 ],
                    [ '*%*/tables_as_xml_files', 6 ],
                    [ '*%*/params.xml', 6 ],
                    [ 'com.txt', 6 ] ]

So, for example, *hklin.mtz* is catagory 1, a scratch file that is
liable to be deleted in many context. But many tasks have either
addition files suitable for deletion or perhaps some that require saving
despite appearing on the above list. Particularly pipelines might have
sub-job files suitable for deletion. The task-specific files should be
listed as an attribute of the plugin class - for example for the
*import_merged* task:

::

       PURGESEARCHLIST = [ [ 'HKLIN*.mtz' , 1 ],
                           ['aimless_pipe%*/HKLOUT*.mtz', 1]
                          ]

*PURGESEARCHLIST* is a list of filename and purge catagory pairs. Here
intermediate input mtz files and the output from *aimless_pipe* (run
just for the analysis) are catagorised at scratch files (catagory 1).
Note that '*' can be used as a wildcard in filenames and that files in
sub-jobs are specified by *taskname*\ %\ *subjob_number*/*filename*. The
*taskname* or *subjob_number* can have '*' as wildcard. The sub-jobs can
be nested. A more complex example is for *buccaneer_build_refine_mr*:

::

       PURGESEARCHLIST = [ [ 'buccaneer_mr%*/XYZOUT.*' , 5 , 'XYZOUT' ],
                         [ 'refmac%*/ABCDOUT.mtz', 5, 'ABCDOUT' ],
                         [ 'refmac%*/FPHIOUT.mtz', 5, 'FPHIOUT' ],
                         [ 'refmac%*/DIFFPHIOUT.mtz', 5, 'DIFFPHIOUT' ]
        ....
                         [ 'prosmart_refmac%*/refmac%*/XYZOUT.pdb', 5, 'XYZOUT' ],
                         [ 'prosmart_refmac%*/refmac%*/COOTSCRIPTOUT.scm', 5, 'COOTSCRIPTOUT' ] ]

Here the files being catagorised are the intermediate sub-job data files
from *buccaneer_mr* and *refmac* subtasks which by default the sub-tasks
record in the project database. When these files are deleted they must
also be removed from the database and the list above also gives the
sub-job parameter name for each file so it can be found in the database
and removed.

The purge mechanism builds a search list for a job by first creating a
search list for sub-jobs composed of the program-wide default list
overridden by the specific purge list for the sub-job task and then
overidden by the parent task purge list. Then, for all items on the
search list in the right catagories for the context, the mechanism looks
for files matching the filename and deletes them.

Appendix
~~~~~~~~

Parents
'''''''

For Qt signals to work *CPluginScript* is sub-classed from a *QObject*.
All *QObject*\ s must have a parent which is another *QObject* or the
*QApplication* (the top-level Qt object). If you are instantiating a
pipeline with no obvious parent, such as in the unittest code, then the
function *CCP4Modules.QTAPPLICATION()* can be used to get the
*QApplication* that can be used as a parent.

File Organisation
'''''''''''''''''

The objective is to have everything for one project in one place and to
allow the pipeline developer as much flexibility as possible but to have
some minimal organisation so 'the system' can find the 'plugins' and
other developers can use the resources.

| I suggest the following (but this is totally up for discussion):
| The wrappers for individual programs are in the *wrappers* directory
  where they are a resource for all developers and we avoid duplication
  of effort.
| The 'pipelines' go in a separate *pipelines* directory that should
  contain sub-directories such as 'MrBump', 'CRANK' etc. Each of these
  projects contain the sub-directories:
| *script*: analogous to the *wrappers* directory this should contain,
  for example, a 'MrBump.py' file containing the MrBump class that is a
  sub-class of *CPluginScript*. The *script* directory could contain
  more than one script module - it is likely one pipeline development
  group might contribute more than one top level script.
| *utils* containing utility code that is part of the project but
  hopefully is available to other developers
| *test_data* please provide all test data (er.. within reason)
| *gui* containing (by analogy with the *tasks* directory) a file
  CTaskMrBump.py containing a CTaskMrBump class that is a sub-class of
  *CTaskWidget*.
| *docs* containing documentation

There is an alternative approach to this organisation: the project could
provide a function or perhaps a 'resources' file that lists what
resources (scripts, guis, docs) it wants to have 'published' in the
overall ccp4i2 system. But, if it is not too onerous on pipeline
developers, I think some standard organisation would help us to work
together.

The Process Manager
'''''''''''''''''''

The process manager currently just starts jobs on the local mchine but
is intended in future to be the interface to access batch queues and
remote resources etc.. There is one instance of a process manager for
each Python process but note that if that process is a pipeline there
may be many different plugins (instances of *CPluginScript*) run in that
one Python process. The developer implementing a plugin should not
normally need to access the process manager directly but the method
*PROCESSMANAGER().getJobData()* may be used to retrieve information on a
running job:

::

   from core import CCP4Modules
   status = CCP4Modules.PROCESSMANAGER().getJobData(pid=self.runningProcessId,attribute='exitStatus')

The pid is an integer and is i2's internal unique id for a process.
*CPluginScript.runningProcessId* is the pid of the last process run by
the plugin. The attributes are listed below - a value of None will be
returned if the attribute is not set.

+-------------+-------------------------+--------------------------+
| command     | string                  | The name of the          |
|             |                         | executable               |
+-------------+-------------------------+--------------------------+
| arglist     | list of strings         | The arguments to the     |
|             |                         | external program         |
+-------------+-------------------------+--------------------------+
| handler     | Python command          | Callback to be run after |
|             |                         | external process         |
|             |                         | completed                |
+-------------+-------------------------+--------------------------+
| jobId       | string                  | The i2 database job id   |
|             |                         | (if provided)            |
+-------------+-------------------------+--------------------------+
| jobNumber   | string                  | The i2 job number in the |
|             |                         | form n[.m[.i]] (if       |
|             |                         | provided)                |
+-------------+-------------------------+--------------------------+
| ifAsync     | Boolean                 | Flag if process is to    |
|             |                         | run asynchronously       |
+-------------+-------------------------+--------------------------+
| timeout     | float                   | The number of seconds to |
|             |                         | wait for completion of a |
|             |                         | job run in blocking mode |
+-------------+-------------------------+--------------------------+
| editEnv     | list of list of strings | A list of edits to       |
|             |                         | subprocess environment   |
+-------------+-------------------------+--------------------------+
| inputFile   | string                  | path to input command    |
|             |                         | file                     |
+-------------+-------------------------+--------------------------+
| logFile     | string                  | path to ouput log file   |
+-------------+-------------------------+--------------------------+
| exitStatus  | integer                 | 0=OK,1=crashed           |
+-------------+-------------------------+--------------------------+
| exitCode    | integer                 | exit code of the last    |
|             |                         | process finished         |
+-------------+-------------------------+--------------------------+
| errorReport | CErrorReport            | list of exceptions       |
|             |                         | caught in running the    |
|             |                         | process                  |
+-------------+-------------------------+--------------------------+
| startTime   | float                   | The machine processor    |
|             |                         | time that process        |
|             |                         | started                  |
+-------------+-------------------------+--------------------------+
| finishTime  | float                   | The machine processor    |
|             |                         | time that process        |
|             |                         | finished                 |
+-------------+-------------------------+--------------------------+
| qprocess    | Qt object               | The QProcess that        |
|             |                         | supported the external   |
|             |                         | process (if this is used |
|             |                         | rather than subprocess)  |
+-------------+-------------------------+--------------------------+

By default external processes are run in blocking mode,
(*PROCESSMANAGER().startProcess()* does not return until the external
process has completed. The *ifAsync* flag is set from the
*CPluginScript.async* flag and if that is True then
*PROCESSMANAGER().startProcess()* will return immediately and the
pipieline could possible start more jobs. The Process Manager will only
start a limited number of asyncronous jobs at one time and will queue
any excess requests to start when other processes finish. The maximum
number of processes is set in the file
*$HOME/.CCP4I2/configs/ccp4i2_configs.params.xml* as
*maxRunningprocesses* parameter. (If the parameter is not in your file
then delete the file so a new one is created.) Note that if the user
starts multiple tasks then there could be many more running processes.

The *editEnv* argument is a list of edits applied to the current
environemnt. Each item in the list is a list. The first item in the
sub-list is an environment variable name. If the sub-list is length 1
(just the environment variable name) then that environment variable is
deleted; if the sublist is length 2 then the environment variable is set
to the second item in the list.

XML Files
'''''''''

All of the XML files used by CCP4i2 have the same basic structure, for
example a parameter file..

::

   <?xml version='1.0' encoding='ASCII'?>
   <ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
     <ccp4i2_header>
       <function>PARAMS</function>
       <ccp4iVersion>0.0.5<ccp4iVersion>
       <comment/>
       <hostName>bilbo.chem.york.ac.uk</hostName>
       <userId>lizp</userId>
       <creationTime>11:06 20/Jun/14</creationTime>
       <projectName>whokshop2</projectName>
       <jobId>82d6e44af86211e39b3a3c0754185dfb</jobId>
       <jobNumber>60</jobNumber>
       <projectId>3a3041bdf16b11e38aeb3c0754185dfb</projectId>
       <pluginName>prosmart_refmac</pluginName>
       <pluginTitle>Refine with Prosmart-restrained Refmac</pluginTitle>
     </ccp4i2_header>
     <ccp4i2_body>
       <inputData>
         <XYZIN>
             ....
     </ccp4i2_body>
   </ccp4:ccp4i2>

Here the root element is *ccp4i2* with an *xmlns* attribute that gives
the namespace as "http://www.ccp4.ac.uk/ccp4ns". The *ccp4i2* element
has two sub-elements: *ccp4i2_header* which contains metadata and
*ccp4i2_body* whose contents will depend of the type of file which is
specified by the *function* element in the header. The header usually
contains details of the program version and source of the file
(host,user, creation time) and, where appropriate, details of the
project and job that the file is associated with. The header is usually
written or read automatically by the *CCP4File.CI2XmlHeader* class and
developers should have little need to work directly with it.

There are some tools for handling XML files in CCP4i2. The Python
library for *lxml* is included and has many tools which are described at
`lxml website <http://lxml.de/index.html>`__. The internal
representation of XML is as an `etree.Element <lxml.de/tutorial.html>`__
for which there is brief (but likely enough for ccp4i2 developers)
documentation within the Python documentation
`here <http://docs.python.org/2/library/xml.etree.elementtree.html#module-xml.etree.ElementTree>`__
(scroll down for Element). See also the comments in the CCP4i2
documentation for `Reports <./reports.html#xml>`__. The *etree*
functionality is used to build up an internal representation of XML
which can be written to a file with the *CCP4Utils.saveEtreeToFile()*
module. An example to create a 'program' XML file in a CPluginScript:

::

         from lxml import etree
         from core import CCP4Utils
         # Create a 'root' element and add a 'mode' subElement
         rootEle = etree.Element('MYPROGRAM')
         modeEle = etree.SubElement(rootEle,'mode')
         modeEle.text = myMode

         # Load data from another XML file and append it to our tree
         try:
           otherXmlTree = CCP4Utils.openFileToEtree(..other-file-name..)
         except:
           # Report failure to read file
           self.appendErrorReport(110,other-file-name,stack=False)
         else:
           rootEle.append(otherXmlTree)

         # Save the xml tree to the 'PROGRAMXML' file
         CCP4Utils.saveEtreeToFile(rootEle,self.makeFileName('PROGRAMXML'))

--------------

Last modified: Tue Sep 13 14:03:56 BST 2016
