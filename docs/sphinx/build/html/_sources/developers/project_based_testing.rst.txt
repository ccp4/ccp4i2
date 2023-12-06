########################################
CCP4I2 Developers: Project Based Testing
########################################

Introduction
~~~~~~~~~~~~

Project based testing provides automatic testing of the non-graphical
aspects of tasks and also provide test data and indications of the
normal use for manual testing of the gui. Ideally we need one
test-project per task that exercises all the main functionality of task.
Some test projects that run through a series of tasks in a realistic
structure solution would also be useful for testing and as examples for
demonstrations and tutorials.

The automatic mechanism tests plugin scripts by rerunning all of the
jobs from a source project and comparing the output from the test jobs
with that from the source jobs. The source project can be set up
graphically in CCP4i2 but the test procedure is normally run
non-graphically by the *testi2sys* script. When the test jobs are run
they take their input files from the source project directory even where
the job input could be taken from a preceeding job in the test project
(but sub-jobs in pipelines do pass files through the pipeline in the
expected fashion).

After a test job has finished an automatic comparison is made of the
outputs. For each job run the test system calls
*CPluginScript.assertSame()* to do a comparison of the output data of
the source and the test job as taken from the *params.xml* file of the
two job runs. This method iterates over the data objects in the
*outputData* container and calls the *assertSame()* method passing it
the 'other' job run data value. This method returns a *CErrorReport*
which will contain a list of difference between the two values. Only the
*assertSame()* for file objects and performance indicator objects do
meaningful tests.

There is usually one performance indicator parameter per task whose main
function is a simple indication of progress shown in the *Job list* (see
`here <./pipelines.html#performance_indicators>`__). The same parameter
turns out to be useful for comparison in testing. The implementations of
*CPerformanceIndicator.assertSame()* usually compares the float,int or
string members of the class - there is an allowed tolerance and it can
be specified that a changed value may be an improvement (e.g. a decrease
in test Rfactor compared to source is reported as a warning not an
error). It is possible to implement a performance indicator class so
that it serves as a test class but does not put anything in the database
or *Job list* (eg. *CAtomCountPerformance*).

The test on output files is that they exist and then dependent on file
type:

+--------------+------------------------------------------------------+
| CDataFile    | compare file checksum and/or optionally file size    |
+--------------+------------------------------------------------------+
| CPdbDataFile | compare file checksum with *ad hoc* editing out of   |
|              | problem meta data                                    |
+--------------+------------------------------------------------------+
| CMtzDataFile | compare cell,spaceGroup,resolutionRange, number of   |
|              | reflections and column types/labels                  |
+--------------+------------------------------------------------------+

File comparison is tricky - for example a checksum test on PDB files
output by Molrep failed because of the date in a HEADER line - this has
been adressed but more, similar gotchas are expected. This is a work in
progress - please talk to Liz about any false calls and suggestions.

The *com.txt* input file for the program run by each task is also
compared.

The results are listed to the terminal and saved to a file in the test
project called *CCP4_test-*\ year-month-day-hour-minute\ *.log*.

Creating Source Projects
~~~~~~~~~~~~~~~~~~~~~~~~

A source project can be created as a normal project within CCP4i2. All
jobs should be reasonably quick to run and complete successfully. The
source project is exported by the *Export* option in the *Manage
projects* window (accessed from the *Projects* menu). This creates a zip
file which can be sent to Liz for inclusion in the automatic test suite.
Before exporting the test project please test it from the *Manage
projects* window by selecting the appropriate project and using the
*Rerun test project* option - you will be required to specify a
directory for output. Note that *Rerun test project* is only available
if you are running in *DEVELOPER* mode (if necessary edit
*.CCP4I2/configs/ccp4i2_cofig.params.xml* and set *developer* to
*True*). A test directory will be created and the test jobs will then
proceed automatically one at a time. The results from the comparison of
the source tasks and the test tasks are shown in real time in a web
browser window. You can also open the project window for the test
project to view progress. The test results are also saved to a file in
the test project directory called
*CCP4_test-*\ year-month-day-hour-minute\ *.log*. The log file should
list each job that is run (only the top level jobs) and the output data
from the job looking like:

::

   Job: 1 molrep_mr
   XYZOUT               -OK- CPdbDataFile:300 Compared objects are the same
   Job: 35 buccaneer_build_refine
   XYZOUT               -OK- CPdbDataFile:300 Compared objects are the same
   FPHIOUT              -ERROR- CMapCoeffsDataFile:308 Files failed checksum comparison/Users/lizp/Desktop/testing/demo8/CCP4_JOBS/job_35/FPHIOUT.mtz : /Users/lizp/Desktop/testing/demo8_rerun/CCP4_JOBS/job_35/FPHIOUT.mtz
   ...
   ...

The items on the line for each item of data are: the parameter name, the
test status (OK,WARNING or ERROR), the data class, the error code (error
code 300 is 'OK') and further details. If this test run fails to run of
gives false fails please discuss with Liz.

Note that if you look at the test project in the *Project Viewer* the
imported files will be absent.

Running Tests Non-Graphically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A test project or a suite of test projects can be run using the
*ccp4i2/bin/testi2sys* script. This script works with a temporary
database and unpacks the source projects into the *outputDirectory* and
loads them into the database. It then reruns the source projects
creating a test project. Each test project is given the name of the
source project with the date and time and '_rerun' appended. The test
project directories created in the *outputDirectory* are given the test
project name. The jobs from the source project are then run one at a
time with progress and errors listed to the terminal and saved to file
*CCP4_test.....log* in the *outputDirectory*. The command line
arguments:

+---------------------+-----------------------------------------------+
| -p *path*           | path of a compressed source project file or a |
|                     | directory containing them                     |
+---------------------+-----------------------------------------------+
| -d *directory*      | directory where the compressed source project |
|                     | is unpacked and temporary database and test   |
|                     | projects are created                          |
+---------------------+-----------------------------------------------+
| -o                  | Overwrite any existing temporary database     |
+---------------------+-----------------------------------------------+
| -runMode *flag*     | 0 :no testing, 1 :run jobs, 2:run tests,      |
|                     | 3:run jobs and tests. Default: 3              |
+---------------------+-----------------------------------------------+
| -verbosity *flag*   | Report 0:fails only, 1:jobs and parameters    |
|                     | tested, 3:jobs, parameters and parameter      |
|                     | values. Default: 1                            |
+---------------------+-----------------------------------------------+
| -testSubJobs *flag* | 0:no, 1:for failed jobs, 2:for all jobs.      |
|                     | Default:1                                     |
+---------------------+-----------------------------------------------+

Last modified: Thu Jan 29 12:22:01 GMT 2015
