######################
CCP4i2 — Customisation
######################

`Brief tutorial <tutorial.html>`__

`Experimental data <data_files.html>`__

`Model data <model_data.html>`__

`Atom selection <atom_selection.html>`__

`Customisation <customisation.html>`__

Introduction
~~~~~~~~~~~~

| There are several tools to enable users to customise CCP4i2:
| `Workflows <#workflows>`__ allow you to save and rerun a series of
  tasks
| `Task parameters and patches <#com_file_patchs>`__ are used to save
  and restore task control parameters and/or to edit command files in
  order to access program features not available in the interface
| `Custom tasks <#custom_tasks>`__ allow you to run programs not
  supported by the interface - you must provide a command line and
  command files

These customisations can be set up and managed from the *Customisation*
pull-down menu. Each customisation has a manager window with options to
create, edit and delete. Each customisation can be exported (as a
compressed file) or imported so it is possible to share with colleagues.

The various customisations are saved in `your profile directory/custom
in the appropriate directory: workflows, comfilepatchs, tasks or
importedjobs. Each customisation has its own directory which is given
the customisation name that you provide and the contents of this
directory will depend on the type of
customisation. <./tutorial.html#profile>`__

Workflows
~~~~~~~~~

Workflows are a mechanism for users to save a series of tasks that they
use frequently and rerun them as one task. The new workflow will have an
auto-generated user interface to select the minimal input data. The
user's preferred control parameters from the template tasks are saved
and reused in each run of the workflow.

| To create a workflow you should run the series of tasks that you want
  to encapsulate as a workflow. Be particularly careful about what
  output data from one task is chosen as input for another task because
  the workflow manager will automatically copy the same data flow in the
  generated workflow. When the tasks have run and completed successfully
  choose *Workflows* from the pull-down menu and select *New* in the
  Workflow Manager window. In the next window you will have to select
  the project and the template jobs and provide a unique name for the
  workflow. The workflow manager will then analyse the template jobs and
  saves the workflow in a directory: *your profile
  directory/custom/workflows/workflow name*. You are then presented with
  window with options to customise the interface to the workflow. The
  same window can be recalled by clicking the *Edit* button in the
  Workflow Manager window. The options are:
| *Short title for workflow task* - enter the text that will appear in
  the task menu and the project viewer. *Labels for input files* - This
  is a list of the input data files required by the workflow - check
  that it matches with your expectation. By default the input gui will
  have labels appropriate for the file type but you can change these to
  give specific hints if necessary.
| *Select overall workflow output files* - This is a list of files
  output by all of the tasks in the workflow and by checking the boxes
  you can control which files are displayed in the project viewer as the
  overall output of the workflow. By default the output files of the
  last task are the overall output files. For each file you can also
  enter a short description that will appear in the project viewer.

The workflows that are present in the directory *your profile
directory/custom/workflows* are listed at the top of the task menu and
clicking on one of this will bring up a task interface for the required
input files. There is a second tab on the interface, *Control parameters
for sub-tasks* and under this is a list of the sub-tasks which can be
clicked to open a window with an interface to edit the control
parameters of the task.

| Workflows work well where there is input and output of just one
  coordinate file but some known limitations..
| If there is output of more than one cordinate file (or any other type
  of file) then only the first file will be used in subsequent tasks.
| Tasks which require addidional input (eg Phaser requiring Identity/RMS
  data for each input ensemble) do not work in a workflow.

Task Parameters and Patches
~~~~~~~~~~~~~~~~~~~~~~~~~~~

| This is a tool to do one or both of two things:
| Save the control parameters used in one job so they can be restored
  and used in another job with one click
| Edit a program command file and save the changes to apply to the
  command file for subsequent runs of the program. This could be useful
  to access program functionality not currently interfaced by CCP4i2 or
  possibly to set multiple parameters with one-click rather than editing
  multiple options in the interface.

To use this tool you should first run the program in question from
CCP4i2 in order to create a template command file that can be edited
later. If the job is long-running and not useful you can kill it without
loosing the command file. From the *Customisation* pull-down choose
*Task parameters and patches* and click the *New* option. In the *Task
Parameter and Patches Manager* window you should first choose the
project and job that will provide the parameters and the template
command file. You should check/uncheck the option to save the task
control parameters and then edit the template command file is necessary.
This command file patch tool uses an approach similar to the *Unix diff*
and *patch* tools that are described at
`Wikipedia <http://en.wikipedia.org/wiki/Diff>`__; these will probably
work best if you make the changes as simple as possible - ideally as
additional lines at the end of the file. You should enter a name and
title (to appear on the interface task menu) for the patch before saving
it.

When there are task parameters and patches available for a task then
will be listed in the interface for that task at the bottom of the
*Input data* tab. Selected task parameters will be applied immediately
to the parameters in the interface.Selected patch will be applied to the
command file when the task is run. You might want to check that the
resultant command file is correct by viewing it in the *Project
directory* on the left-hand side of the Project Viewer window; the
command file is in the appropriate job directory and normally called
*com.txt*.

Custom Tasks
~~~~~~~~~~~~

This is a mechanism for you to run non-CCP4 programs from CCP4i2 with a
simple interface, similar to the standard CCP4i2 task interfaces,
enabling users to select input data and set important parameters. The
custom task jobs with their input and output data are recorded in the
database and shown in the project viewer the same as any other job. To
set up a custom task you must provide the text of a command line to run
a program and also, if it is used, the text of a command file. The text
can include variables (such as input file names) as parameter names. You
must then provide more information on each variable parameter - the data
type (such as file type) is essential. Based on the information
provided, a simple interface is auto-generated for the user to input the
variables such as input file names.

| The *Custom task* interface to set up a custom task is accessed from
  the *Customisation* menu. The interface to enter a new custom task has
  two main sections - one to enter the command line and, optionally,
  command file text and below that a section to define the variable
  parameters that appear in the command line or file. You should first
  enter the command line with variable names beginning with a hash tag -
  something like:
| *myprogram #HKLIN #PDBIN -ncycles #NCYCLES*
| Here the program *myprogram* has command line parameters: HKLIN - an
  MTZ experimental data file, PDBIN - structure coordinates and NCYCLES
  - the number of cycles. The parameter name begins with a hash tag
  though and alternative tag can be specified at te top of the
  interface. The command file text can be entered in similar style. By
  default when the custom task is run a command file will be written to
  the task's working directory and will be called *com.txt*. It is also
  possible to specify the name of the command file if the program
  expects files with a specific name.

After entering the command line/file text click the *Update parameters*
button to fill out the parameters section with the parameters from the
command line/file section. For each parameter you will need to choose:

-  a label that will appear in the interface - this is optional and
   defaults to the usual label for the data type
-  the data type (i.e. probably the type of file)
-  function (input,output..).
-  if the parameter is obligatory (check box)
-  if the parameter is recorded in the database (check box)

It is important to set the 'Record in database' option correctly.
Usually coordinate models and experimental data such as new phases
should be recorded in the database so the flow of data can be followed
and output data can be passed on to subsequent programs. Other data such
as log files or program-specific 'library' files or program control
parameters should not be recorded in the database but they are still
saved in parameter files.

For parameters that are files there is an optional field to enter a file
name. This name is relative to the work directory where the job is run.
This filename only needs to be entered if the program requires an input
file with a specific name (in which case the user-selected file input
file will be copied to that filename) or produces an output file whose
filename is not specified in the command line or command file (in which
case the file will be copied from the filename to a CCP4i2 standard
filename). For the output files the filename can be a filename wildcard
(a.k.a. 'pathname expansion' or 'glob') which can include a '*' to match
any character string. There are further wildcard options described
`here <#pathname_expansions>`__.

The most difficult issue may be handling
`mini-MTZs <./data_files.html#MTZ>`__: CCP4i2 handles experimental data
as separate observed data, phases, FreeR flags or map coefficients but
presently most programs expect input with these combined to a single
MTZ. The way to handle this for custom tasks is to have the single input
MTZ parameter in the command text and this has the data type
*CMtzDataFile*. You then create additional parameters in the parameter
list with one parameter for each mini-MTZ that constitutes the
monster-MTZ. (To add a parameter to the parameter table click the '+'
symbol by the parameter table.) The extra parameter should be assigned a
mini-MTZ datatype (CObsDataFile,CPhsDataFile, CMapCoeffsDataFile or
CFreeRDataFile) and then an extra line appears in the interface with an
option to *Merge to or split from monster MTZ*; this has a combo box to
select from any defined *CMtzDataFile* file. The mini-MTZs that you
define will appear in the user interface for the task and the system
will automatically merge them to create the monster-MTZ that is the
program input file. To handle splitting an output monster-MTZ you also
need to enter the names of the MTZ columns that are split out into the
mini-MTZs. There is another possible complication if the program will
only accept observed data or phases in a limited set of representations
(see explanation `here <./miniMTZs.html#representations>`__). If you set
a data type *CObsDataFile* or *CPhsDataFile* with *input* function then
a list of the possible representations will appear and you can un-check
any which are unacceptable to the program. The custom task run mechanism
will attempt to convert data to an appropriate representation
automatically.

Distributed with CCP4i2 there is an example custom task for running the
model building program `ARP-wARP <http://www.embl-hamburg.de/ARP>`__.
The ARP-wARP program itself is not automatically included in the CCP4
distribution and you may need to download it from
`CCP4 <http://www.ccp4.ac.uk/download>`__ or the
`ARP-wARP <http://www.embl-hamburg.de/ARP>`__ site. To import the
example custom task choose *Custom task* from the *Customisation* menu
and choose *Import* in the *Custom task Manager* dialog box and select
the file *Your CCp4i2
directory/test/custom_tasks/arp_trace.ccp4_customtask.tar.gz*. You can
then double-click the newly-imported *arp_trace* task to edit and see
the contents. At the top of the *Create custom task* window you can see
the name and the title (text that will appear in the interface).

The ARP-wARP autotrace script is run with a simple command line that
includes the name of a monster-MTZ file that must include observed data
and phases. In the *Create custom task* window the *Command line* field
is the path to the ARP-wARP autotrace script (which you will probably
need to edit appropriately for your computer) and keyword, 'datafile'
and then the parameter '#HKLIN'. There is no command file for running
ARP-wARP so that field is empty. The list of parameters is important for
ensuring the correct behaviour of the interface. HKLIN is the one actual
input file to ARP-wARP and is a monster-MTZ file of type *CMtzDataFile*
but it is made up from OBSIN (observed data with data type CObsDataFile)
and PHSIN (phases with data type CPhsDataFile) which have the *Merge to
MTZ* value set to 'HKLIN' which is the cue to the system to create the
monster-MTZ HKLIN from those files. For the input observed data ARP-wARP
will only accept a limited set of representations and so the
*Representation* field is set to 'Fmean'.

The output files from ARP-wARP are PDBOUT which has type CPdbDataFile
and HKLOUT, a monster-MTZ of type CMtzDataFile. For both of these the
*File path* field is filled in with a path to the expected output file.
ARP-wARP creates a directory for output files with a default name
devised from the current data and time so we can not predict it but we
can use a '*' wildcard in the file path and, though it may look odd in
this case, the output filenames include the name of the input MTZ file
which will be 'HKLIN'. There are two further mini-MTZ output files
defined, MAPOUT and DIFFMAPOUT that will contain map coefficients so
have type 'CMapCoeffsDataFile' - these are split out from the program
output HKLOUT file so have 'HKLOUT' in the *Merge to MTZ* field and have
entries in the *List of column labels* field. Finally note the *Record
in database* column: the two input and output monster-MTZ files are not
recorded as they are regarded as temporary intermediate files but the
input and output mini-MTZs and the output PDBOUT file are recorded.

Once the example ARP-wARP custom task is loaded the task title will
appear in the *Task menu* in the *Custom tasks* section which is near
the top of the list of tasks. If you click on this then the interface
appears - the input requires just two files for the observed data and
phases. The ARP-wARP script takes some time to run but you can follow
progress by viewing the log file (right-mouse click on the job in the
left project view window to get a context menu with a *View* option.)
You can also see what files are created by the job by choosing the
*Project directory* view from the tab options on the left side of the
Project Viewer window and opening the directory folder for the
appropriate job.

Appendix - pathname expansions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some examples of pathname expansions:

-  **./*.pdb** finds files in the current directory with any filename
   and the file extension '.pdb'
-  **./*/*.pdb**\ finds files in any sub-directory with any filename and
   the file extension '.pdb'
-  **./[0-9]???.pdb** finds files in the current directory with a
   filename beginning with a character in the range 0 to 9 (i.e any
   number), followed by three characters that may be anything and the
   extension '.pdb'

Last modified: Fri Feb 20 16:17:24 GMT 2015

.. |CCP4 web logo| image:: ccp4i2.png
   :width: 75px
