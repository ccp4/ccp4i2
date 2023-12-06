##########################################################
Non-programming requirements for implementing CCP4i2 tasks
##########################################################

Supporting 'What next?'
~~~~~~~~~~~~~~~~~~~~~~~

After each job the GUI presents the user with a list of reasonable tasks
to perform next and a web page with more discussion of the options. It
is the responsibility of the task author to provide the information for
these facilities. Note that if the user opts to follow from a sub-job in
a pipeline then it is the parent pipeline that is queried to provide the
'what next' information.

There are two mechanisms to provide the information; the simplest is as
a WHATNEXT attribute of either the wrapper class (CPluginScript) or gui
class (CTaskWidget). So, for example, for buccaneer the WHATNEXT is set
as a list of names of tasks..

::

   class CTaskBuccaneer(CCP4TaskWidget.CTaskWidget):
     TASKNAME = 'buccaneer'
     ...
     WHATNEXT = ['refmac_simple']

The more flexible approach is to have a function whatNext() in either
the wrapper module or the gui module (note it is a function and not a
member of the class). This is passed four arguments:

+----------------+----------------------------------------------------+
| jobId          | The jobId of the job (or the jobId of the parent   |
|                | pipeline job if a sub-job is selected)             |
+----------------+----------------------------------------------------+
| taskName       | The taskname of the job (or of the sub-job if a    |
|                | sub-job was selected)                              |
+----------------+----------------------------------------------------+
| childJobNumber | The job number of the sub-job if a sub-job was     |
|                | selected                                           |
+----------------+----------------------------------------------------+
| projectName    | The project name                                   |
+----------------+----------------------------------------------------+

This function should return either a list of tasks that can be run next
or a list of task name and text to appear on the GUI. For example for
buccaneer_build_refine:

The WHATNEXT parameter or the return from *whatNext()* is normally
expected to be a list of tasknames but the list item could be a list
containing two items: the taskname and the filepath of a params file
that provides default parameters to the next task. The params file path
should be of the form

::

   $CCP4I2/pipelines/bucref_mr/scripts/bucref_after_experimental.params.xml

i.e. the file is expected to be found somewhere in the CCP4I2 directory.
The *bucref_after_experimental.params.xml* is used to set one input
parameter in the *Autobuild protein* gui so that it expects input from
experimental phasing.

In principle this mechanism could be used to pick up a project or
task-specific file if the *whatNext()* function returns a list including
a reference to a file in the current project directory of the form

::

   $PROJECT/CCP4_JOBS/job_5/for_whatever_next_task.params.xml

This is interpreted in the context of the project viewer that
understands $PROJECT to mean the current project directory. (Note this
is not tested.) It is the responsibility of either the wrapper or the
function to provide the params file in the project directory.

Each task should provide a web page with discussion of what to do next.
This should be $CCP4I2/docs/whatnext/*taskName*.html.

Bibliographic References
~~~~~~~~~~~~~~~~~~~~~~~~

| Task reports will normally include a section of bibliographic
  references for software used by the task. This is presented as text
  with options to link to the original paper (if available on-line) and
  to download the reference in either MedLine (suitable for EndNote) or
  BibTeX format. The task implementer must provide two files containing
  the references in the appropriate formats. These files should be
  called *wrappers/mytask/script/mytask.medline.txt* and
  *wrappers/mytask/script/mytask.bibtex.txt*. The first can be created
  and downloaded from the PubMed site:
| http://www.ncbi.nlm.nih.gov/pubmed/advanced
| After selecting the references use the *Send to* tool at the top left
  of the page to create a file in MEDLINE format. To get the BibTeX file
  input the PUbMed Ids (PMIDs) (or other suitable search query) here:
| http://www.bioinformatics.org/texmed/

The MedLine file is used by CCP4i2 to generate the information in the
report and if a copy of the paper is available online the MedLine entry
should be edited to append the a line such as:

::

   URL - http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3069738/

If there are multiple references then the first in the file will be
treated as the one for which citation is preferred.

Note that the 'MEDLINE' format file can be imported into *EndNote* as
*PubMed(NLM)* format (found in the *Other filters* menu of the *Import*
tool.

| These references are accessed from three places in the gui:
| References for the current task can be accessed from the *Help* menu
  in the project viewer. These references are autogenerated from the
  contents of the *medline* file for the task and recursively for
  sub-tasks that are specified by the script's SUBTASKS parameter.
| References are also automatically generated and appended to a job
  report. But if the tasks used by a pipeline are only determined at
  run-time then it may be worthwhile specifying references in the report
  explicitly see `here <./reports.html#references>`__.
| PLACEHOLDER - not yet implemented: autogenerated for all tasks used in
  determining a data file seelcted by the user.

Task documentation
~~~~~~~~~~~~~~~~~~

Each task should have a page of documentation named
$CCP4I2/docs/tasks/*taskName*/index.html and any images or other
material for the page should be in the same directory. There is a
template html file: *ccp4i2/docs/tasks/resources/template.html*. See the
*parrot* documentation for example of required style.

In order to support installation of third party created pipelines with
all resources in one directory the task documentation can as be in
$CCP4I2/pipelines/*taskName*/docs/index.html.

Creating documentation images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There is a customised screen grab tool in the *Developer tools* sub-menu
of the *Utilities* menu. The tool will grab either the task input frame
or report as it appears in the Project viewer. You can resize the window
and use the scroll bar to set up the required frame size. Beware with
the task input you should resize the window to eliminate any space below
the input widgets and for reports try to keep all screenshots the same
width for a neater presentation in the web page.

| For the first screenshot of a new task you will be prompted to select
  a directory to contain the screenshots and these will be automatically
  named:
| *directory/taskname\_*\ task\ *\_1.png* or
| *directory/taskname\_*\ report\ *\_1.png* etc.

To standardise documentation screenshot appearance please use the images
of numbers in *ccp4i2/docs/resources*. One way to add these to the
screenshot is using *Inkscape*: open the screenshot image then use
*File->Import* to import each label *svg* file and then drag the label
to correct place. Save by using the *File->Export PNG image*.

Last modified: Fri Oct 2 15:01:56 BST 2015
