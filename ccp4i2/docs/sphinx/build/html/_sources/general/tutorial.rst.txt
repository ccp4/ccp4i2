#######################
CCP4i2 — Brief tutorial
#######################

`Brief tutorial <tutorial.html>`__

`Experimental data <data_files.html>`__

`Model data <model_data.html>`__

`Atom selection <atom_selection.html>`__

Creating and Managing Projects
------------------------------

The CCP4 interface requires you to group your work into projects. There
are no rules on what goes into one project but, as a guideline, one
project would usually result in one structure. Corresponding to each
project is a directory that contains all the files associated with
project. By default, the directory for a new project is created in the
master *CCP4I2_PROJECTS* directory (in your home area) and given the
same name as the project but you have the option to create or use a
different directory. Please do not move the project directory or any of
the files, except using tools in the *Project Manager* window, otherwise
CCP4i2 will then be unable to find them! The earlier CCP4i had a
database file for each project kept in the project directory but the new
CCP4i2, by default, has one database for all of the projects of one
user. CCP4i2 maintains the database (and other important data) in a
profile directory usually called .CCP4I2 in your home area.

It is possible to order projects in an hierarchy - this is new in CCP4i2
- do this by drag-and-dropping the projects in *Project Manager* window.
Note that if you do change the logical heirachy of projects in CCP4i2
the project directories are not moved. CCP4i2 does not allow you to
create project directories within other project directories as this
complicates managing the project.

Starting CCP4i2 for the first time you should see a web browser page
with options for `creating a project <./welcome.html>`__. Follow these
to create either a *General* project or a specific project. Follow the
instructions and you should then see a project window for that project.

Project Windows
---------------

On the left side of the window is a list of jobs and on the right is a
view of one selected job. To open a new job click the **Task menu**
button from the toolbar and select a task from the list that appears on
the right-hand side (double click a line in the list).

The Task
--------

The task frame on the right hand side of the project window has three
tabs:

-  **Input**: Fuller options for task input. Most tasks will have
   multiple tabs including **Input data**. If you prefer the 'old style'
   folder layout to tabs then this can be set in the **Preferences**.
-  **Results**: A task report appears when the job completes. If there
   is no report available then a log file or error diagnostic may
   appear.
-  **Comments**:A place for you to add any comments on the job.

Entering data for a job
-----------------------

Most tasks only require you to select the appropriate input data and
they can then run with reasonable default parameters. The quick way to
select input data is selecting the previous job to 'follow from'. The
task input will be set as the output from that job, or if there is no
appropriate output then the input or the output from preceding jobs. To
select a file you can:

-  Select from the drop down list of appropriate data from preceding
   jobs
-  Click **Browse files** to use a file browser
-  Drag a file icon from the job list or a job report onto the file icon
-  Drag a job icon from the job list onto the file icon

If input data is missing or invalid this is indicated by a red border to
the data widget. Click the **Run** or **Save** button to run the task.
There is a warning message if data is missing or invalid.

The **Follow from** facility
----------------------------

You will frequently want to take the data output by one job directly
into the next task. Typically this data will be experimental data or a
model. To make this easier the interface to each task has a **Follow
from** option which allows you to choose which preceding job to take
data from. If the *follow from* job has not output any data of the
appropriate type then CCP4i2 will search back through the input to the
*follow from* job and preceding jobs to find data of the appropriate
type. The selected data is then displayed in the interface and you can
change the individual items of data independently if necessary.

The Job list
------------

The job list shows a job number and name with an icon indicating its
current status (pending, running, failed, finished). A job be opened to
show sub-jobs and output files. A right mouse click on the job list will
give a context menu with options to Delete, Clone or View a job or View,
Edit label or Export a file. If you wish to re-run a job then you should
clone it - all the data and parameters that you have set will be copied
into the new task.

The next column in the job list shows the 'follow from' job as green
up-arrows. By default the data output by this job is input to the next
opened job. The 'follow from' status can be changed by right mouse
clicking the green arrow or dots and selecting the new status from the
menu. The subsequent columns in the table are the job rating and summary
that you may have entered in the **Status** frame.

The Project directory
---------------------

This is an alternative view to the job list. It is a standard directory
and file listing of the project directory where there is one directory
(called job\_\ *n*) for each job. This is most useful for
trouble-shooting - ideally users should not need to do this!

The Profile directory
---------------------

CCP4i2 maintains a directory in your home area which contains your
personal data and preferences. Particularly it contains the database for
all of your projects and so should be treated with care and backed up!
On Linux and Mac systems the profile directory is .CCP4I2
(*dot*\ CCP4I2) in your home area and on Windows it is CCP4I2 (no dot)
in ??.

On Linux and Macs it is possible to set an alternative profile directory
in the environment variable CCP4_LOCAL_HOME. This may be most useful in
a workshop situation where multiple users have the same account. The
users will have to start CCP4i2 from the command line with something
like:

::

   prompt> CCP4_LOCAL_HOME="/whereever/master_dir/my_dir" ccp4i2

where each user has their own *my_dir* in a *master_dir*.

Running jobs and processes
--------------------------

A CCP4i2 job is started immediately when you click the *Run* button on
the toolbar. A 'Running man' icon and gold colouring on the *Job list*
indicate that the job is running. For long-running tasks intermediate
'running' reports are displayed and you may also see *sub-tasks* and
intermediate data files listed in the *Job list*. The job is run as a
separate, background process that is also called 'ccp4i2' and this may
start one or more processes running CCP4 programs. When the background
process finishes it updates the CCP4i2 database with the job status
(usually 'Finished' or 'Failed'). The graphical CCP4i2 window regularly
checks the job status in the database and will update the 'Job list'
when the job finishes.

The *Utilities* menu item *List running jobs and processes* will list
all of the jobs that are running according to the CCP4i2 database and
for each job will show the 'ccp4i2' background process and any program
processes that are 'children' of the ccp4i2 process. If no 'ccp4i2'
process is listed for a job then it appears that the job has failed
without updating it's database status and it is recommended that you use
the *Mark as failed* button to correct the database.

This tool also lists processes running on another
`server <./servers.html>`__ machine via the *Run remote* button.

.. |CCP4 web logo| image:: ccp4i2.png
   :width: 75px
