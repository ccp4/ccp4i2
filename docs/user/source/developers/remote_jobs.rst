##############################
CCP4I2 Developers: Remote Jobs
##############################

Introduction
~~~~~~~~~~~~

Presently (May 2013) CCP4i2 has a mechanism to create a compressed file
containing all data files, command files and database info needed to run
a single job. The *runTask* process can be given the name of the
compressed file as an argument and it will unpack the file and run the
job and then repack the resultant data files and database status to
another compressed file which can be loaded into the original project
directory and database. Currently there is no mechanism to move the
compressed files or start a job on a remote machine. I would expect that
to be handled by a separate non-graphical process that supports all CCP4
software and provides access to various services. What does exist in
CCP4i2 is primarilly a demonstration of exporting and importing data for
a single job (or a job with sub-jobs) that is related to *Project
Management* tools for import and export of a project or part of a
project. To test the 'run remote' facilities the *IFTESTREMOTE* flag in
*qtgui/CCP4ProjectViewer.py* must be set *True* so the *Run* button on
the project gui is converted to a menu with options to 'Remotely - test
output' and 'Remotely - import result'. Selecting the first of these
will create a 'setup' compressed file for remote running and will report
the name of the file and give a command line that can be used to run the
job. You should use the command line to run the job and create another
'final' compressed file. You can then use the 'Remotely - import result'
to import from the final compressed file.

The remote-job compressed file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The files used to move data for remote-running jobs are a special case
of the files used to move projects. Currently they are *tar.gz* files
supported by the Python *tarfile* module. The 'setup' compressed file
contains a fragment of the project directory and maintains exactly the
same directory structure as the source project directory - this includes
the *CCP4_JOBS* directory and may include the *CCP4_IMPORTED_FILE* and
the *CCP4_PROJECT_FILES* directory. The *CCP4_JOBS* directory will
contain the directory for the running job and may contain other job
directories which contain input files. The 'final' compressed file
contains only the job directory for the run job (and any sub-job
directories) since the job should only output files to its own
directory.

The compressed file also contains an XML representation of a fragment of
the database with just the data pertaining to running the single job. In
order to run the job independently the *runTask* process creates a new,
temporary sqlite database and loads it with the data from the XML file.
This database is updated by the running process and, when the job is
complete, the entire updated database is written out to XML and packaged
with the 'final' compressed file.

Last modified: Thu May 16 13:07:03 BST 2013
