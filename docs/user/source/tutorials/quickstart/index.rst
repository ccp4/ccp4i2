##########################################
CCP4i2 — Quickstart Introduction to CCP4i2
##########################################

|CCP4 web logo|

Creating a Project
------------------

.. container::
   :name: quickstart_step

   |image1|

.. container::
   :name: quickstart_annotation

   | You should organise your work into projects - usually one project
     per structure solution but ther are no rules.
   | A new project is created by selecting **New project** option on the
     **Projects** pull-down menu.
   | Each project has an associated directory for the files created in
     the project and you will have to create or choose an existing
     directory and select a name for the project. A *Project viewer*
     should then open. You can open the project by selecting it from the
     **Projects** menu at any time.

--------------

The project window - create a job
---------------------------------

   |image2|

   | If you need some test data go to *Download test data* on the
     *Utilities* menu and choose the *rnase* data.
   | Double click on a task in the *Task menu* to open a job - try
     *Refinement* first.

--------------

Input data for a job
--------------------

   |image3|

   | 1 The new job is shown in the *Job list*
   | 2 The *Input Data* tab is open for you to select data - other tabs
     have other parameter input
   | 3 Missing or invalid essential data is highlighted in red
   | 4 Click on browser button to open file browser - note that in the
     file browser you can download data from EBI or other databases.

--------------

Importing a file
----------------

   |image4|

   | Whenever you select a new file CCP4i2 will make a copy to the
     project directory so nothing is lost and ..
   | 1 this window asks you for information on the source of the file to
     save to the CCP4i2 database. The progam provides a sensible default
     which you may edit.
   | 2 or you can click to make the bothersome window go away!.
     Unchecking 3 will stop these windows appearing when you load files
     in future.

--------------

Importing experimental data - the tricky bit
--------------------------------------------

   |image5|

   The four types of experimental data: reflections, phases, freeR sets
   and map coefficients are put in separate files in CCP4i2 (`more
   discussion <../general/data_files.html>`__) but if you have an
   old-style 'monster' MTZ file you will need to select the right
   columns from that file so it can be imported into a new-style 'mini'
   MTZ. By default CCP4i2 will import the richest source of data
   (e.g.anomalous intensities have more information than a structure
   factor and sigma so they are preferred). Don't worry - these will be
   converted automatically for any program that can not handle a
   particular form of data. **You must also load the same file into the
   Free R set input for this task. CCP4i2 does not do this
   automatically!**

--------------

Running the job
---------------

   |image6|

   | 1 After selecting the required data click the *Run* button.
   | 2 A running man icon on the *Job list* and in job status bar
     indicates the job is running.
   | 3 For many tasks the *Results* tab opens to display a constantly
     updating report.

--------------

The job has finished
--------------------

   |image7|

   | 1 The *Job list* shows the job, any sub-jobs, imported files and
     output files - right-mouse click on these for pop-up menu with
     options to view job log files, database info or view files in Coot
     or CCP4mg.
   | 2 After looking at the report choose your next job - click on the
     *Task menu* button to get back to the list of tasks or
   | 3 under the report are one-click buttons to open probable next
     jobs. Click on *COOT* to create a Coot task.

--------------

*Manual model rebuilding* input data
------------------------------------

   |image8|

   | When you open a new job CCP4i2 tries to fill-in the input data by
     taking the output files (or sometimes the input files or even
     earlier files) from the *follow-on* job.
   | 1 The job used to provide the input files is indicated at the top
     of the job *Input data* - there is a menu to select the job.
   | 2 The icon by any data can be right-mouse clicked for a menu with
     options to *View* a file and other useful tools. It is also
     possible to drag a data object from the job list to this data input
     icon.
   | 3 The data may have been selected automatically but you do need to
     click the *Run* button to start *Coot*.

--------------

Saving models from Coot
-----------------------

   |image9|

   | 1 You should do your model building in *Coot* as usual but save the
     model by the *Save to CCP4i2* option on the *CCP4i2 extensions*
     menu or the *File* menu. It is not necessary to do the usual *Save
     Coordinates* in *Coot*.
   | 2 The *Manual model rebuilding* job and the model data are saved in
     the CCP4i2 database and displayed in the *Job list*. Note: if you
     have not saved any data to *CCP4i2* when you close *Coot* the job
     will be removed from the *Job list*.

See the `Overview <../roadmap/index.html>`__ for a summary of the data
and functionality in *CCP4i2*.

.. |CCP4 web logo| image:: ../../general/ccp4i2.png
   :width: 75px

.. |image1| image:: ./new_project.png

.. |image2| image:: ./proj_win.jpg

.. |image3| image:: ./input_data.png

.. |image4| image:: ./provenance.png

.. |image5| image:: ./mtz_import.png

.. |image6| image:: ./running.png

.. |image7| image:: ./job_done.png

.. |image8| image:: ./coot_data_input.png

.. |image9| image:: ./coot_save.png
