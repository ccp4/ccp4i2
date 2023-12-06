######################################
Merge experimental data objects to MTZ
######################################

   CCP4i2 now uses `experimental data
   objects <../../general/data_files.html>`__ where each object is one
   'mini' MTZ file. Each experimental data object is one instance of
   either reflections, phases, map coefficients or free R sets. This
   task will recreate an old-style 'monster' MTZ from experimental data
   objects. Note that there may be an easier way to get the data you
   require: on the CCP4i2 *Job list* for many jobs the context menu
   (right-mouse click) has an *Export* sub-menu which will include an
   option to export an old-style MTZ of the data associated with that
   job.

Input
=====

   |image1|
   The input is a list of experimental data objects with the column
   names to be given to the columns in the output file. Select the
   experimental data object **(1)** and edit or remove the tag to be
   added to each column name.\ **(2)** The default output column names
   are shown as a comma-separated list of names compounded from the
   column names in the input file and the tag. These can be edited
   manually if necessary.\ **(3)**

   The list of data objects can be extended by clicking the 'Add
   input..' button.\ **(4)**

   Enter a name for the output file.\ **(5)**

.. |image1| image:: ./mergeMtz_task_1.png
