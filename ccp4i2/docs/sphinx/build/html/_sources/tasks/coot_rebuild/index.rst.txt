#######################
Manual model rebuilding
#######################

   |image1|

   The manual model rebuilding task will lauch the Coot graphical model
   building software, and automatically load files from the CCP4i2
   database. It is commonly launched as a follow-on task to an existing
   task.

   The task is also commonly used to view or examine results. If you use
   Coot and do not save any changes through the CCP4i2 extensions menu,
   the task will be removed from the CCP4i2 task list on completion. If
   you do save a modified model or dictionary, the task will remain in
   the task list as a record of those changes, and the saved objects
   will be recorded in the CCP4i2 database.

   *Note: for changes to be saved in the CCP4i2 database, you must save
   your changes through the CCP4i2 extensions menu*\ **(1)**\ *, not
   through the normal File/Save menu in Coot. If you use File/Save by
   mistake, you can re-import the saved model, but you will have no
   record of the task.*
   
Input
=====

   |image2|

   A selection of models, maps, difference maps, dictionary and saved
   state may be provided to Coot.

   One or more models can be selected **(1)**. Each of these will
   receive an entry in the Coot display manager.

   One or more maps can be selected **(2)**.

   One or more difference maps can also be selected **(3)**. These will
   be displayed with positive and negative contour levels.

   A ligand dictionary for any non-standard monomer types may also be
   provided **(4)**.

Results
=======

   .. container:: tab-content
      :name: tab-content2

      |image3|
      
      The report details all the data objects which were passed to Coot,
      and also any files which where saved from Coot into the CCP4i2
      database.

.. |image1| image:: coot_rebuild_task_0.png
.. |image2| image:: coot_rebuild_task_1.png
.. |image3| image:: coot_rebuild_report_1.png
