###################
Estimate AU content
###################

   This task estimates the number of copies of a structure in a crystal
   cell based on the volume of the cell and the estimated structure
   molecular weight. It reports the number of copies of the structure
   that would give a crystal solvent content consistent with usually
   observed values. The estimated number of copies can be useful
   guidance to molecular replacement programs.
   `Details <../../../../html/matthews_coef.html#matthews_number>`__

Input
=====

   |image1|
   The task interface requires an MTZ experimental data file **1** from
   which it can read the cell parameters and calculate the cell volume.
   The expected AU contents can be entered as a *AU content* object
   which contains a list of seqeunces\ **2**, or the number of amino
   acid residues in the structure or the structure molecular weight.

   The results are calculated for the data provided and listed in a box
   at the bottom of the input tab **3**. If you 'Run' the task a
   permanent report of the result is created.

Results
=======

.. |image1| image:: matthews_task_1.png
