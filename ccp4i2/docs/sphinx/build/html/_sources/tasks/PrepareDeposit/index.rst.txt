############################
Prepare and validate files for deposition
############################

   The final stage of the crystallographic structure solution is the
   deposition of the final structure and the observations in the Protein
   Data Bank (PDB). The PDB require that the files be provided in a
   standard form, specifically the observations and the atomic model
   must be in CIF format. This task will automatically prepare the
   required files in the appropriate formats.

Input
=====

   |image1|

   Firstly, the contents of the asymmetric unit must be supplied **(1)**.
   This must contain the complete sequences of all macromolecule chains
   present in the coordinate file i.e. including parts present in the
   crystallized species but not modelled.

   Next add the coordinates from the final refinement job **(2)**.
   At this point all the information from this job will be harvested
   from the database and the related files section will be populated.

   You may choose not to send the data and coordinates to the wwPDB
   Validation Service **(3)**. In this case the task will prepare the
   files but you will not receive a validation report.

   The option to add statistics from data reduction **(4)** is
   described below.

   If you have run the data reduction task in this CCP4i2 project this
   task should be able to identify the correct job that produced the
   data that were used in refinement and add data reduction statistics
   to the coordinate mmCIF file. 
   
   There are two reasons why this automatic detection does not work.
   Firstly, data reduction was not done within CCP4i2. However, if data
   reduction was done within CCP4i2 but the space group used differs
   from that of the refined structure the detection will fail. In both
   cases this pop-up will appear.

   |image2|

   How to proceed depends on whether data reduction was done within
   CCP4i2. 

   |image3|

   If data reduction was done within CCP4i2 you can locate the output
   XML from Aimless that contains the data reduction statistics **(5)**.
   Alternatively you can choose to not add data reduction statistics to
   the coordinate mmCIF file **(4)** and instead enter this information
   into the wwPDB OneDep system. 

Results
=======
   If you chose to send the coordinates and data to the wwPDB
   Validation Service the task will take at least five minutes (but
   often much longer) to run as it waits until  the Validation Service
   returns the validation report. 
   
   |image4|

   If you chose to add data reduction statistics a summary of these is
   provided and/or chose to send the coordinates and data to the wwPDB
   Validation Service a summary is shown **(1)**.

   |image5|

   You can view and export the coordinates **(2)**, reflection data
   **(3)**, and validation report **(4)**. You can also run a Coot job
   **(5)** to inspect issues highlighted in the validation report.


.. |image1| image:: PrepareDeposit_task_1.png
.. |image2| image:: PrepareDeposit_task_2.png
.. |image3| image:: PrepareDeposit_task_3.png
.. |image4| image:: PrepareDeposit_report_1.png
.. |image5| image:: PrepareDeposit_report_2.png
