########
splitMtz
########

Input

.. container:: tab-content
   :name: tab-content1

   CCP4i2 now uses `experimental data
   objects <../../general/data_files.html>`__ where each object is one
   'mini' MTZ file. Each experimental data object is one instance of
   either reflections, phases, map coefficients or free R sets. This
   task will split an old style, 'monster' MTZ, containing mixed data,
   into 'mini' MTZs.

   .. rubric:: Input
      :name: input

   |image1|
   Select the old style, 'monster' MTZ file. **(1)**

   The column data in the file is analysed and the new style
   experimental data objects are listed. **(2)** Beware some redundant
   data in the old-style files does not correspond to any experimental
   data type (e.g. *ISYM* columns). Check the boxes to select the
   experimental data objects you want to import to CCP4i2.

.. |image1| image:: splitMtz_task_1.png
