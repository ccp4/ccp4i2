##################
Define AU contents
##################

   This task allows you to import the sequences expected in the
   structure from a variety of sequence file formats or a coordinate
   file. You can enter the expected number of copies - it is not
   critical that this is correct. The information is saved to a *AU
   content* object that is used in any subsequent task requiring
   sequence information. The *AU content* object corresponds to a file
   with the extension *.asu.xml*. You can use this task to edit an
   existing *AU content* - select the object to load it and then make
   any changes.

Input
=====

   |image1|
   You can select an existing *AU content* object to edit\ **(1)**.

   You may need to click on a *Show list* button to enter more that one
   sequence **(2)** and then use the +/- buttons **(3)** to edit the
   list of sequences.

   A sequence may be entered by cut-n-paste into the window or loading
   from a sequence or coordinate file **(4)**. Any white space and most
   unexpected characters in the input will be removed before saving but
   if the sequence contains invalid letters (BJOUXZ) this will be
   highlighted as they may be from comments or some other issue that
   requires your attention. If a loaded file contains multiple sequences
   you will have the chance to select those required and the sequence
   list will be extended for more than one sequence input.

   If you already have a sequence loaded using earlier versions of
   CCP4i2 then you can add the seqeunce object to a AU content.

Results
=======

      |image2|
      The imported sequences and input and output files are reported.

.. |image1| image:: ProvideAsuContents_task_1.png
.. |image2| image:: ProvideAsuContents_report_1.png
