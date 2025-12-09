###########################################
Import or select subset from coordinate set
###########################################

   This tool enables you to select a limited set of atoms in a
   coordinate file and save them to a new file.

Input
=====

   .. container:: images

      |image1|

   After selecting a coordinate file for this task you can see its
   contents using the *View* option on the icon menu **(1)** (right
   mouse click on the icon) to show either a text listing of the file or
   a summary of the content.

   You then need to type into the *Atom selection* line the command to
   select the required chains, residues or atoms. The selection language
   is explained `here <../../general/atom_selection.html>`__; this page
   can also be accessed from the *Help* button **(2)**. Next to the text
   input the number of atoms currently selected is reported.

   There is also an option **(3)** to specify whether the coordiantes
   are a model, a homolog, a fragment of the structure or heavy atoms -
   this information will be used to ensure the new file is only used in
   the appropriate context.

   An alternative, more visual, approach to selecting atoms is with
   CCP4mg.

.. |image1| image:: ./coordinate_selector_task_1.png
