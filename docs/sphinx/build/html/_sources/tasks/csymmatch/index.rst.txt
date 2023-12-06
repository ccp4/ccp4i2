#################################
SYMMETRY MATCH MODEL TO REFERENCE
#################################

   The model matching task is a tool to enable the simple comparison of
   atomic models. It uses the crystal symmetry to move chains, ligands
   and waters of an atomic model, in order to place them in the same
   part of the unit cell as an existing model. It can optionally be used
   to change the origin on which the model is defined, for example for
   comparing different molecular replacement results.

Input
=====

   |image1|

   Two atomic models must be selected. The first atomic model **(1)** is
   the model to be moved. This model will be broken up into connected
   fragments, and each fragment will be moved using the crystal symmetry
   to best overlap the second model. A distance criterion is used to
   determine what constitutes a connected fragment - by varying this
   criterion the behaviour of the program can be adjusted.

   The second atomic model **(2)** is the model to be matched. This
   model is fixed, and is not moved or altered in any way.

   Two options may be specified **(3)**: an optional change of origin
   and/or hand may be applied to best match the models. A change of hand
   is meaningless for most biological structures, but is useful when
   comparing heavy atom substructures which are typically subject to
   hand ambiguity. Note that the origin shift is not restricted obey
   symmetry constraints. 
   
   The connectivity criterion (radius to use for stitching) may also be
   specified: the model is normally moved a chain at a time, however if
   monomers within the same chain are separated by more than the
   connectivity distance (default 2A) they are moved separately. This allows chains
   to remain intact while moving solvent molecules and ligands
   individually.

Results
=======

      |image2|
      
      The report describes any origin or hand shift which has been
      applied to the model, and the number of connected fragments into
      which the model was divided for matching. Check for a change of
      hand - a hand change is meaningless for anything other than a
      heavy atom substructure, and usually indicates a problem in
      matching the molecules. Also check whether any change of origin is
      consistent with the spacegroup symmetry.

      The model fragments and associated symmetry transformations are
      also reported.

.. |image1| image:: csymmatch_task_1.png
.. |image2| image:: csymmatch_report_1.png
