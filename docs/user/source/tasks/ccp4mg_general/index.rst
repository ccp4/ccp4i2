#############################################################
Molecular graphics vizualization and figure creation - CCP4MG
#############################################################

Input
=====

   CCP4mg is the quick and easy way to create beautiful publication
   quality images and movies. It is also the easy way to superpose and
   analyse structures.

   The range of model drawing styles include ribbons, surfaces, spheres
   and ball-and-stick. Surfaces can be coloured by electrostatic
   potential calculated by the built-in Poisson-Boltzmann method.
   Nucleic acids can be drawn with ribbons, base blocks or base pair
   sticks. CCP4mg also displays electron density as cylinders or 2D
   slices, text with fine control over colour and font and any vectors
   in a range of styles. You can set up complex pictures quickly using
   the Picture Wizard and fine tune images with the display table and
   customised drawing styles.

   Full details of the program can be found at `the CCP4MG
   website <http://www.ccp4.ac.uk/MG/>`__.

   |image1|
   
   In the simplest case, CCP4MG is started with a structure to load
   specified in the task input. CCP4MG may then be used to select a
   subset of atoms from the loaded structure and then save these to the
   CCP4i2 database. Once the CCP4MG session is finished CCP4i2 will
   present a gallery of the structures which have been saved. The atom
   selection may be done using any of CCP4MG's selection tools:- simple
   menu selections, Selection Browser, Sequence Viewer. The user has
   full control of which atoms are in the selections to be saved. `The
   CCP4MG documentation" <http://www.ccp4.ac.uk/MG/>`__ should be
   consulted for further information.

   More than one structure may be loaded at a time by adding entries to
   the *Atomic model* list (Show list, followed by "+"). Also one or
   more electron density or difference density mtz files may be loaded,
   and one or more sequences.

   Once the CCP4MG session is started, structures, maps, sequences, etc.
   may be loaded into CCP4MG by its usual means. Any saved selections or
   ensembles will become output files for the task and appear in the
   program report. Details of saving selections and ensembles may be
   found in the `"CCP4MG Edit Search
   model" <../ccp4mg_edit_model/index.html>`__ task documentation.

Results
=======

      The results tab presents a gallery of the structures saved to the
      CCP4i2 database during the CCP4MG session. These may be viewed in
      Coot or CCP4MG again.

.. |image1| image:: ccp4mg_general_task_1.png
