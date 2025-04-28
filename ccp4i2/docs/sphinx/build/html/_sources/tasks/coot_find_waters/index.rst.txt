#####################
Find waters with COOT
#####################

   The 'Find Waters' task searches an electron density map for any
   water-like peaks which are not already accounted for by atoms in an
   atomic model. Water molecules are added to account for any peaks
   which appear to be waters on the basis of size and distance from
   existing features.

   The map is masked by the masking atoms and a search is made of
   features in the map about the electron density cut-off value. Waters
   are added if the feature is approximately water-sized and can make
   sensible hydrogen bonds to the protein atoms. The new waters are
   added to the atomic model as a new chain.

Input
=====

   |image1|

   The input data objects are an atomic model describing the known parts
   of the structure, and a set of electron density map coefficients. The
   map should be a normal rather than difference map, i.e. weighted
   F\ :sub:`obs` coefficients from experimental phasing or density
   modification, or 2mF\ :sub:`o`-DF\ :sub:`c` coefficients from
   refinment.

Options
-------

   |image2|
   
   The minimum and maximum allowable distances between new waters and
   the masking molecule (usually the protein) may be set. Defaults are
   2.4 and 3.2 Å.

Results
=======

      The results tab reports the number of waters added to the
      structure, along with the input and output files.

.. |image1| image:: coot_find_waters_task_1.png
.. |image2| image:: coot_find_waters_task_2.png
