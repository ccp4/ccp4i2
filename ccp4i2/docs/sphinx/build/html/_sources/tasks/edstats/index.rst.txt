#####################################################
Analyse agreement between model and density (edstats)
#####################################################

Input
=====

   |image1|

   Fill **(1)** with 2mFo-DFc and mFo-DFc map coefficients and the model
   you want to analyse.

Options
-------

   |image2|

   It is possible to set different resolution cut-offs to the ones
   automatically detected from the map coefficient files by using
   **(3)**. You can set more strict or looser rejection criteria by
   using **(4)**. Those residues whose metrics lie outside these values
   will be posted directly to Coot through it's scripting interface.
   This can be achieved by selecting the 'Manual model building'
   follow-on task.

Results
=======

   |image3|
   
   This graph shows the different metrics, separated by main chain and
   side chain. Those residues whose metrics lie outside these values
   will be posted directly to Coot through it's scripting interface. Th
   is can be achieved by selecting the 'Manual model building' follow-on
   task.

.. |image1| image:: edstats_input.png
.. |image2| image:: edstats_options.png
.. |image3| image:: edstats_report1.png
