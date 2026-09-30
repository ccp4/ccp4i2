########################
Autobuild RNA (nautilus)
########################

   The "Autobuild RNA" task is used to either build an atomic model into
   an electron density map, or to rebuild a molecular replacement model
   to better fit the observations. The minimum required inputs are a the
   protein sequence, a set of observations, a set of Free-R flags, and
   either a set of phases or a molecular replacement model.

Input
=====

   |image1|

   First the experimental data must be selected. A set of observations,
   a set of phases and a set of Free-R flags are required\ **(1)**.

   Next the contents of the asymmetric unit must be specified in the
   form of one or more chain sequences. Each sequence need only be
   provided once, however many times it appears in the asymmetric
   unit\ **(2)**.

   A partially built model may optionally be provided as a starting
   point for building\ **(3)**.

   Non-nucleotide features from the initial input model are preserved
   automatically. Nucleotide features will be preserved unless a more
   compelling trace is obtained for the same density.

Options
-------

   |image2|

   The first group of options controls the model building pipeline. The
   only one of these you might normally change is the first: the number
   of build/refine iterations to be performed. If an initial model
   building job gives an incomplete model, rerunning with more cycles
   may help. While there is very rarely a benefit in using more than 20
   cycles, we have seen rare examples where hundreds of cycles lead to a
   poor map being solved.

   The final option - anisotropy correction - applies a correction to
   the data to remove the effect of the data being sharper in one
   direction than another. This should normally be on.

Advanced options
----------------

   The "Resolution" option allows the use of higher resolution data, at
   a cost of run time.

Results
=======

   |image3|
   
   The report summary describes how well the calculation has worked.
   Nautilus builds fragments of nucleotide chain, and tries to match
   those with the expected sequence. If the model is good, it will
   usually consist of a few long fragments which are grouped into as
   many chains as are present in the asymmetric unit. The chains will be
   largely complete, and all nucleotides built will be part of one of
   the chains. The refinement statistics are also an important
   indication of model accuracy.

   Better models typically consist of more nucleotide units built, and
   few fragments (unless there is a lot of NCS). The results table
   provides a comprehensive analysis of the state of the final model.

   Follow-on tasks include manual model building - use this if the model
   is nearly complete or if building has stalled.

.. |image1| image:: nautilus_build_refine_1.png
.. |image2| image:: nautilus_build_refine_2.png
.. |image3| image:: nautilus_build_refine_report_1.png
