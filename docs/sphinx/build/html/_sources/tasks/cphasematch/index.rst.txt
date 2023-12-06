##############################
Phase comparison (cphasematch)
##############################

   The phase comparison task allows phase estimates from different
   sources to be compared.

Input
=====

   |image1|

   A set of reflection data and two sets of phase probabilities must be
   provided. The phases will be compared. The reflection data are used
   to calculate map correlations for maps calculated using the two phase
   sets.

Results
=======

      |image2|
      
      The mean phase error between the datasets is reported. A table of
      addition statistics, along with graphs against resolution, are
      provided.

      Weighted mean phase error includes a figure-of-merit to weight the
      phase comparison. This avoids the phase error being inflated by
      reflections for which no phase error is present, and also provides
      a suggestion as to the meaningfulness of the figure of merit.
      Weighting may be done using the figure of merit from either set of
      phase probabilities.

      Map correlations further introduce the reflection magnitude into
      the weighting calculation. This avoids the phase error being
      inflated by reflections for which there is little or no intensity.
      F map correlation tends to be dominated by phase errors at low
      resolution, whereas E-map correlation is dominated by phases at
      high resolution.

.. |image1| image:: cphasematch_task_1.png
.. |image2| image:: cphasematch_report_1.png
