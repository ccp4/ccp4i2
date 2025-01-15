#####################
Generate a Free R set
#####################

   A set of Free-R reflections may be selected at random for a given set
   of observations. The proportion of the reflections assigned to the
   Free-R set defaults to 5%, but can be altered. The Free-R set are
   selected in the highest symmetry point group consistent with the
   current point group, in case the spacegroup has been wrongly
   determined. The Free-R set is described by a data object containing
   flags which identify which reflections belong to the Free-R set.

Input
=====

   |image1|
   To generate a Free-R set for a given set of observations, all that is
   required is the reflection data object for those observations.

   Sometimes it is necessary to extend an existing set of Free-R flags
   to cover additional reflections, for example if a new set of
   observations is available at a higher resolution. In this case the
   same Free-R set should be retained for the existing reflections to
   avoid biasing the Free-R factor. In this case the 'complete an
   existing set of FreeR flags' checkbox should be checked, and the
   existing Free-R set object selected.

Options
-------

   |image2|
   There is one option which allows the proportion of data in the free-R
   set to be selected. By default this is set to 5% of the data (a value
   of 0.05), however by checking the box a different proportion may be
   selected. For very small structures, a larger value may be required
   to provide sufficient reflections in the Free-R set in any given
   resolution shell.

.. |image1| image:: freerflag_task_1.png
.. |image2| image:: freerflag_task_2.png
