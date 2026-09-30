#####################
Generate a Free R set
#####################

A set of Free-R reflections may be selected at random for a given set of
observations. The proportion of the reflections assigned to the Free-R
set defaults to 5%, but can be altered. The Free-R set is selected in the
highest symmetry point group consistent with the current point group, in
case the space group has been wrongly determined. The Free-R set is
described by a data object containing flags which identify which
reflections belong to it.

Data reduction (*Aimless*) and *Import merged data* make a Free-R set
already, so this task is needed only for data that arrived without one,
or to extend an existing set to new data.

Input
=====

.. figure:: freerflag_input.png
   :alt: Figure 1: Generate a Free R set input

   Figure 1: Generate a Free R set input

To generate a Free-R set for a set of observations, all that is required
is the reflection data object for those observations **(1)**.

Sometimes it is necessary to extend an existing set of Free-R flags to
cover additional reflections, for example if a new set of observations is
available at a higher resolution. In that case the same Free-R set should
be kept for the existing reflections, to avoid biasing the Free-R factor:
choose *complete an existing set of freeR flags* **(2)** and select the
existing Free-R set.

The proportion of the data in the Free-R set **(3)** is 5% by default (a
value of 0.05). For very small structures a larger value may be needed, to
give enough Free-R reflections in each resolution shell. The set can also
be limited to a high resolution **(4)**.

Results
=======

.. figure:: freerflag_report.png
   :alt: Figure 2: Generate a Free R set report

   Figure 2: Generate a Free R set report

The report says what was done, and the new Free-R set is listed among the
outputs with its space group, resolution and cell, which identify the
data it belongs with.
