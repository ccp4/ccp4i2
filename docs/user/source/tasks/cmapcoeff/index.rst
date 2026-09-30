##############################################
Calculate unusual map coefficients (cmapcoeff)
##############################################

Refinement produces the usual 2mFo-DFc and mFo-DFc map coefficients; this
task makes the others. An *anomalous difference* map, from one dataset of
anomalous pairs and a set of phases, has peaks at the anomalous
scatterers: it identifies metals, halides, sulphurs or, as here, xenons,
and confirms their positions. A *weighted difference* map, from two
datasets and one set of phases, shows what changed between them, a bound
ligand or a heavy atom. A *weighted F* map uses one dataset.

The pictures on this page calculate the anomalous difference map of the
xenon-derivative data of the *gamma* demo that comes with CCP4i2.

Input
=====

.. figure:: cmapcoeff_input.png
   :alt: Figure 1: cmapcoeff input

   Figure 1: cmapcoeff input

Choose the kind of map **(1)**. For an anomalous difference map the first
dataset **(2)** must be anomalous pairs (intensities or amplitudes), with
phases; the second **(3)** is used only for a difference between two
datasets.

The coefficients can be sharpened (a negative B-factor) or blurred (a
positive one), or cut at a different resolution **(4)**, and a map file
can be written as well **(5)**.

Results
=======

The report lists the output map coefficients. Open them in Coot or
Moorhen, or search them for peaks. Here the two strongest peaks, at 41 and
18 times the map's RMS, are xenon sites (they match the demo's
heavy_atoms.pdb, which is on an origin shifted by ½ along a); the next
peaks are below 6 RMS and at no site. Peaks that stand this far above the
rest are anomalous scatterers; a peak of 4-6 RMS may be a weaker one, a
sulphur, say, or noise.
