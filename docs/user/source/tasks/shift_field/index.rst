======================
Shift field refinement
======================

The task performs model morphing by shift field refinement. This is a very fast form of preliminary refinement which can be applied as a precursor to conventional refinement. Shift field refinement can be applied at any resolution (at low resolution is particularly suited to correcting large shifts in secondary structure elements in molecular replacement) and used to refine B factors.

-----
Input
-----

A reflection file and model coordinate file must be provided as a task's input. **Reflection file** should contain the data for the unknown, work structure. The required columns are F and phi or F, sigF, and a set of HL coefficients from a phasing, phase improvement, or refinement program.


For coordinate refinement, the calculation would normally start at low resolution (6Å), and then increase to 3 or 4Å before moving to regular refinement. 

.. image:: image_1.png
        :scale: 50 %
        :align: center

Depending on the speed and fit required, b-factor refinement may be usefully performed at higher resolutions. 

.. image:: image_2.png
        :scale: 50 %
        :align: center

1. Automatically set the radius in proportion to the resolution for the current cycle. The resolution is multiplied by this factor to get the radius (Default 4.0)

2. Select to perform pseudo-regularization

3. A number of cycles to perform. For coordinates, 10-20 is probably reasonable, B-factor refinement is faster (Default 1)

------
Output
------

The program outputs table of R-factors (and, if a flag is provided, free R-factor) at the start of each cycle. Note that these are based on the reflections used for a given cycle, which may change if the resolution changes from cycle to cycle. Values are therefore not comparable from cycle to cycle unless the resolution is unchanged. As a result, the values are only really useful for evaluating radius settings using the same resolution settings.

The free-R factor is very noisy at low resolutions and degrades performance, so the use of a free set is probably inadvisable.


----------------
Acknowledgements
----------------

This article uses materials kindly provided by **Professor Kevin Cowtan**, whose help is greatly appreciated.


Read more `here <http://www.ysbl.york.ac.uk/~cowtan/sheetbend/csheetbend.html>`_


`Macromolecular refinement using shift field optimization and regularization talk by K. Cowtan <https://www.youtube.com/watch?v=V9EmwP0mqUY&t=135s>`_


**References**

`Cowtan, K., Metcalfe, S. & Bond, P. (2020). Acta Cryst. D76, 1192-1200. <https://journals.iucr.org/d/issues/2020/12/00/di5041/>`_