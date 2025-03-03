=======================
SAD phasing with Phaser
=======================


.. image:: phaser_ep_task_1.png
      :scale: 70 %


When solving a structure using anomalous data, it can be useful to look for (additional) anomalous scatters which can be used to improve the phasing calculation, or in order to locate anomalous scattering atoms such as sulfurs which may help in building the structure. Anomalous scatterers may be located by using the anomalous differences in combination which phases from a partial structure, for example from molecular replacement.

Input
=====

In order to generate an optimal anomalous map, the reflection data (which must include anomalous pairs of reflections) and a set of atomic coordinates (to provide phase information) are required.

Before carrying out SAD phasing in Phaser, you need an atomic model. Several model options are available in CCP4i2.


**Find heavy atoms**

The minimal Phaser input requires specification of the atom type that provides the anomalous scattering (e.g. Se, S, I, Fe, Zn etc), the number of atoms can be specified if known.

Alternatively, you can provide a partial model of the structure of the macromolecule, with or without a set of anomalous scatterers.

**Partial protein/nucleic acid model**

For this mode, partial model must be provided. Phaser will phase SAD data and completes the structure from log-likelihood gradient maps using a partial structure as the initial atomic substructure. Only a non-anomalous signal from the partial structure is used. 

It's possible to provide an electron density map as partial model by choosing **Map coeeficents for protein\nucleic acid** option

A set of **resolution limits** may optionally be specified to reduce the noise due to missing low resolution terms or larger phase errors at high resolution.

The **composition of the unit cell** must also be specified, to allow for the contribution of unaccounted for scattering matter to be evaluated. This may be specified in terms of molecular weights or sequences. First, select the type of information to be provided, and then provide the required information in the boxes below.

**Destiny modification** is performed by `Parrot <../parrot/index.html>`_. It improves the initial phases using prior information about protein electron density maps.

Structure rebuilding is usually required after phasing, SAD phasing with Phaser task uses the `Buccaneer <../buccaneer_build_refine_mr/index.html>`_ as a **model building** program.

Read more about Phaser on `Phaser wiki <https://www.phaser.cimr.cam.ac.uk/index.php/Experimental_Phasing>`_

Keywords
========

   |image2|
   For documentation of the advanced options, see the `Phaser
   documentation <http://www.phaser.cimr.cam.ac.uk/index.php/Keywords>`__.

.. |image2| image:: phaser_EP_LLG_task_2.png


**References**

`A.J. McCoy, R.W. Grosse-Kunstleve, P.D. Adams, M.D. Winn, L.C. Storoni & R.J. Read. "Phaser Crystallographic Software". J. Appl. Cryst. (2007). 40, 658-674 <https://doi.org/10.1107/S0021889807021206>`_
