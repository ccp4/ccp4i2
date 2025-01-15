=================
Arcimboldo Borges
=================

Arcimboldo `Borges <http://chango.ibmb.csic.es/BORGES/>`_, uses tertiary-structure searching in the PDB to solve the crystallographic phase problem. It combines the use of composite secondary-structure elements as a search model with density modification and tracing to reveal the rest of the structure when both steps are successful.
Before the use of BORGES it is important to perform a secondary structure prediction analysis in order to define an initial hypothesis as to local folds present in the structure.
BORGES may automatically perform an initial assessment with subsets of libraries to prioritize starting hypotheses according to Molecular Replacement figures of merit (FOM) or follow the order of the most frequent folds.

    **FOM** is a measure of phase quality used for evaluating initial experimental phases. FOM ranges is from 0 to 1, where 1 is the best possible score. 

-----
Input 
-----

++++++++++++++
Borges Library
++++++++++++++

Borges libraries previously created with ALEPH are distributed with CCP4 for use as input search models in Arcimboldo Borges.

The internal nomenclature **U** (up) and **D** (down) is used to describe the relative orientations of the fragments composing the fold thus, ``UUU`` means three parallel fragments and ``UDU`` means antiparallel.

A table of libraries is available `here <https://journals.iucr.org/d/issues/2020/03/00/ba5305/index.html#TABLE2>`_

++++++++++++++++++
Phaser GYRE option
++++++++++++++++++

The maximum-likelihood GYRE method implemented in Phaser for optimizing the orientation and relative position of rigid-body fragments of a model after the orientation of the model has been identified, but before the model has been positioned in the unit cell. It performs the refinement of rigid-body groups against the rotation function target.

Atoms with different chain identifiers within an ensemble will be treated as independent rigid groups, refining their rotation and relative  translation. At each orientation during gyre refinement, the amplitudes (but not the phases) of the structure factors of the symmetry-related copies of the molecular-replacement fragments oriented (but not positioned) in the unit cell can be calculated, giving a set of normalized structure-factor amplitudes for the rotating components.

In gyre refinement, an initial RMSD parameter is chosen as a tradeoff between convergence radius and sensitivity to coordinate accuracy, iterating refinement and decreasing the RMSD parameter estimation sequentially. The goal is to improve and select among the many possible models those with a true r.m.s.d. of below 0.6 A˚ , and thus susceptible of being expanded to the full solution in the density-modification and auto-tracing step. The chain definition also changes between cycles in order to increase the number of fragments and thus the degrees of freedom for model refinement, as predefined in the template partitioning step.

After model translation, packing filtering and refinement, models are either disassembled into predetermined rigid groups and refined (gimble refinement)

++++++++++++++++++++
Phaser GIMBLE option
++++++++++++++++++++

GIMBLE refinement is used for gimble the refinement of rigid-body fragments of the model after positioning. it is simply a re-implementation of Phaser's rigid-body refinement developed for ease of scripting. It performs the refinement of rigid-body groups against the translation function target.


------
Output
------

Arcimboldo Borges report components are:

* **Summary information** (spacegroup, cell dimensions, resolution, number of unique reflections)

* **Rotation clustering bar chart and table** illustrate the number of rotation clusters and their FOMs. 


**LLG** (Log-likelihood-gain)
 displays how well given model explains data in comparison with random atoms. The LLG has to be positive (as high as possible), and it should increase as the solution progresses. 
 
.. Note::
        The LLG allows to compare of different models against the same data set, but the LLG values for different data sets should not be compared with each other.

.. toggle-header::
        :header: **Guide to expected LLG values**

           Top solution correct?   

            <25      -no 

            25-36     -unlikely

            36-49     -possibly

            49-64     -probably

            >64       -definitely


**Zscore**
       The signal of placement is indicated by its translation-function Z-score (TFZ). It is computed by comparing the LLG values from the rotation or translation search with LLG values for a set of random rotations or translations. signal-to-noise is judged 

.. toggle-header::
        :header: **Guide to expected translation-function Z-score (TFZ-score)**

                    TFZ-score 

                    less than 5 	 **not** a solution

                    5 - 6	        **unlikely** a solution

                    6 - 7	       **possibly** a solution

                    7 - 8	        **probably** a solution

                    more than 8* 	**definitely** a solution

                        * 6 for 1st model in monoclinic space groups	


* A **search and expansion** table provides the characteristics and FOMs for each step and cluster

* In **backtracking** information user may find the links to the PDB of the best solution and its map

Once a solution is found (SHELXE traced mainchain CC > 25%), it stops after some recycling steps to improve it.


**CC** correlation coefficient 

    **CC more than 25% indicates the solution is found**

**CC** calculated on all data is reliable when atomic resolution data are available but at lower resolution, all random collections of a large enough number of unconstrained atoms show equally high CC values.


For Borges tutorial visit `www.chango.ibmb.csic.es/tutorial <http://chango.ibmb.csic.es/tutorial>`_


**References**

`Sammito, M. D., Millán, C., Rodríguez, D. D., de Ilarduya, I. M., Meindl, K., De Marino, I., Petrillo, G., Buey, R. M., de Pereda, J. M., Zeth, K., Sheldrick, G. M. & Usón, I. (2013) Exploiting tertiary structure through local folds for crystallographic phasing. Nature Methods, 10, 1099–1101 <https://doi.org/10.1038/nmeth.2644>`_

`McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D., Storoni, L.C., Read R.J. (2007) Phaser Crystallographic Software. J. Appl. Cryst. 40: 658-674 <doi:10.1107/S0021889807021206>`_

`McCoy, A. J., Oeffner, R. D., Millán, C., Sammito, M., Usón, I., Read, R. J. (2018) Gyre and gimble: a maximum-likelihood replacement for Patterson correlation refinement. Acta Cryst. D74, 279-289, <https://doi.org/10.1107/S2059798318001353>`_

`Uson, I., Sheldrick, G.M. (2018) An introduction to experimental phasing of macromolecules illustrated by SHELX; new autotracing features. Acta Cryst. D74: 106-116 <doi:10.1107/S2059798317015121>`_

`Millán C, Sammito M. D, McCoy A. J, Nascimento A. F. Z, Petrillo G, Oeffner R. D, Domínguez-Gil T, Hermoso J. A, Read R. J, Usón I (2018) Exploiting distant homologues for phasing through the generation of compact fragments, local fold refinement and partial solution combination. Acta Cryst. D74: 290-304 <doi:10.1107/S2059798318001365>`_


