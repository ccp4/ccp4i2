===============
Arcimboldo Lite
===============

`ARCIMBOLDO_LITE  <https://doi.org/10.1107/S1399004715010846>`_ runs in two macrocycles. In the first macrocycle, all copies of the search fragments are placed with Phaser and optimized with `SHELXE <../shelx/shelx.html>`_. Similar solutions are clustered. A number of substructure solutions equal to the number of physical cores minus one is subject to density modification and autotracing with SHELXE, prioritizing according to top LLG and CC. Arcimboldo Lite runs on all cores of a single machine, without requiring a grid, supercomputer or database.

-----
Input
-----

*  The default value of the rmsd for model helices is 0.2 Å (this value can be increased).

*  It is mandatory to define the number of helices to search and their length (length needs to be specified according to secondary-structure prediction or previous knowledge about the target structure).

*  The molecular weight, number of components and estimated rmsd should be defined for the PHASER search.

++++++++++++++++
Coil coiled mode
++++++++++++++++

    *represent a search algorithm to probe and verify alternative helix directions*

Coil coiled mode entails:

* VRMS calculation in the refinement step to optimize the RMSD parameter in order to maximize the LLG.

* Activation of PHASER’s packing filter at translation, so at least one translated solution will pass the packing check.

* Generate and probe reverse helices. At resolution below 2 Å it was frequently observed that placement of the first helices occasionally took place in the correct position but in reversed direction and at low resolution are not distinguishable from the PHASER figures of merit.

* Final verification step, an additional step that generates perturbations of the substructure leading to the best solution and compares their scores before and after extension. As the reverse helices this will be only activated at a resolution worse than 2 Å.

* SHELXE with helical sliding, which improves the autotracing of the coiled-coil structures.


        adapted from http://chango.ibmb.csic.es/tutorial_coiled.

------
Output
------

Arcimboldo report window summarizes your job results. Search and expansion tables are updated in real-time and can help to judge parameterization.

A search and expansion tables show:

**LLG** (Log-likelihood-gain)
 displays how well the given model explains data in comparison with random atoms. The LLG has to be positive (as high as possible), and it should increase as the solution progresses. 
 
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

 

**CC** correlation coefficient 

    **CC more than 25%  indicates the solution is found**

**CC** calculated on all data is reliable when atomic resolution data are available but at a lower resolution, all random collections of a large enough number of unconstrained atoms show equally high CC values.




For online documentation and tutorials visit `www.chango.ibmb.csic.es <http://chango.ibmb.csic.es/>`_

ARCIMBOLDO_LITE tutorials can be found `here <http://chango.ibmb.csic.es/tutorial_arc>`_


**References**

`McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D., Storoni, L.C., Read R.J. (2007) Phaser Crystallographic Software. J. Appl. Cryst. 40: 658-674 <https://doi.org/10.1107/S0021889807021206>`_

`Rodríguez DD, Grosse C, Himmel S, González C, Martínez de Ilarduya I, Becker S, Sheldrick GM & Usón I (2009) Crystallographic ab initio protein structure solution below atomic resolution. Nat Methods 6, 651– 653. <https://www.nature.com/articles/nmeth.1365>`_

`Sammito, M., Millan, C.,Frieske, D., Rodriguez-Freire, E., Borges, R. J., & Uson, I.(2015) ARCIMBOLDO-LITE: single-workstation implementation and use  Acta Cryst. D71, 1921-1939. <https://doi.org/10.1107/S1399004715010846>`_

`Millán C, Sammito M, Usón I (2015) Macromolecular ab initio phasing enforcing secondary and tertiary structure. IUCrJ 2, 95-105. <https://doi.org/10.1107/S2052252514024117>`_

`Millán C, Sammito M. D, McCoy A. J, Nascimento A. F. Z, Petrillo G, Oeffner R. D, Domínguez-Gil T, Hermoso J. A, Read R. J, Usón I (2018) Exploiting distant homologues for phasing through the generation of compact fragments, local fold refinement and partial solution combination. Acta Cryst. D74: 290-304 <https://doi.org/10.1107/S2059798318001365>`_

`Caballero I, Sammito M, Millán C, Lebedev A, Soler N, Usón I (2018) ARCIMBOLDO on coiled coils. Acta Cryst D74, 194-204, <https://doi.org/10.1107/S2059798317017582>`_

`Uson, I., Sheldrick, G.M. (2018) An introduction to experimental phasing of macromolecules illustrated by SHELX; new autotracing features. Acta Cryst. D74: 106-116 <https://doi.org/10.1107/S2059798317015121>`_
