==========================================
MR with model splitting using Slice-n-Dice
==========================================

Slice-n-Dice is an automated MR pipeline designed to pre-process AlphaFold2 and RosettaFold models by removing low confidence regions and converting confidence scores into predicted B-factors. It also slices predicted models into distinct structural units then automatically places the slices using Phaser. The slicing step can use AlphaFold2’s predicted aligned error (PAE), or can operate via a variety of Cɑ atom clustering algorithms, extending applicability to structures of any origin. 

Slice-n-Dice helps to deal with inaccuracy in domain-domain orientations of predicted models. 

**Structure prediction software puts pLLDT scores (AlphaFold) or rmsd estimates (RoseTTAFold) in place of B-factors**. This gives a wrong perception regarding the regions with low or high reliability, and, therefore, the **scores must be recalculated as B-factors**.  That's why important to **select B-factor treatment option** correctly. The default **pLDDT threshold** for AlphaFold2 models is **70** and the default **RMS threshold** for RosettaFold models is **1.75**.


-----
Slice
-----

Slice option in CCP4i2 uses `Birch clustering algorithm <https://scikit-learn.org/stable/modules/clustering.html#birch>`_ based on the coordinates of the Cɑ atoms.


The **number of splits** (clusters) can be specified in the Options tab. MIN and MAX splits are the minimum and the maximum number of splits to make to the search models. For example, a MIN split of 1 and a MAX split of 3 will result in 3 MR jobs:

#. using an unsplit model, 

#. a model split in two and 

#. a model split three ways. 

.. note:: The large target protein needs to be modeled in several chunks. An increasing number of splits allows slice and dice to test a range of different splits.


----
Dice
----

Dice method performs molecular replacement on the individual slices produced by Slice by `Phaser <../phaser_pipeline/index.html>`_. Phaser automatically places as many slices as possible. A number of cores for Phaser can be specified in the Options tab.
After that `Refmac <../prosmart_refmac/index.html>`_ is used to assess if the placed slice has improved the solution. A number of refinement cycles in Refmac can be specified in the Options tab.

**ACKNOWLEDGEMENTS**

This article uses materials kindly provided by Dr. Adam Simpkin and Dr. Ronan Keegan, whose help is greatly appreciated.

**References**

`Simpkin, A. J., Elliott, L. G., Stevenson, K., Krissinel, E., Rigden, D. J., Keegan, R. M. (2022) Slice’N’Dice: Maximising the value of predicted models for structural biologists, bioRxiv 2022.06.30.497974; <https://doi.org/10.1101/2022.06.30.497974>`_

`Murshudov, G.N., Skubak, P., Lebedev, A.A., Pannu, N.S., Steiner, R.A., Nicholls, R.A., Winn, M.D., Long, F., and Vagin, A.A. (2011)  Acta Cryst. D67: 355-367; <https://doi.org/10.1107/S0907444911001314>`_

`McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D., Storoni, L.C., Read R.J. (2007) Phaser Crystallographic Software. J. Appl. Cryst. 40: 658-674; <https://doi.org/10.1107/S0021889807021206>`_

