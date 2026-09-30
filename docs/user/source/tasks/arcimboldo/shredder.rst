===================
Arcimboldo Shredder
===================

`Arcimboldo Shredder <http://chango.ibmb.csic.es/shredder>`__ uses fragments derived from distant homologs based on the experimental data, instead of by sequence similarity. In this approach, the initial template is systematically shredded, and fragments are scored against each unique solution of the rotation function. Results are combined into a score per residue and the template is trimmed accordingly.

The Arcimboldo Shredder performs parallel and systematic molecular-replacement searches using a library of small structure motifs derived from a homologous structure and analyses the results to extract information from the persistence of solutions for different fragments among the noisy rotation-function results from Phaser.


Potential molecular-replacement solutions are passed to SHELXE for density modification and model building, with the prospect that any correctly placed fragments can be expanded into a full structure.

-----
Input
-----

+++++++++++++
Shredder mode
+++++++++++++

In the **sequential** mode, the initial template is systematically shredded and fragments are scored against each unique solution of the rotation function. Results are combined into a score per residue and the template is trimmed accordingly. Models are selected by optimization of PHASER’s rotation function and trimming in SHELXE. Sequential shredding can improve models where the high average deviation from the target is owing to dissimilar or flexible regions.

The **spherical** mode generates a model, starting from a distant homolog template that is annotated in terms of secondary and tertiary structure, a set of compact search models used as a library in ARCIMBOLDO BORGES, profiting from specific options in PHASER and SHELXE for model improvement. Spherical shredding can improve models where deviations of the template from the target are evenly distributed among a similar fold.

+++++++++++++++++++++++
Perform gyre refinement
+++++++++++++++++++++++

The maximum-likelihood GYRE method implemented in Phaser for optimizing the orientation and relative position of rigid-body fragments of a model after the orientation of the model has been identified, but before the model has been positioned in the unit cell. It performs the refinement of rigid-body groups against the rotation function target.

Atoms with different chain identifiers within an ensemble will be treated as independent rigid groups, refining their rotation and relative translation. At each orientation during gyre refinement, the amplitudes (but not the phases) of the structure factors of the symmetry-related copies of the molecular-replacement fragments oriented (but not positioned) in the unit cell can be calculated, giving a set of normalized structure-factor amplitudes for the rotating components.

In gyre refinement, an initial RMSD parameter is chosen as a tradeoff between convergence radius and sensitivity to coordinate accuracy, iterating refinement and decreasing the RMSD parameter estimation sequentially. The goal is to improve and select among the many possible models those with a true r.m.s.d. of below 0.6 Å , and thus susceptible to being expanded to the full solution in the density-modification and autrotracing step. The chain definition also changes between cycles in order to increase the number of fragments and thus the degrees of freedom for model refinement, as predefined in the template partitioning step.

After model translation, packing filtering, and refinement, models are either disassembled into predetermined rigid groups or refined (gimble refinement)

+++++++++++++++++++++++++
Perform gimble refinement
+++++++++++++++++++++++++

GIMBLE refinement is used for gimble the refinement of rigid-body fragments of the model after positioning. it is simply a re-implementation of Phaser's rigid-body refinement developed for ease of scripting. It performs the refinement of rigid-body groups against the translation function target.


++++++++++++++++++++++++++
Perform LLG-guided pruning
++++++++++++++++++++++++++

LLG-guided pruning is used to trim the model of residues that are not contributing signals to the LLG at the target rmsd value.

To find out more please read the article `here <https://dx.doi.org/10.1107%2FS2059798318004357>`__


    more Shredder tutorials you can find on `here <http://chango.ibmb.csic.es/tutorial_shredder>`__ and `here <http://chango.ibmb.csic.es/tutorial_shredder_spherical>`__



**References**

`McCoy, A.J., Grosse-Kunstleve, R.W., Adams, P.D., Winn, M.D., Storoni, L.C., Read R.J. (2007) Phaser Crystallographic Software. J. Appl. Cryst. 40: 658-674 <https://doi.org/10.1107/S0021889807021206>`_

`Sammito, M. D., Millán, C., Rodríguez, D. D., de Ilarduya, I. M., Meindl, K., De Marino, I., Petrillo, G., Buey, R. M., de Pereda, J. M., Zeth, K., Sheldrick, G. M. & Usón, I. (2013) Exploiting tertiary structure through local folds for crystallographic phasing. Nature Methods, 10, 1099–1101 <https://doi.org/10.1038/nmeth.2644>`_

`Sammito, M., Meindl,K., Ilarduya,I. M., Millán, C., Artola-Recolons, C., Hermoso, J. A., Usón, I. (2014) Structure solution with ARCIMBOLDO using fragments derived from distant homology models. FEBS J. 281, 4029–4045. <https://doi.org/10.1111/febs.12897>`_

`McCoy, A. J., Oeffner, R. D., Millán, C., Sammito, M., Usón, I., Read, R. J. (2018) Gyre and gimble: a maximum-likelihood replacement for Patterson correlation refinement. Acta Cryst. D74, 279-289, <https://doi.org/10.1107/S2059798318001353>`_

`Uson, I., Sheldrick, G.M. (2018) An introduction to experimental phasing of macromolecules illustrated by SHELX; new autotracing features. Acta Cryst. D74: 106-116 <https://doi.org/10.1107/S2059798317015121>`_

`Oeffner, R.D, Afonine, P.V, Millán, C, et al. (2018) On the application of the expected log-likelihood gain to decision making in molecular replacement. Acta Crystallogr D Struct Biol.;74(Pt 4): 245-255 <https://doi.org/10.1107/S2059798318004357>`_
