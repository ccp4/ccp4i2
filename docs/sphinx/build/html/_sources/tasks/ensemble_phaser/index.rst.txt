============================
Build an ensemble for Phaser
============================

Phaser must be given the models it will use for molecular replacement. MR Ensembles are superposed sets of sufficiently similar MR models. Using MR Ensembles instead of individual MR models usually leads to a higher success rate of Molecular Replacement. 

An ensemble corresponds to a unit of structure that can be positioned in an MR run. It generally constitutes a "rigid body" that is likely to be common between the search model(s) provided and the structure for which the data have been collected. An ensemble should only be defined once, even if there are several copies of the molecule in the asymmetric unit. 

If you construct a model by homology modeling, remember that the RMS error you expect is essentially the error you expect from the template structure. Specify the sequence identity of the template, not of the homology model.


read more om `Phaser wiki <https://www.phaser.cimr.cam.ac.uk/index.php/Molecular_Replacement#How_to_Define_Models>`_

`Phaser Keywords <https://www.phaser.cimr.cam.ac.uk/index.php/Keywords#ENSEMBLE>`_
