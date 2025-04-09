======
DIMPLE
======

`DIMPLE <https://journals.iucr.org/a/issues/2013/a1/00/a50958/a50958.pdf>`_ (DIfference Map PipeLinE)Automated difference map generation pipeline to allow for quick assessment of crystal data to see if a ligand has bound to the structure. It is a joint development by the CCP4 software group and the Diamond light source. It is an automated software pipeline for rapidly processing crystals that contain a known protein and possibly a ligand bound to this protein. 

**DIMPLE** takes the already known “apo” structure for the target protein and compares it with the electron density map for the same, crystallised structure (and possibly bound ligand) interpreted from the X-ray diffraction images. It then presents the user with a set of snapshot images illustrating the regions where there are large areas of electron density unaccounted for by the protein structure model. Within a few seconds, an experienced user can decipher from these images if any of this unaccounted for density indicates a bound ligand. A positive result can mean that the rest of the batch of crystals for this particular ligand can be discarded allowing for a more efficient use of the beam time and saving the user from lots of potentially unnecessary processing.

The main goal is to present a user with a quick answer to the question of whether or not they have a bound ligand or drug candidate in their crystal. 

**DIMPLE** can work only with reflection datasets that include mean intensities.
 
Please visit the `official site <http://ccp4.github.io/dimple/>`_ to get more information.

    
**References**

`Wojdyr, M., Keegan, R., Winter, G., Ashton, A. (2013) DIMPLE - a pipeline for the rapid generation of difference maps from protein crystals with putatively bound ligands. Acta Cryst. A69: s299 <https://doi.org/10.1107/S0108767313097419>`_

