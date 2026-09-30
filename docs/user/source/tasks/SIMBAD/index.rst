###############################
Detecting contaminants (SIMBAD)
###############################

About
=====

SIMBAD (Sequence-Independent Molecular replacement Based on Available Databases) is used for both contaminant testing and brute-force approaches for molecular replacement. It is designed to check for the presence of a contaminant in a crystal data set. The pipeline is also useful in cases of the misidentification of a crystal as well scenarios where no obvious homologue is available as a search model or the most suitable search model is not among those most highly ranked by sequence comparisons.

In CCP4i2 SIMBAD contains two steps:
* lattice-parameter search against the entire Protein Data Bank, rapidly determining whether or not a homologue exists in the same crystal form.
* contaminant search, a not uncommon occurrence in macromolecular crystallography. To catch for this eventuality, SIMBAD rapidly screens the data against a database of known contaminant structures. This database is compiled to include entries from ContaBase and Dimple.


Input
=====

SIMBAD can be found under the Molecular Replacement menu in the CCP4i2 GUI:
|image1|

Opening SIMBAD will bring you to the following menu:
|image2|

To run SIMBAD all you require is an MTZ file. SIMBAD will run a lattice parameter search followed by search against a database of common contaminants.
The user can also set the contaminant search to only search for contaminants from a specific organism (e.g. the expression organism).

Additionally the following advanced options can be selected:
|image3|


Lattice Results
===============

When SIMBAD finishes, the user is presented with a summary tab. This provides an overview of the best solution found by SIMBAD:
|image4|

Clicking on "Lattice Parameter search results" tab will display two tables, the results from the lattice parameter search and the molecular replacement results:
|image5|


Contaminant Results
===================

When SIMBAD finishes, the user is presented with a summary tab. This provides an overview of the best solution found by SIMBAD:
|image6|

Clicking on "Contaminant search results" tab will display two tables, the results from the contaminant rotation search and the molecular replacement results:
|image7|



.. |image1| image:: ccp4i2_simbad.png
.. |image2| image:: ccp4i2_simbad_run.png
.. |image3| image:: ccp4i2_simbad_run_2.png
.. |image4| image:: ccp4i2_lattice_summary.png
.. |image5| image:: ccp4i2_lattice_results.png
.. |image6| image:: ccp4i2_contaminant_summary.png
.. |image7| image:: ccp4i2_contaminant_results.png
