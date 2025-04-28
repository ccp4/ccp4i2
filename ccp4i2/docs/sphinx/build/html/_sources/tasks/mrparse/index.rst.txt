########################################################################
MrParse: Finding search models in the PDB and the EBI Alphafold database
########################################################################

About
=====

MrParse is a program designed to aid in the decision making process in Molecular Replacement (MR).


Input
=====

MrParse can be found under the  menu in the CCP4i2 GUI:
|image1|

Opening MrParse will bring you to the following menu:
|image2|

To run MrParse all you require is a sequence file. If a reflection file is provided, MrParse will also calculate the eLLG for each of the hits in the PDB.
eLLG has been shown to be a better indicator of whether a search model will work in MR than sequence identity and so we recommended providing a reflection file if one is available.
Additionally you can select to perform sequence classification, to get the most out of this however, you may wish to install some additional software.


Additional Software
===================

MrParse attempts to classify the sequence according to its secondary structure and whether any regions are expected to be Coiled-coil or Transmembrane.

* Secondary structure classification is currently carried by submitting jobs to the [JPRED](http://www.compbio.dundee.ac.uk/jpred/) server.
* Coiled-Coil classification is carried out with [Deepcoil](https://github.com/labstructbioinf/DeepCoil). This needs to be installed locally.
* Transmembrane classification is carried out with [TMHMM](https://github.com/dansondergaard/tmhmm.py). This needs to be installed locally.

After installing Deepcoil and TMHMM, you can configure MrParse to find the relevant executable files by changing the executable paths in the MrParse configuration file ($CCP4/share/mrparse/data/mrparse.config).


Results
=======

When MrParse finishes, an HTML page will pop up showing the results of the search:
|image3|

The sections of the MrParse report page are highlighted in different colours:

* In red is information on the input reflection file, including resolution, space group and crystal pathology.

* In teal is information about the PDB entries identified by Phmmer and visualisations of the matches.

* In purple is the protein classification report, this includes a secondary structure prediction, a coiled-coil prediction, and a transmembrane prediction.

* In blue is information about the AlphaFold models identified by Phmmer and visualisations of the matches coloured by pLDDT on an orange to blue scale, where orange indicates very low confidence in the model and blue indicates very high confidence in the model.

Individual models can be downloaded by clicking on the name of the model in the report page. They can also be found by navigating to the homologs and models directories within the MrParse run directory.


.. |image1| image:: ccp4i2_mrparse.png
.. |image2| image:: ccp4i2_mrparse_run.png
.. |image3| image:: ccp4i2_mrparse_results.png
