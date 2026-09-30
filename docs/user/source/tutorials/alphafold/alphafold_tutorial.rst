CCP4i2 MR Tutorial: Preparing AlphaFold/RoseTTAFold search models for MR. (Gamma example)
========================================================================================

This tutorial explains how to handle AlphaFold or RoseTTAFold predicted models in CCP4i2 for use as search models in 
Molecular Replacement. In addition, we will demonstrate how to search for homolgous predicted models in the EBI-AFDB database using the CCP4mg-MrBUMP graphical interface. We will use the Gamma tutorial example for the demonstration. 

We assume that the user is familiar with the basic concepts of file importing and running other tasks in CCP4i2 such as MR with Phaser or refinement with Refmac.

#. **Generate predicted models using AlphaFold and RoseTTAFold**

   Use one of the AlphaFold Colab servers
   (e.g. `AlphaFold.ipynb <https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb>`_) 
   to generate a prediction for our target by giving it the sequence. 
   Additionally, you can create a RoseTTAFold model by providing the sequence to the server (`Robetta server <https://robetta.bakerlab.org>`_).
   The sequence for our Gamma example is:
      
      MHHHHHHLVPRGSHMIPSITAYSKNGLKIEFTFERSNTNPSVTVITIQASNSTELDM
      TDFVFQAAVPKTFQLQLLSPSSSVVPAFNTGTITQVIKVLNPQKQQLRMRIKLTYNH
      KGSAMQDLAEVNNFPPQSWQ 

   It may take some time for the servers to generate their predictions (typically an hour or more). 
   You can proceed with steps 2 to 5 and then 7 and subsequent steps, returning to step 6 when 
   the servers have produced predictions. 

   (**Note**: RoseTTAFold predictions will contain non-standard lines in the header of the pdb file. 
   Use an editor to remove these lines  before importing the coordinates into CCP4i2. 
   Remove all lines above the "MODEL    1" line.)

#. **Create project and import tutorial data**

   From the Files/Projects menu create a new project called "AF-Gamma" or something similar.
   Go to the Utilities menu and choose the "Copy demo data" option and select the "Gamma" data.
   This will copy the Gamma files to the project folder.

#. **Import Sequence(s) Task**

   Import the gamma.pir sequence from the project folder using the "Import sequence(s)" task.

#. **Import merged reflection data Task**

   Import the "merged_intensities_Xe.mtz" reflection file from the project folder using the "Import merged reflection data" task.

#. **Define AU Contents Task**

   Load the imported gamma sequence and the merged data. The Task should
   indicate the most likely AU content is one molecule in the asymmetric unit.

#. **Import a coordinate set Task**

   (*Skip to step 7 if still awaiting server results*) 
   
   When you have downloaded the predicted models from AlphaFold and/or RoseTTAFold, 
   import them into CCP4i2 using the "Import coordinate set" task.
   
   **Important!** The model predictions will contain pLDDT scores (AlphaFold) or rmsd 
   estimates (RoseTTAFold) in the B-factor column. Before we use the models in MR, we 
   need to convert these values to B-factors. Do this using the "Edit B-factor" 
   task in the "Coordinate data tools" menu. Choose the option based on the prediction 
   model type you are using. Once this conversion has been done, the models are ready 
   for use in Molecular Replacement (Phaser or Molrep).

      |Edit B-factors Gui|

#. **Search the PDB and EBI-AFDB using MrParse**

   The `EBI-AFDB <https://alphafold.ebi.ac.uk>`_ database from DeepMind and the EMBL-EBI, provides
   pre-calculated AlphaFold predictions of proteins, giving us access to many potentially homologous models to our
   target structure suitable for MR. 
   
   The MrParse task in the Bioinformatics menu can be used to 
   search for potential homologues for our target based on sequence alignment. It takes the AU definition as well
   as the reflection data as input and produces a html report displaying the list of
   hits from the PDB and EBI-AFDB as well as various associated scores for each hit. Users
   can then retrieve any of the models from the databases for use in MR (following the procedure in step 6). 
   
   The figure below shows the results of MrParse for our Gamma target. Colours for the PDB hits correspond to 
   sequence identity (red for high, green for low). For the EBI-AFDB hits, the colours indicate the predicted model 
   quality (pLDDT). They range from blue to indicate high confidence to orange for low confidence residues.

      |MrParse results|

#. **Interactive model preparation - Using the CCP4mg-MrBUMP Task to search and prepare models from the EBI-AFDB**

   -  The Interactive model preparation CCP4mg-MrBUMP task can be found in the Bioinformatics menu. Provide the 
      sequence loaded above and use default settings. By default, the task will identify homologues from the PDB 
      and the EBI-AFDB and display them in CCP4mg (takes a couple of minutes). You can search only the PDB or 
      EBI-AFDB by unticking the relevant box in the input menu. For AFDB models, you can adjust the
      pLDDT cut-off. Higher values mean only the most confidently predicted residues are 
      retained in the search models (**Note**: the CCP4mg-MrBUMP task will automatically convert
      pLDDT scores to B-factors for EBI-AFDB models).

      |CCP4MG-MRBUMP Input|

   -  After a couple of minutes the results of the MrBUMP sequence-based search should appear
      in the graphical interface. Note that the models will have been truncated to better
      match the target structure by referencing the sequence alignment between the model
      and the target. The `Sculptor <https://ccp4i2.gitlab.io/rstdocs/tasks/sculptor/index.html?highlight=sculptor>`_ 
      application is used to perform the truncation. The resulting search models are structurally aligned 
      using `Gesamt <https://ccp4i2.gitlab.io/rstdocs/tasks/gesamt/index.html?highlight=gesamt>`_. 

      |CCP4MG-MrBUMP results|

   -  The File ‑> Save all to CCP4i2 database menu option, will save all
      that is currently visible in CCP4mg and make it an "ensemble
      model" for potential use in a molecular replacement calculation
      within CCP4i2. 

   -  File ‑> Save all to CCP4i2 database will save all models, if all
      models are visible. To save a single model, make all other models
      invisible. To do this in CCP4mg, right-click (or Ctrl-click /
      Cmd-click on a laptop) on the coloured icon next to the model's
      name in the "Display Table" (right-hand panel of CCP4MG window.
      All the icons next to all the molecule names, except for the 
      selection from the AlphaFold prediction, `AF-P22892-F1 <https://alphafold.ebi.ac.uk/entry/P22892>`_, 
      have been right-clicked in the example below. The icon turns grey and
      the model becomes invisible in the graphics window. Right-clicking
      again makes the model visible again.

      |CCP4MG Display Table|

   -  CCP4mg offers many more sophisticated methods for choosing which
      happens to display, but this is beyond the scope of this tutorial.

   -  **Quit CCP4mg**

      The list of homologues will be displayed in the results
      window. You can select one with lower identity if you want to
      challenge Buccaneer.

#. **Molecular Replacement and Refinement PHASER Task**

   Load the data and Free reflection set, specify the AU contents (output from task above)
   and provide one of the the search model coordinates (models from AlphaFold, 
   RoseTTAFold or from the PDB/EBI-AFDB through CCP4MG-MrBUMP). 
   Search for one copy of the search model. For predicted models we recommend
   using an RMS of 1.2 Angstroms for the similarity to target setting.

#. **Refinement REFMAC Task**

   Run 30 cycles of refinement on the model from MR.

#. **Autobuild protein Buccaneer Task**

   Optionally, you can run the Buccaneer task to perform model building using
   the resulting map/phases from refinment. It should build an almost complete model.
   The output can be inspected in COOT and further cycles of refinement run.


.. |Edit B-factors Gui| image:: editbfactors.png
   :class: center
   :height: 350px
.. |CCP4MG-MRBUMP Input| image:: mg-input.png
   :class: center
   :height: 350px
.. |CCP4MG Display Table| image:: mg-gamma-displaytable.png
   :class: center
   :height: 350px
.. |CCP4MG-MrBUMP results| image:: mg-gamma-65.png
   :class: center
   :height: 350px
.. |MrParse results| image:: mrparse-output.png
   :class: center
   :height: 350px
