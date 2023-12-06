CCP4i2 MR Tutorial: Apo-CeuE (PDB:3ZKW)
=======================================

A medium sized structure with a MW of ~30kDa which can be solved by
molecular replacement. There are three independent monomers of apo-CeuE
in the P1 asymmetric unit. The solution is slower than for the Getting
Started set, but shows the power of Buccaneer in building larger models.
The structure was published by Raines, D.J., Moroz, O.V., Wilson, K.S.,
and Duhme-Klair, A. "Interactions of a Periplasmic Binding Protein with
a Tetradentate Siderophore Mimic". Angew Chem Int Ed Engl (2013) 52,
4595-8", PubMed ID
`23512642 <http://europepmc.org/article/MED/23512642>`__, DOI
`10.1002/anie.201300751 <http://dx.doi.org/10.1002/anie.201300751>`__.

Screenshots are not extensively included in this script - examples can
be found in the Getting Started scripts if needed.

#. **Data Reduction - Aimless Task**

   Load the unmerged XDS data apo-ceue-sad-sweep1.hkl from the CeuE Demo
   Data folder. The structure was deposited with data to 1.7Å, but from
   the statistics the data to 1.6Å can be used in the tutorial.

#. **Import Sequence(s) Task**

   Import the 3zkw.seq sequence from the CeuE Demo Data folder.

#. **Define AU Contents Task**

   Load the imported 3zkw sequence and the merged data. The Task should
   indicate the most likely AU content is three molecules of CeuE.

#. **Interactive model preparation - CCP4mg and MrBUMP Task**

   -  Provide the sequence loaded above and use default settings. The
      Task will identify homologues from the PDB and display them in
      CCP4mg. (Takes a couple of minutes). **Do not select the identity
      structure, 3zkw**. 4jcc or 3gfv are good search models with
      identities just over 40%.

   -  The File ‑> Save all to CCP4i2 database menu option, will save all
      that is currently visible in CCP4mg and make it an "ensemble
      model" for potential use in a molecular replacement calculation
      within CCP4i2.

   -  File ‑> Save all to CCP4i2 database will save all models, if all
      models are visible. To save a single model, make all other models
      invisible. To do this in CCP4mg, right-click (or Ctrl-click /
      Cmd-click on a laptop) on the coloured icon next to the model's
      name in the "Display Table" (right-hand panel of CCP4MG window.
      All the icons next to all the molecule names, except 1a3h have
      been right-clicked in the example below. The icon turns grey and
      the model becomes invisible in the graphics window. Right-clicking
      again makes the model visible again.

      |CCP4MG Display Table|

   -  CCP4mg offers many more sophisticated methods for choosing which
      happens to display, but this is beyond the scope of this tutorial.

   -  **Quit CCP4mg**

      and the list of homologues will be displayed in the results
      window. You can select one with lower identity if you want to
      challenge Buccaneer.

      |CCP4MG / MrBUMP results|

      Alternatively, you can skip this step and use the *Import a
      Coordinates set Task* to download either the 4jcc or 3gfv
      coordinates directly from the EBI-PDB - remember to select Protein
      and Chain A. This is not quite the same as using CCP4mg which
      sculpts and aligns the model to match your sequence as best as
      possible.

#. **Molecular Replacement and Refinement MOLREP Task**

   Load the data and Rfree, the search model coordinates and the CeuE
   sequence. This should find three copies in the AU and run 10 cycles
   of refinement. You can alternatively run the *Basic Molecular
   Replacement - Phaser Task*. For Phaser, if you have downloaded the
   sculpted model from the CCP4mg-MrBUMP Task you can directly provide
   this to Phaser (or indeed an ensemble of models). If you are using a
   complete model from the PDB, use the *Truncate search model -
   Sculptor Task* to trim the model. Either MOLREP or Phaser will solve
   the structure easily.

#. **Refinement REFMAC Task**

   Run 25 cycles of refinement on the model from MR.

#. **Autobuild protein Buccaneer Task**

   Buccaneer will now build an almost complete model for all three
   chains. The output can be inspected in COOT and further cycles of
   refinement run.

#. **ACORN - Phase refinement with Dynamic Density Modification Task**

   This can be run before the Autobuild Task since the data extend to
   1.6Å resolution. ACORN will produce a set of phases which can be
   treated as experimental phases in Buccaneer and give an even better
   starting point than the refined MR model. ACORN also eliminates any
   residual model bias. Now run Buccaneer.

.. |CCP4MG Display Table| image:: image1.png
   :class: center
   :height: 350px
.. |CCP4MG / MrBUMP results| image:: image2.png
   :class: center
   :height: 350px
