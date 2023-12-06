###########################
CCP4i2 Command Line - i2run
###########################

i2run is a mechanism for running CCP4i2 tasks and pipelines from the command line.
Tasks can be run so that they are recorded into a user's project database, recorded into *some other* database or executed without reference to a database

Why would I use i2run ?
=======================

A few scenarios

- Because you\'re in a hurry and can\'t be bothered to fire up CCP4i2
- Because you want to automate and do the same things many times
- Because you want to use the CCP4i2 database from a third party program (e.g. coot ?)

I'm interested...show me how easy it is
=======================================

Thank you.
Here is text that will run *refmac*, capturing the output in a (pre-existing) project called "GammaDemo" in your CCP4i2 database

.. code::

   i2run prosmart_refmac \
        --F_SIGF $CCP4/share/ccp4i2/demo_data/gamma/merged_intensities_Xe.mtz \
        --XYZIN $CCP4/share/ccp4i2/demo_data/gamma/gamma_model.pdb \
        --projectName GammaDemo
        
To create a new project (*e.g. GammaDemo3* ) to run your job, you will additionally have to specify a directory where the project should be created

.. code::

    i2run prosmart_refmac \
         --F_SIGF $CCP4/share/ccp4i2/demo_data/gamma/merged_intensities_Xe.mtz \
         --XYZIN $CCP4/share/ccp4i2/demo_data/gamma/gamma_model.pdb \
         --projectName GammaDemo3 \
         --projectPath /Users/martin/CCP4I2_PROJECTS/GammaDemo3
         
But why not just run (e.g.) *refmac* ?
======================================
***SO*** many reasons.

**Firstly** because the prosmart_refmac pipeline allows you to do much more.  In the same command line you can, for example, provide a reference PDB to which prosmart restraints will be calculated, and trivially add waters afterwards

.. code::

    # Generate some mdm2 data
    i2run aimless_pipe \
        --UNMERGEDFILES \
            crystalName=hg7 \
            dataset=DS1 \
            file=$CCP4/share/ccp4i2/demo_data/mdm2/mdm2_unmerged.mtz \
        --projectName hg7 \
        --projectPath /Users/martin/CCP4I2_PROJECTS/hg7
    #
    #Re-refine 4hg7.pdb, removing waters from start model, restrained against MDM2 for 4qo4.pdb
    i2run prosmart_refmac \
        --F_SIGF fileUse="-1.HKLOUT[0]" \
        --XYZIN \
            fullPath=$CCP4/share/ccp4i2/demo_data/mdm2/4hg7.pdb \
            selection/text="not (HOH)" \
        --prosmartProtein.TOGGLE True\
        --prosmartProtein.REFERENCE_MODELS \
            fullPath=$CCP4/share/ccp4i2/demo_data/mdm2/4qo4.cif \
        --projectName hg7

**Secondly** because your job is now organised without book keeping: one command line has sorted out how to store and name the output files so that they can readily be reused. It has also provided you with a nice digested report that is (arguably) even easier to analyse than a refmac log file !

         
Column labels
=============

Okay, so what if I want to run a job using an input MTZ file that does not use the standard mtz column labels used by CCP4i2 ?

Well, we\'ve got your back.  For example, if observed reflection data are stored in an MTZ file with  non canonical column labels, you just add a "columnLabels=" subdirective in the readily understood clipper syntax. In this case you specify the path to the input reflection file using the "fullPath=" directive

.. code::

    i2run phaser_simple \
        --F_SIGF \
            fullPath=$CCP4/share/ccp4i2/demo_data/beta_blip/beta_blip_P3221.mtz \
            columnLabels="/*/*/[Fobs,Sigma]" \
        --F_OR_I F \
        --XYZIN \
            $CCP4/share/ccp4i2/demo_data/beta_blip/beta.pdb \
        --projectName BetaBlip
        
*N.B.* The phaser task normally expects to find a set of intensities to use, unless the "F_OR_I" input parameter is set to the value "F"

Discovering the i2run command to achieve a particular run
=========================================================

Keywords for a task
--------------------

To find the keyword options available for a particular task or pipeline, simply issue an i2run command, specifying the name of the pipeline and giving the "-h" flag

.. code::

    i2run prosmart_refmac -h
    
From the CCP4i2 GUI
-------------------

Because the syntax to define a particular set of parameters for an i2 executable cna be complicated 
(which is, after all, why we provide a GUI), there is a mechanism to allow you to configure a run in the GUI
and discover the corrsponding i2run command.  Because this is a sophisticated option it is not exposed by default
and requires a user to "Right click" on the toolbar to control which buttons are visible.  Enable the button
"Show i2run command".  Now, your toolbar should include an extra button (Fig. 1)

   .. figure:: ShowI2runCommandButton.png
      :width: 50%

      Figure 1. Show i2run command button in CCP4i2 toolbar

With this option now available, you can configure a task to run in the GUI and, if the corresponding task is 
selected in the task list, clicking the "Show i2run command" button will pop up a modal dialog showing the exact
i2run command which would run the equivalent job (Fig.2)

   .. figure:: ShowI2runCommandResults.png
      :width: 50%

      Figure 2. i2run command for a task configured in the GUI

This text can be copied and pasted into a script to run the command from the command line. It can also be edited 
to modify, e.g. selecting different input files, or different control parameter values.

.. code::

    i2run phaser_simple --projectName GammaDemo \
        --F_SIGF \
            "dbFileId=fd96a5dc3ff511e99e2d6a00025fd280" \
            "contentFlag=1" \
            "baseName=merged_intensities_Xe_-imported_1.mtz" \
            "project=84e59a573ff511e99d256a00025fd280" \
            "subType=1" \
            "annotation=Reflections imported from merged_intensities_Xe_-imported.mtz by job 1" \
            "relPath=CCP4_IMPORTED_FILES" \
        --ASU_PROTEIN_MW 20000.0 \
        --ASU_NUCLEICACID_MW 1.0 \
        --ENSEMBLES \
            "use=True" \
            "number=1" \
            "label=Ensemble_1" \
        --FREERFLAG \
            "dbFileId=2d824287cd8011ea98779801a7ae2c3b" \
            "contentFlag=1" \
            "baseName=FREERFLAG.mtz" \
            "project=84e59a573ff511e99d256a00025fd280" \
            "annotation=45 FreeR - Spg:P 21 21 21;Resln:1.81A;Cell:34.1,54.8,68.0,90.0,90.0,90.0" \
            "relPath=CCP4_JOBS/job_45" \
        --RUNSHEETBEND True \
        --RUNREFMAC True \
        --XYZIN \
            "dbFileId=fd9b35de3ff511e98d276a00025fd280" \
            "contentFlag=1" \
            "baseName=gamma_model_1.pdb" \
            "project=84e59a573ff511e99d256a00025fd280" \
            "subType=0" \
            "annotation=Search imported from gamma_model.pdb by job 1" \
            "relPath=CCP4_IMPORTED_FILES" \
        --TITL SET] \
        --MUTE True \
        --MACA_MACR_REF_ANIS True \
        --MACA_MACR_REF_BINS True \
        --MACA_MACR_REF_SOLK True \
        --MACA_MACR_REF_SOLB True \
        --MACM_MACR_REF_ROTA True \
        --MACM_MACR_REF_TRAN True \
        --MACM_MACR_REF_BFAC True \
        --MACM_MACR_REF_VRMS True \
        --MACM_MACR_REF_CELL True \
        --MACM_MACR_REF_OFAC True \
        --MACM_MACR_LAST_ONLY True \
        --PURG_RNP_NUMB 0 \
        --RESC_ROTA True \
        --RESC_TRAN True \
        --SGAL_SELE all \

The command is relatively verbose, and the example above can easily be abridged to give the following `snippet`_.

.. _snippet:

.. code::

    i2run phaser_simple --projectName GammaDemo \
        --F_SIGF \
            "dbFileId=fd96a5dc3ff511e99e2d6a00025fd280" \
        --ASU_PROTEIN_MW 20000.0 \
        --ASU_NUCLEICACID_MW 1.0 \
        --ENSEMBLES \
            "use=True" \
            "number=1" \
            "label=Ensemble_1" \
        --FREERFLAG \
            "dbFileId=2d824287cd8011ea98779801a7ae2c3b" \
        --RUNSHEETBEND True \
        --RUNREFMAC True \
        --XYZIN \
            "dbFileId=fd9b35de3ff511e98d276a00025fd280" \
        --TITL SET] \
        --MUTE True \
        --MACA_MACR_REF_ANIS True \
        --MACA_MACR_REF_BINS True \
        --MACA_MACR_REF_SOLK True \
        --MACA_MACR_REF_SOLB True \
        --MACM_MACR_REF_ROTA True \
        --MACM_MACR_REF_TRAN True \
        --MACM_MACR_REF_BFAC True \
        --MACM_MACR_REF_VRMS True \
        --MACM_MACR_REF_CELL True \
        --MACM_MACR_REF_OFAC True \
        --MACM_MACR_LAST_ONLY True \
        --PURG_RNP_NUMB 0 \
        --RESC_ROTA True \
        --RESC_TRAN True \
        --SGAL_SELE all \

Input data objects
==================

Exotic dataobject types
=======================

CCP4i2 uses some special data classes, which present a challenge to specifying on the command line.  An example is the CEnsemble data structure which represents one or more models, each formed of one or more related sets of coordinate, to be used in moleculr replacement.  i2run handles these challenges by using subdirectives to address the different elements of the underlying data structure.  While hard to say, this is relatively easy to illustrate.  Here is an i2run command which specifies a set of ensembles to use in the phaser expert task to solve beta blip

.. code ::

    i2run phaser_pipeline \
        --F_SIGF \
            fullPath=$CCP4/share/ccp4i2/demo_data/beta_blip/beta_blip_P3221.mtz \
            columnLabels="/*/*/[Fobs,Sigma]" \
        --ENSEMBLES \
            use=True \
            pdbItemList/identity_to_target=0.9 \
            pdbItemList/structure=$CCP4/share/ccp4i2/demo_data/beta_blip/beta.pdb \
        --ENSEMBLES \
            use=True \
            pdbItemList/identity_to_target=0.9 \
            pdbItemList/structure=$CCP4/share/ccp4i2/demo_data/beta_blip/blip.pdb \
        --RUNREFMAC False \
        --RUNSHEETBEND False \
        --F_OR_I F \
        --RESOLUTION_HIGH 3.0 \
        --ASUFILE \
            seqFile=$CCP4/share/ccp4i2/demo_data/beta_blip/beta.seq \
            seqFile=$CCP4/share/ccp4i2/demo_data/beta_blip/blip.seq \
        --projectName BetaBlip

Specifying file inputs
======================

We have already introduced the option of specifying input files by their full path, *e.g.*

.. code::

    i2run coordinate_selector \
        --XYZIN fullPath=$CCP4/share/ccp4i2/demo_data/CDK1CyclinBCKS2/1jst.pdb \
        --projectName CDK1CyclinBCKS2
    
or more succinctly:

.. code::

    i2run coordinate_selector \
        --XYZIN $CCP4/share/ccp4i2/demo_data/CDK1CyclinBCKS2/1jst.pdb \
        --projectName CDK1CyclinBCKS2

But it is also possible to use a file that already exists within the database by using a reference to a job number and parameter name where it was created/used with the syntax *jobNumber*. *paramName* :

.. code::

    i2run coordinate_selector \
        --XYZIN fileUse=1.XYZIN \
        --projectName CDK1CyclinBCKS2

A *negative* jobNumber indicates a *relative* job number, *e.g.* a jobNumber of -1 refers to the previous job in the project

If the file was an element of an output list (as is the case for tasks which can produce multiple outputs, such as multiple datasets output from AIMLESS), then you can use array-style indexing to specify which element to choose:

.. code ::

    i2run import_merged \
        --HKLIN fullPath=$CCP4/share/ccp4i2/demo_data/gamma/merged_intensities_Xe.mtz \
        --projectName GammaDemo
    
    i2run aimless_pipe --UNMERGEDFILES \
        crystalName=Nat1 dataset=DS1 file=$CCP4/share/ccp4i2/demo_data/gamma/gamma_native.mtz \
        --projectName GammaDemo
        
    i2run shelx \
        --SHELXCDE True \
        --ATOM_TYPE Xe \
        --EXPTYPE SIRAS \
        --NATIVE True \
        --F_SIGFanom fileUse="-2.OBSOUT" \
        --F_SIGFnative fileUse="-1.HKLOUT[0]" \
        --SEQIN seqFile=$CCP4/share/ccp4i2/demo_data/gamma/gamma.pir \
        --USER_MBREF_BIGCYC True \
        --MBREF_BIGCYC 2 \
        --projectName GammaDemo

Or by using the files database UUID, which may be discovered by inspecting the database entry of the file in an interactive CCP4i2 session

Coordinates and selections
==========================

A subset of CCP4i2 tasks \(*e.g.* list\) allow for the task to run on a specified subset of an input coordinate file.  To specify this subset on the command line use the "selection/text=" subdirective on a coordinate file specification. Selections are provided in the mmdb syntax:

.. code::

    i2run coordinate_selector \
        --XYZIN \
            fullPath=$CCP4/share/ccp4i2/demo_data/CDK1CyclinBCKS2/1jst.pdb \
            "selection/text=A/ or (CYS)" \
        --projectName CDK1CyclinBCKS2
        
Asymmetric unit content files
=============================

Several CCP4i2 tasks \(*e.g.* phaser\, buccaneer, parrot) use a special file (file.asucontent.xml) to specify the expected contents of the asymmetric unit. These files can be created by a dedicated task using i2run...

.. code::

    i2run ProvideAsuContents \
        --ASU_CONTENT \
            description=Beta-lactamase \
            sequence=HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW \
            source/baseName=beta.seq \
            source/relPath=$CCP4/share/ccp4i2/demo_data/beta_blip \
            nCopies=1 \
            polymerType=PROTEIN \
            name=Beta-lactamase \
        --ASU_CONTENT \
            description="Beta-lactamase inhibitory protein" \
            sequence=AGVMTGAKFTQIQFGMTRQQVLDIAGAENCETGGSFGDSIHCRGHAAGDYYAYATFGFTSAAADAKVDSKSQEKLLAPSAPTLTLAKFNQVTVGMTRAQVLATVGQGSCTTWSEYYPAYPSTAGVTLSLSCFDVDGYSSTGFYRGSAHLWFTDGVLQGKRQWDLV \
            source/baseName=blip.seq \
            source/relPath=$CCP4/share/ccp4i2/demo_data/beta_blip \
            nCopies=1 \
            polymerType=PROTEIN \
            name=Beta-lactamase-inhibitor \
        --projectName BetaBlip \
        --projectPath /tmp/martin/BetaBlip \
        --dbFile /tmp/martin/BetaBlip.sqlite

Alternatively, i2run provides a mechanism to generate them "on the fly" in tasks that use them:

.. code::

    i2run parrot\
        --F_SIGF fullPath=$CCP4/share/ccp4i2/demo_data/gamma/merged_intensities_Xe.mtz \
        --ABCD fullPath=$CCP4/share/ccp4i2/demo_data/gamma/initial_phases.mtz \
        --ASUIN \
            seqFile=$CCP4/share/ccp4i2/demo_data/gamma/gamma.pir \
        --projectName GammaDemo


