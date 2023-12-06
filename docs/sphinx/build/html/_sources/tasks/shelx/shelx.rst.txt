####################################
SHELXC/D/E experimental phasing task
####################################

Input Data
==========

   The "Automated structure solution - Shelxcde" experimental phasing
   task builds a model automatically from intensities or amplitudes from
   single or multiple-wavelength anomalous diffraction (SAD/MAD) or
   single isomorphous replacement (SIRAS) experiments. The user can
   optionally start or end the automated pipeline at each of the steps
   below:

   #. Substructure detection with
      `SHELXD <http://shelx.uni-ac.gwdg.de/SHELX/shelxd_mm_keywords.php>`__
      using FA values obtained from
      `SHELXC <http://shelx.uni-ac.gwdg.de/SHELX/shelxc_keywords.php>`__
   #. Density modification and poly-alanine tracing with
      `SHELXE <http://shelx.uni-ac.gwdg.de/SHELX/shelxe_keywords.php>`__
   #. Model building with
      `Buccaneer <http://www.ccp4.ac.uk/html/cbuccaneer.html>`__, or
      `ARP/wARP <http://www.embl-hamburg.de/ARP/Manual/UserGuide7.5.html>`__
      iterated with model refinement with
      `REFMAC <http://www.ccp4.ac.uk/html/refmac5.html>`__
   #. Model refinement with REFMAC

Input
-----

   .. figure:: shelxcde-input-labels.png
      :alt: Figure 1: Main input

      Figure 1: Main input

   A pipeline can be started or stopped at any stage **(1.1)** For
   example, if an anomalously scattering substructure has been
   determined already, it can be input in PDB format by selecting to
   start the pipeline at Density modification and poly-Ala tracing.

   Inputting a protein sequence **(1.2)** is beneficial for model
   building. If it is not available, click off the Input protein
   sequence option and input the number of protein residues per monomer.

   Input the substructure atom element symbol (i.e. Se for Selenium) and
   the number of substructure atoms in the asymmetric unit **(1.3)**.
   The number of substructure atoms in the asymmetric unit is the number
   of substructure atoms in your protein multiplied by the number of
   protein molecules that are present in the asymmetric unit. For
   sulfur-SAD phasing of structures containing disulfide bonds, the
   sulfur atoms that make such bonds should be counted separately in the
   number of substructure sites, and the minimum distance between sites
   should be given as 1.5A in the "Advanced options" section (3.3). If
   the anomalous resolution is worse than about 2.1A, the disulfides
   will not be resolved, so the number of disulfides (i.e. "S-S pairs")
   should be specified and the minimum distance increased to 3.5A, but
   the number of sites is still the total number of anomalously
   scattering atoms.

   Input the anomalous reflections data in a merged MTZ file containing
   the Bijvoet/anomalous intensities or structure factor amplitudes
   **(1.5)**. Alternatively, unmerged or merged XDS (XDS_ASCII.HKL),
   HKL2000 (.sca) or SHELX (.hkl) files can be input **(1.4)**.
   Inputting unmerged files is recommended if SHELXC is used since it
   can make use of it for improved estimates of FA values and reporting
   of statistics such as the anomalous correlation coefficient
   (CCanom1/2). If non-MTZ reflection data is input, fields for crystal
   cell dimensions and spacegroup appear and are automatically filled by
   the ccp4i2 interface if the information exists in the input file. If
   they have not been filled in, please input them. For MAD data, the
   reflection data can be input for all the other wavelengths collected
   **(1.6)**.

   If an isomorphic native data set is available, it should be input
   **(1.7)** as it typically extends to higher resolution. The native
   dataset is always used for building and refinement but it can be also
   used for SIRAS phasing: If only a single anomalous dataset is input
   together with the native data, an option to switch between the SIRAS
   and SAD phasing can be selected under the "Important Paramaters" tab.
   Finally, a previously defined cross-validation or "free" set can be
   input **(1.8)** or the cross-validation can be turned off, otherwise
   the default is to define a new set

Important Options
-----------------

   .. figure:: crank2-important-options.png
      :alt: Figure 2: Important options

      Figure 2: Important options

   This page (and the next one) offers options that can improve results
   if default (automatically determined) values are not optimal. From
   the number of protein residues per monomer, the interface will obtain
   a guess for the number of molecules in the asymmetric unit and the
   solvent content based on a Matthew's coefficient analysis. These
   values are only a guess and (re-)running the task with 'correct'
   values can significantly improve results.

   .. figure:: shelxcde-advanced-options-labels.png
      :alt: Figure 3: Advanced options

      Figure 3: Advanced options

Advanced Options
----------------

   The third pane offers advanced options, ones you will not normally
   need to change for the first run. If you suspect that substructure
   detection was not successful, consider running different resolution
   cutoffs **(3.1)** and more trials **(3.2)** if the results indicate
   that there is an anomalous signal. For sulfur-SAD phasing, as
   mentioned above, if your anomalous resolution is greater than 2.1,
   adjust the minimum distance between atoms **(3.3)** to 1.5. Another
   useful option that can improve results is increasing the number of
   cycles of model building in SHELXE\ **(3.4)** or adjusting the
   solvent content and inputting the correct NCS copies in the
   "Important options" section. If the substructure is correct, it is
   faster to re-run SHELXE with different options from the "SHELX
   phasing and building" pipeline starting at Density modification and
   poly-Ala tracing. The partial model and initial phases can also be
   input to the "CRANK2 phasing & building" pipeline and starting at
   Substructure improvement. For all steps, the gui provides the ability
   to input custom keywords for individual programs (i.e. **(3.5)**,
   **(3.6)**). For example, inputting -t5 to SHELXE for predominantly
   beta sheet proteins can provide improved tracing, but also increase
   run time. Please look at all the documentation for individual
   programs (links shown above) to see supported keywords. Multiple
   options can be inputted, separated by a comma.

Results
=======

   While the task is running, graphs are updated to indicate the
   performance for each step. After and typically also while each step
   is being performed, a summary with statistics indicating the
   performance is provided. For many graphs, a shaded gray area is shown
   to indicate a region with poor statistics. For some steps, the user
   may click on a button on the bottom right panel of the Results
   section to stop the current step. If the button is pressed, the step
   will be stopped after finishing the current iteration and the next
   step of the pipeline will follow. This can be used to save time if
   the user is happy with the results of the current step as indicated
   by the Results report.

FA estimation
-------------

   .. figure:: danosigdano.png
      :alt: Figure 4: Anomalous difference to noise versus resolution

      Figure 4: Anomalous difference to noise versus resolution

   Besides generating the FA values, this step may also provide useful
   information about the input datasets. Figure 4 shows the graph of
   Dano/Sigdano - a measure of the anomalous signal strength - as a
   function of resolution. The larger the number, the greater the signal
   and a value of about 0.8 indicates zero signal. If unmerged data was
   input, the graph of CCAnom1/2 - the anomalous difference correlation
   within a dataset - versus resolution (Figure 5)

   .. figure:: anom-cchalf.png
      :alt: Figure 5: Anomalous difference correlation within dataset

      Figure 5: Anomalous difference correlation within dataset

   is shown. Again the larger the number, the greater the signal and
   values over approximately 30% indicate a significant anomalous
   signal. Evaluating these graphs is useful if the substructure was not
   correctly identified and the user should re-run the task by manually
   inputting a resolution cutoff, or even several jobs varying the
   cutoff for difficult cases with small anomalous signals. Similarly,
   in the case of a MAD experiment, the signed anomalous correlation
   between data sets is shown (Figure 6)

   .. figure:: signed-anom-diff-correl-MAD.png
      :alt: Figure 6: Anomalous difference correlation between MAD datasets

      Figure 6: Anomalous difference correlation between MAD datasets

   If there is a significant anomalous signal for different
   wavelengths, a high anomalous difference correlation for them will be
   present (above 30%).

Substructure detection
----------------------

   .. figure:: shelxd-occ-cc-ccweak.png
      :alt: Figure 7: Summary, Occupancy and CC versus CCweak for the trials graphs

      Figure 7: Summary, Occupancy and CC versus CCweak for the trials
      graphs

   .. figure:: cfom-trials.png
      :alt: Figure 8: Distribution of CFOM from trials

      Figure 8: Distribution of CFOM from trials

   A good indication that the substructure has been correctly determined
   is if the threshold (i.e. CFOM in the case of SHELXD) has been
   reached before the maximum number of cycles has been performed: the
   report summary for the Substructure detection step will state this. A
   division of trials in the Distribution of CFOM (Figure 7) is also an
   encouraging indicator: the graph shows the separation of incorrect
   (lower CFOM) trials from correct (higher CFOM) solutions. This
   division can also be observed as clustering in a graph of CC against
   CCweak where the trials (shown in circles) in the the top right
   indicate correct solutions with higher values of CC and CCweak.
   However, some datasets, especially those with small anomalous
   signals, show little discrimination between the CFOM of solutions and
   of wrong substructures. The graph of occupancies versus atom number
   details the number of significant peaks that have been found and can
   be used cautiously to indicate the number of substructure atoms and
   protein molecules in the asymmetric unit.

Density modification and poly-alanine trace
-------------------------------------------

   .. figure:: shelxe-summary-cc.png
      :alt: Figure 9: Summary and correlation coefficient of enantiomorphs

      Figure 9: Summary and correlation coefficient of enantiomorphs

   The results and summary of running density modification and initial
   model building for both substructure enantiomorphs is presented. A
   significant difference in between both hands for the correlation
   coefficient in model building (Figure 9),

   .. figure:: shelxe-contrast.png
      :alt: Figure 11: Contrast between enantiomorphs

      Figure 11: Contrast between enantiomorphs

   contrast (Figure 11) and

   .. figure:: shelxe-residues.png
      :alt: Figure 10: Residues built for enantiomorphs

      Figure 10: Residues built for enantiomorphs

   residues built (Figure 10),

   .. figure:: shelxe-chain-length.png
      :alt: Figure 12: Chain lengths of enantiomorphs

      Figure 12: Chain lengths of enantiomorphs

   chain length (Figure 12) are promising sign. The most important
   number output by shelxe is the CC for the traces against the native
   data. If this is greater than 25%, the structure is probably solved.
   This is still true even if a totally incorrect protein sequence has
   been provided!

Model building
--------------

   .. figure:: residues-per-cycle-and-rcomb.png
      :alt: Figure 13: Summary and graphs of Residues, FOM and Rfactors

      Figure 13: Summary and graphs of Residues, FOM and Rfactors

   The model building summary (Figure 13) will indicate whether a
   majority of the model was correctly built. In the graphs, an R-factor
   under 40% and a large number of residue built in a small number of
   fragments are good indication of success.

Model refinement
----------------

   .. figure:: refine-refmac.png
      :alt: Figure 14: Summary and graphs of Rfactors

      Figure 14: Summary and graphs of Rfactors

   Finally, the model refinement summmary (Figure 14) quotes the
   R-factors versus the cycle number: lower R-factors (underneath 40%)
   are a good indication and a small difference between the working and
   free R-factors is positive. In general, the R-values should plateau -
   if this hasn't occurred, and the R-factors are still going down, more
   refinement steps can further improve the model.
