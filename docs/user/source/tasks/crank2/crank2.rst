################################
Crank2 experimental phasing task
################################

   The "Automated structure solution - Crank2" experimental phasing task
   builds a model automatically from intensities or amplitudes from
   single or multiple-wavelength anomalous diffraction (SAD/MAD) or
   single isomorphous replacement (SIRAS) experiments. The user can
   optionally choose to input a partial model and start or end the
   automated pipeline at each of the steps below:

   #. Substructure detection with
      `SHELXD <https://shelx.uni-goettingen.de/shelxd_mm_keywords.php>`__
      or PRASA (to be released in 2016) using FA values obtained from
      `SHELXC <https://shelx.uni-goettingen.de/shelxc_keywords.php>`__,
      `AFRO <http://www.ccp4.ac.uk/html/afro.html>`__ or
      `ECALC <http://www.ccp4.ac.uk/html/ecalc.html>`__
   #. Substructure improvement & phasing with
      `REFMAC <http://www.ccp4.ac.uk/html/refmac5.html>`__,
      `SHELXE <https://shelx.uni-goettingen.de/shelxe_keywords.php>`__,
      or `BP3 <http://www.ccp4.ac.uk/html/bp3.html>`__
   #. Hand determination with
      `SOLOMON <http://www.ccp4.ac.uk/html/solomon.html>`__ and MAPRO or
      SHELXE
   #. Density modification with SHELXE or
      `PARROT <http://www.ccp4.ac.uk/html/parrot.html>`__ or SOLOMON
      with phase combination with REFMAC or MULTICOMB
   #. Model building with
      `Buccaneer <http://www.ccp4.ac.uk/html/cbuccaneer.html>`__,
      `ARP/wARP <http://www.embl-hamburg.de/ARP/Manual/UserGuide7.5.html>`__
      or SHELXE iterated with model refinement with REFMAC
   #. Model refinement with REFMAC

Input Data
==========

   .. figure:: crank2-input-labels.png
      :alt: Figure 1: Main input

      Figure 1: Main input

   A pipeline can be started or stopped at any stage **(1.1)** For
   example, if an anomalously scattering substructure has been
   determined already, it can be input in PDB format by selecting to
   start the pipeline at "Substructure improvement" (SAD experiment) or
   "Substructure phasing" (SIRAS,MAD).

   A partial protein model can be input in PDB format **(1.2)**, to be
   rebuilt using the SAD data. If the partial model is biased, the
   option to exclude the non-anomalous atoms from the partial model can
   give better results. It is recommended to try both approaches.

   Inputting a protein sequence **(1.3)** is beneficial for model
   building. If it is not available, click off the Input protein
   sequence option and input the number of protein residues per monomer.

   Input the substructure atom element symbol (i.e. Se for Selenium) and
   the number of substructure atoms in the asymmetric unit **(1.4)**.
   The number of substructure atoms in the asymmetric unit is the number
   of substructure atoms in your protein multiplied by the number of
   protein molecules that are present in the asymmetric unit. For
   sulfur-SAD phasing of structures containing disulfide bonds, the
   sulfur atoms that make such bonds should be counted separately in the
   number of substructure sites, and the minimum distance between sites
   should be given as 1.5A in the "Advanced options" section **(3.3)**.
   If the anomalous resolution is worse than about 2.1A, the disulfides
   will not be resolved, so the number of disulfides (i.e. "S-S pairs")
   should be specified and the minimum distance increased to 3.5A, but
   the number of sites is still the total number of anomalously
   scattering atoms.

   Input the anomalous reflections data in a merged MTZ file containing
   the Bijvoet/anomalous intensities or structure factor amplitudes
   **(1.6)**. Alternatively, unmerged or merged XDS (XDS_ASCII.HKL),
   HKL2000 (.sca) or SHELX (.hkl) files can be input **(1.5)**.
   Inputting unmerged files is recommended if SHELXC is used since it
   can make use of it for improved estimates of FA values and reporting
   of statistics such as the anomalous correlation coefficient
   (CCanom1/2). If non-MTZ reflection data is input, fields for crystal
   cell dimensions and spacegroup appear and are automatically filled by
   the ccp4i2 interface if the information exists in the input file. If
   they have not been filled in, please input them.

   Finally, the Crank2 pipeline requires the anomalous scattering
   coefficients f' and f" for the substructure atom at the wavelength of
   the data collected. The interface can automatically fill in the
   theoretical values for the scattering factors but it is very
   important to check these, especially for wavelengths close to the
   absorption edge. The values from a fluorescent scan are always better
   and recommended to use. For MAD data, the reflection data and f' and
   f" values can be input for all the other wavelengths collected
   **(1.7)**.

   If a native data set is available, it can also be input **(1.8)**.
   The native dataset is always used for building and refinement but it
   can be also used for SIRAS phasing: If only a single anomalous
   dataset is input together with the native data, an option to switch
   between the SIRAS and SAD phasing can be selected under the
   "Important Paramaters" tab. Finally, a previously defined
   cross-validation or "free" set can be input **(1.9)** or the
   cross-validation can be turned off, otherwise the default is to
   define a new set

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

   .. figure:: crank2-advanced-options-labels.png
      :alt: Figure 3: Advanced options

      Figure 3: Advanced options

Advanced Options
----------------

   The third pane offers advanced options, ones you will not normally
   need to change for the first run. If you suspect that substructure
   detection was not successful, consider running with different
   resolution cutoffs **(3.1)** and more trials **(3.2)** if the results
   indicate that there is an anomalous signal. For sulfur-SAD phasing,
   as mentioned above, if your anomalous resolution is greater than 2.1,
   adjust the minimum distance between atoms **(3.3)** to 1.5. Another
   useful option that can improve results is more iterations in model
   building.\ **(3.4)** If you wish to improve a partial model, it is
   generally better to "continue" the previous job by inputting the
   partially built model and phases from it by using the "CRANK2 phasing
   & building" button that appears below the report from a completed
   job. For all steps, the gui provides the ability to input custom
   keywords for individual programs (i.e. **(3.5)**, **(3.6)**). Please
   look at all the documentation for individual programs (links shown
   above) to see supported keywords. Multiple options can be inputted,
   separated by a comma.

Results
=======

   While the task is running, graphs are updated to indicate the
   performance for each step. After and typically while each step is
   being performed, a summary with statistics indicating the performance
   is provided. For many graphs, a shaded gray area is shown to indicate
   a region of poor statistics. For some steps, the user may click on a
   button on the bottom right panel of the Results section to stop the
   current step. If the button is pressed, the step will be stopped
   after finishing the current iteration and the next step of the
   pipeline will follow. This can be used to save time if the user is
   happy with the results of the current step as indicated by the
   Results report.

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
^^^^^^^^^^^^^^^^^^^^^^

   .. figure:: shelxd-occ-cc-ccweak.png
      :alt: Figure 7: Summary, Occupancy and CC versus CCweak for the trials graphs

      Figure 7: Summary, Occupancy and CC versus CCweak for the trials graphs

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

Substructure improvement and phasing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   .. figure:: phasing-summary.png
      :alt: Figure 9: Phasing summary

      Figure 9: Phasing summary

   A Crank2 job from SAD data will attempt to automatically improve the
   substructure obtained from phasing, using anomalous likelohood
   gradient maps for detection of additional anomalous scatterers. The
   phasing report summary shows the number of anomalous atoms that have
   been added or removed and quotes the estimated mean figure of merit
   (FOM) of phases. For substructure phasing, density modification and
   model building and refinement, the mean FOM is the estimated mean
   cosine of the phase error: thus, a number closer to 1 is better than
   a number closer to zero. Values greater than 0.3 after phasing often
   indicate good phases. Yet, the average FOM is just an estimation of
   the quality of the phases which may be skewed and occasionally models
   can be successfully built with a low estimated FOM after phasing: for
   example, this data set has a figure of merit of 0.1 after phasing,
   but still leads to successful model building.

Hand determination
^^^^^^^^^^^^^^^^^^

   .. figure:: hand-det-crank2.png
      :alt: Figure 10: Hand determination summary

      Figure 10: Hand determination summary

   Reciprocal space information from substructure phasing can not, in
   general, resolve the correct enantiomorph. However, in many cases
   comparison of the quality of resulting experimental maps for both
   hands can be used to distinguish the correct hand. The quality of the
   experimental map can be estimated by its skewness or a related
   parameter of correlation with local deviation of the map (CLD).
   Furthermore, statistics from density modification of both maps, such
   as FOM and contrast, can be used to determine the hand, with the
   assumption that the right hand can be improved more than the wrong
   one.
   The Crank2 pipeline uses a combined score from comparison between the
   CLD from phasing and FOM from iterations of a quick density
   modification. In vast majority of cases, the right hand is correctly
   identified by the scoring and a large difference between the scores
   indicates that a correct hand has been chosen (thus, it may be also
   used as an indirect indicator of success of the substructure
   determination). However, in case of very weak experimental phases it
   may happen that the map based indicators fail and the wrong hand is
   picked. Such case can be often recognized by small discrimination
   between the scores and the CLD values of both hands around 0.
   Although such cases are rate, it is recommended to rerun the job with
   the other hand
   The SHELX pipeline does not attempt to determine the hand; instead,
   it runs SHELXE density modification and building with both hands and
   the correct hand is simply the one that leads to a model built.
   However, the contrast of the map for both hands in the initial
   density modification is reported in a plot (Figure X+2); a
   significantly better contrast such as shown in the Figure X+2 is a
   good indicator of the correct hand.

Density modification
^^^^^^^^^^^^^^^^^^^^

   .. figure:: denmod-summary-crank2.png
      :alt: Figure 11: Density modification summary

      Figure 11: Density modification summary

   Figure 11 shows the summary from density modification, reporting if
   non-crystallographic symmetry (NCS) operators were detected by Parrot
   and the figure of merit (FOM) and an estimate of the bias that was
   detected in the density modification (beta parameter). A value of FOM
   underneath 0.35 indicates either very weak phases or a wrong
   solution, ie the substructure (or rarely its hand) was not determined
   correctly. A value of FOM above 0.5 indicates good phases with high
   chances of model building success in the next step (except of
   datasets with resolution lower than ~4A where the resolution may
   become the model building bottleneck).

Model building
^^^^^^^^^^^^^^

   .. figure:: residues-per-cycle-and-rcomb.png
      :alt: Figure 12: Summary and graphs of Residues, FOM and Rfactors

      Figure 12: Summary and graphs of Residues, FOM and Rfactors

   The model building summary (Figure 12) will indicate whether a
   majority of the model was correctly built. In the graphs, an R-factor
   under 40% and a large number of residue built in a small number of
   fragments are good indication of success.

Model refinement
^^^^^^^^^^^^^^^^

   .. figure:: refine-refmac.png
      :alt: Figure 13: Summary and graphs of Rfactors

      Figure 13: Summary and graphs of Rfactors

   Finally, the model refinement summmary quotes the R-factors versus
   the cycle number: lower R-factors (underneath 40%) are a good
   indication and a small difference between the working and free
   R-factors is positive. In general, the R-values should plateau - if
   this hasn't occurred, and the R-factors are still going down, more
   refinement steps can further improve the model.
