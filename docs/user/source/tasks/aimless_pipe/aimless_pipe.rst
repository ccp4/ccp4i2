####################################
The AIMLESS Data reduction pipeline
####################################

   The Data Reduction task takes reflection intensities from an
   integration program such as Mosflm, XDS or DIALS, and converts them
   to scaled merged intensities and structure amplitudes. This consists
   of a number of stages

   #. Read one or more input files,  combining them and putting them on a consistent indexing scheme if necessary (program POINTLESS)
   #. Determine the crystal symmetry, ie the Laue group and if possible the space group, by inspecting the internal symmetry of the intensities (program POINTLESS)
   #. Scale together symmetry-related reflections (program AIMLESS)
   #. Analyse the merged data for information content, twinning, translational non-crystallographic symmetry (program PHASER)
   #. Examine the intensity statistics to detect twinning, and generate amplitudes (F) from intensities (I) (program CTRUNCATE)
   #. Generate or extend a set of FreeR reflections for later use in refinement

   If the option to automatically cut the resolution range is selected
   (see below **(1.8)** ), then the results from Phaser are examined, and
   if necessary AIMLESS is rerun with a new resolution cutoff, and then
   the pipeline continues with PHASER etc as before

   Detailed documentation for individual programs is available here:

   `POINTLESS <https://ccp4.gitlab.io/ccp4docs/pointless.html>`__
   `AIMLESS <https://ccp4.gitlab.io/ccp4docs/aimless.html>`__
   `CTRUNCATE <https://ccp4.gitlab.io/ccp4docs/ctruncate.html>`__

How to read the Results_


Input
=====

.. figure:: aimless_pipe_task_1.png
   :alt: Figure 1: Main input

   Figure 1: Main input

.. rubric:: Reflection input
   :name: reflection-input

Input reflection files are specified here, either from earlier jobs,
or from outside the GUI. They should contain unmerged reflections, ie
the output from an integration program **(1.1)** Click "Show list" to
allow input of multiple files. File formats recognised are: unmerged
MTZ, from Mosflm or DIALS; files from XDS, XDS_ASCII.HKL (scaled) or
INTEGRATE.HKL (unscaled); .raw files from SAINT (Bruker); Scalepack
files (.sca); and ShelX files. The last two cannot be properly scaled
by AIMLESS as they lack diffraction geometry information (only the
CONSTANT scale option is allowed, only useful to scale together
multiple files). Merged MTZ files can also be used, but again can only
be scaled to merge multiple files (use SCALE CONSTANT, see **(3.2)** below).

Each file must be identified by a crystal name and dataset name:
these should be short and not include spaces, and must be set here if
they haven't been set by the integration program **(1.2)** . They are
used to identify the dataset in later structure solution tasks (as
/crystal/dataset/). Files with the same dataset name will be merged
into a single output dataset: different dataset names (eg for MAD
data) will lead to multiple output sets. Parts of the data may be
omitted either by specifying batch ranges for exclusion **(1.3)**, or
by setting resolution limits **(1.4)**. Typically you would do that
after running the task first and inspecting the output. It is often
useful to integrate the images to slightly higher resolution than you
think the pattern extends to, then cut back the resolution following
a first run of Aimless (this can be done automatically, see later).
And parts of the data showing serious
radiation damage may be omitted, provided that completeness is not
compromised (you need to check that yourself). 

The default option is to try to determine the Laue group (the
rotational point group), and if possible the space group. However if
you know the space group you can specify it explicitly, selecting the
option from the pull-down menu **(1.5)**. The space group may also be
set by giving a reference file with the correct space group
**(1.6)**. The reference file is also used to resolve indexing
ambiguities in space groups where there are non-equivalent but valid
alternative indexings (the same cases where merohedral or
pseudo-merohedral twinning is possible). The reference file may be a
previous isomorphous dataset (intensities or Fs), or atomic
coordinates from which structure factors are calculated. If multiple
input files are given, then if necessary they will be checked for
consistent indexing, against the reference data if given, or against
the first file.

The last item on this page can be used to specify a pre-existing
FreeR set for adoption and extension, provided it is compatible. If
this is not given a new FreeR set is generated, but you should use
the same set for all isomorphous datasets in a project **(1.7)** .
Generation of the FreeR set takes into account potential merohedral
or pseudo-merohedral twinning, just in case, ie reflections which are
related by a potential twin operator are either all in the free set or
all in the working set.

An automatic resolution cutoff may be selected  **(1.8)** . Phaser
generates an estimate of the average information content as
bits/reflection as a function of resolution, and the cutoff is then
set to the point where the information content falls below the set
threshold, default 0.1 bits/reflection. If you use this automatic
cutoff, you should examine other indicators carefully to check that
the cutoff point is sensible (see resolution estimation later)

.. rubric:: Options
   :name: options

.. figure:: aimless_pipe_task_2.png
   :alt: Figure 2: Important options

   Figure 2: Important options

| This page (and the next one) set various options which are not
  usually required (and are not all documented here: for further
  information see the program documentation). You can override
  most of the default parameters used in Pointless and Aimless, for
  instance those for estimating the maximum effective resolution by
  various criteria **(2.1)** : note that these parameters do not
  control the actual resolution cutoff. Another possibly useful option is to
  change parameters for the optimisation of the estimated sigma(I)
  **(2.2)** (see Results page). The pull-down allows you to change the
  default option, which is to use the same correction" factors for all
  runs: in some circumstances, you may wish to have different paramaters 
  for each run, either separate ("individual") or restrained to be 
  "similar". The refinement may also be
  stabilised if necessary by fixing the sdB parameters.

.. rubric:: Additional options
   :name: additional-options

.. figure:: aimless_pipe_task_3.png
   :alt: Figure 3: Additional options

   Figure 3: Additional options

   The third pane offers further options, again ones you will not
   normally need to change.
   The first pull-down **(3.1)** allows you to override the default of
   IUCr standard settings for primitive 
   orthorhombic and centred monoclinic spacegroups: this convention
   sets a<b<c for all primitive orthorhombic, eg it allows P
   2\ :sub:`1` 2 2\ :sub:`1` rather than the "reference" setting P
   2\ :sub:`1` 2\ :sub:`1` 2, and also chooses the alternative I2 in
   place of C2 if it gives a smaller β angle. The default is
   recommended. The second pull-down **(3.2)** allows more detailed
   specification of the scaling model: options other than the default
   open additional panes (see eg **(4.4)** to set the interval on
   crystal rotation angle Phi (ROT) for the scaling). You can also change the
   criteria for outlier rejection **(3.3, 4.3)**, and define 'Runs'
   explicitly by batch range **(4.6)**, to override the automatic
   setting, and also set different high resolution cutoffs for each
   run.
   If the data has already been scaled by XDS or DIALS, you may wish
   to turn off scaling here ("No scaling only merge") and also turn off
   optimisation of the SD correction **(2.2)** , and just use this
   pipeline to get some statistics.

|
	
.. figure:: aimless_pipe_task_4.png
   :alt: Figure 4: Additional options expanded

   Figure 3: Additional options expanded

.. _Results:

Results
=======

The output from this task is quite detailed, and some of it is only
relevant if there are problems. You can drill down to the details of
the process, much of it presented as graphs. Based on the output, you
need to make some judgements about your data.

- What is the real resolution? Should you cut the high-resolution data?
- Are there bad batches (individual bad batches or ranges of batches)?
- Was the radiation damage such that you should exclude the later parts?
- Should you exclude some files or runs?
- Is there any apparent anomalous signal?
- Is the outlier detection working well?
- What is the overall quality of the dataset?
- How does it compare to other datasets for this project?

| Contents of the report description
| `1. Key Summary <#the-report-key-summary>`__
| `2. Overall Summary <#overall-summary>`__
| `3. Details of symmetry determination <./symmetry.html>`__
| `4. Details of scaling and merging <./scaling_and_merging.html>`__

.. rubric:: The report: key summary
   :name: the-report-key-summary

.. figure:: aimless_pipe_report_1.png
   :alt: Figure 1: Report headlines

   Figure 1: Report headlines

The report is divided into sections, and clicking in the top line
**(1.1)** will jump to the selected part of the report. The first
section **(1.2)** summarises the main conclusions from the task: (a)
choice of space group, and the confidence in this choice; (b)
various resolution estimates (c) key
statistics from merging, merging R-factors, I/sigma(I), CC(1/2) etc;
(d) any indications of translational NCS or twinning. Warnings are
colour-coded as orange or red (severe) - these may need more
investigation. There are also links to the log files for those who
like them **(1.3)** 

.. figure:: aimless_pipe_report_1a.png
   :alt: Figure 1a: Headline warnings

   Figure 1a: Report headline warnings

Warnings which should be noted include space group
ambiguity, from axial reflections missing or obscured by lattice
centring (eg **1a.4** ), and cases where alternative indexing
schemes are possible **(1a.5)** 


.. rubric:: Overall summary
   :name: overall-summary

.. container::

   .. figure:: aimless_pipe_report_2.png
      :alt: Figure 2: Main summary

      Figure 2: Main summary

   The Key summary is followed by the Overall summary tables and
   graphs. Panels **(2.1 and 2.2)** shows more detail of the space group
   determination. It may be worth looking a finer details later in
   the report, particularly if the confidence is low. Panel **(2.3)**
   is the classic "Table 1" suitable (more or less) for inclusion in
   a publication. The "Download" button writes this table as a CSV
   file which can be imported into word processing programs.
   Panel **(2.4)** shows the
   main graphs of statistics against resolution, which illustrate
   where a resolution cutoff might be applied, and also whether a
   significant anomalous signal is present. The example shown has
   strong data to the edge of collection (red and green lines of CC(1/2)
   values), and a significant anomalous signal at low
   resolution. Note that placing the mouse over the graphs
   brings up a brief description (in a yellow box). Other graphs can by
   accessed from the pull-down menu at the top, and a seperate graph
   window can be created using the arrow icon at the top right.
   Panel **(2.5)** shows
   graphs against "Batch number", ie image number within each Run or
   sweep of data. These graphs can be used to assess whether later
   parts of the data are suffer from serious radiation damage, and
   whether they may be omitted without losing too much completeness
   (see the Cumulative %completeness graph).
   Other useful plots against resolution presented here are the average
   "information content" (from Phaser) in bits/reflection **(2.6)** , and the
   Wilson plot **(2.7)**

   .. rubric:: Space group determination
      :name: space-group-determination

   | Remember that the space group is only a hypothesis until the
      structure is satisfactory solved and refined. 

   | `Details of symmetry determination are here <./symmetry.html>`__

   .. rubric:: Scaling and merging
      :name: scaling-and-merging

   | `Details of the results of scaling and merging are
      here <./scaling_and_merging.html>`__
