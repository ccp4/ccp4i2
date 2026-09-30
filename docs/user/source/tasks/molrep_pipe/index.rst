##############################
Molecular replacement - MOLREP
##############################

Input
=====

   MOLREP is a quick and simple-to-use molecular replacement Task using
   a homologous search structure, (and optionally) the sequence of the
   target model. The best solution is subsequently submitted to REFMAC
   for 5 cycles of refinement.

   Full documentation can be found at:
   http://www.ccp4.ac.uk/html/molrep.html

   The Task only searches for a single or multiple copies of one
   component in each run. The expected number of copies can be input by
   the user, or estimated by the program, assuming a solvent content of
   about 50%

   For multi-component structures, MOLREP will attempt to place the full
   complement of the first component(s), and then a second MOLREP task
   to search for the next component should be carried out with the first
   solution model defined as a 'fixed' structure.

   If a sequence is provided, MOLREP will perform a sequence alignment
   of target sequence to that of a SINGLE CHAIN of the search model, and
   the input chain is then modified by renaming and pruning
   non-conserved residues. Therefore, if you wish to use a search model
   that is a multimer with several copies of a single sequence, it is
   essential to use CHAINSAW or SCULPTOR to edit all chains of the
   model.

Essential input items
---------------------

   |image1|
   
   Reflection data **(1)**, Free R set **(2)**, Input model **(3)**.
   Atom Selection **(4)**: it is possible to select only part of the
   search model coordinates; e.g. only peptide, or a single chain, or a
   span of residues.

   Prior to running MOLREP, sophisticated model preparation can use the
   Bioinformatics tools provided. The Interactive Model preparation
   using MrBump & CCP4mg is comprehensive and produces an ensemble of
   aligned pruned models. Alignment of the template and target sequence
   can be obtained from tools such as HHPred, PROMALS3D or FFAS; and
   these alignments used in the SCULPTOR or CHAINSAW Tasks to edit the
   search model.

Optional input items
--------------------

   If no sequence is input, the search model is used as is.

   Sequence of Target Model **(5)**. The sequence of one component only
   can be provided either from the "Define AU contents" task, selecting
   the single sequence to be used, or from a Sequence file giving the
   sequence of a single chain. MOLREP will then perform a sequence
   alignment of the first chain of the search model and target sequence,
   and modify the input model by renaming and pruning non-conserved
   residues. The names of the retained atoms will be changed to match
   the target, while conserved residues are left unchanged. The residues
   in the modified model are numbered to be consistent with the target
   sequence, i.e. if a residue in the target corresponds to a gap in the
   model, the residue numbers in the output model will contain a gap,
   but if a residue in the model corresponds to a gap in the target, the
   residue numbers in the output will be consecutive.

   If the *AU contents* contains multiple sequences, then a single
   sequence must be selected for each MOLREP pass.

   Number of monomers to search for **(6)**. The number of copies of the
   given model expected can be specified by the user, or automatically
   estimated by MOLREP. MOLREP will take the best results for a first
   search, fix this result, and search for second, third, and so on,
   copies.

   Fixed Model: **(7)** If some of the structure is already known from a
   previous MR run those coordinates should be entered as the 'fixed'
   structure.

Basic Options
-------------
 
   |image2|

   Space group Selection. **(1)**. Where there is ambiguity, MOLREP will
   attempt the first component search in more than one space group.
   After identifying the most likely space group for the first copy,
   subsequent searches for multiple copies will be carried out in that
   space group..

   Resolution range **(2)**. MOLREP uses data to the specified
   resolution limit. The refinement step uses all input data.

Results
=======

      |image3|


.. |image1| image:: molrep_pipe_task_1.png
.. |image2| image:: molrep_pipe_task_2.png
.. |image3| image:: molrep_pipe_report_1.png
