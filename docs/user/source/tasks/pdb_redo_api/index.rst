################################
PDB-REDO web service
################################

.. warning::

   This task **sends your model and your reflection data to the PDB-REDO
   server** (pdb-redo.eu, in the Netherlands) and waits for the result. Use it
   only for data you are allowed to share with an outside service. Nothing
   leaves your machine until you press Run.

PDB-REDO re-refines a model with REFMAC5 and then rebuilds and optimises it
(rotamers, peptide flips, loops, waters, TLS groups and so on), on the
PDB-REDO group's servers rather than on your own computer. Use it when you
want to see what an independent, automated re-refinement and rebuild does
to a model: before deposition, or to compare it with the model you have
refined yourself. The refinement tasks in the Refinement category
(:doc:`../prosmart_refmac/index`, :doc:`../servalcat_pipe/index`) run
locally and keep your data on your machine; this one does not. The
:doc:`../dnatco_pipe/index` page compares a deposited model with its PDB-REDO
re-refinement.

This page was written from the task's code and an unrun job. The task was
not run for the pictures, since running it would send the data to the
service, so what is said about the results is taken from what the task does
with them, not from a result seen.

Before you run: the token
=========================

The service needs an API token, which is free from
`pdb-redo.eu/token <https://pdb-redo.eu/token>`_. A token is a pair, a
*token ID* and a *token secret*. It is not a task parameter: it is not kept
in the job, the project database or an exported project. Set it from the
banner at the top of the task's interface (a *Set token...* button, with a
*Test* and a *Replace* button once one is set), or under
Preferences, Credentials. It can also be given by the environment variables
``PDB_REDO_TOKEN_ID`` and ``PDB_REDO_TOKEN_SECRET``, which win over a stored
one. The secret is used to sign requests on your machine.

Without a token the job fails validation, with an error naming the missing
token. When you press Run, the task also asks the service whether the
token is accepted, so a revoked or mistyped token is found before anything
is uploaded.

Input
=====

.. figure:: pdb_redo_input.png
   :alt: Figure 1: PDB-REDO input

   Figure 1: PDB-REDO input

The model to be refined **(1)**, the observed data **(2)** and the free R
set that goes with them **(3)** are required. The observed data and free R
set are combined into one MTZ file for upload. The sequence **(4)** and a
geometry dictionary for any ligand **(5)** are optional and are uploaded
too if given. As with any refinement, use the free R set the model was
last refined against.

**Paired refinement (6)** asks the service to do a paired refinement, to
help judge the resolution limit (compare :doc:`../pairef/index`).

The **Rebuilding Options** switch parts of PDB-REDO's rebuilding off: loop
completion, peptide flips, side-chain rebuilding or completion, deleting
poor waters, rebuilding carbohydrates, or all rebuilding.

.. figure:: pdb_redo_advanced.png
   :alt: Figure 2: PDB-REDO advanced options

   Figure 2: PDB-REDO advanced options

The Advanced Options tab: update the model even if the initial refinement
is poor; force isotropic or anisotropic B-factors; do not use TLS;
and tighter or looser restraints, each 0, 1 or 2 (the two values are
subtracted, so setting both to 2 has no net effect). Each option is sent
to the service as a flag; what the service does with it is described in
PDB-REDO's own documentation.

The run
=======

The task uploads the files, then asks the service about the run every few
seconds, less often as the wait grows (every 60 seconds at most), and
reports "still running" in the log every ten minutes. A run can take a long
time, and the task sets no limit on the wait for a run that is going well.
Brief network failures are retried; after twelve failures in a row (about
ten minutes) the task gives up waiting, and says which PDB-REDO run number
to look for at pdb-redo.eu, since the run itself is probably still going.
If the token stops being accepted, or the service reports the run stopped,
the job fails with the run number; after a stop the task tries to bring
back PDB-REDO's logs for diagnosis.

Results
=======

When the run ends, the task downloads the results and unpacks from them:

- the fully optimised model and its map coefficients (density and
  difference map);
- the "conservatively optimised" model (the best-TLS result) and its maps;
- the service's logs, including the final REFMAC5 log, shown in folds in
  the report;
- R and R-free from the run's own summary, recorded as the job's
  performance.

The report links to PDB-REDO's results page, and says the service keeps
the run's own report on its website for 21 days. The results page loads
its content from pdb-redo.eu, so opening it contacts the service again.

Which model to take is for you to judge: the fully optimised model has
had the most rebuilding; the conservatively optimised one the least change.
Compare either with your own model by R-free and, more importantly, by
what changed (rotamers, flipped peptides, deleted waters) and by looking
at the maps in Moorhen before accepting any change. An automated rebuild
is a suggestion, not a verdict.

**Reference**

Joosten, R. P., Long, F., Murshudov, G. N. & Perrakis, A. (2014). The
PDB_REDO server for macromolecular structure model optimization. IUCrJ 1,
213-220. https://doi.org/10.1107/S2052252514009324
