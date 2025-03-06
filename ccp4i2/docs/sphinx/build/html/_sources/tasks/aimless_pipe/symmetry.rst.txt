 .. include:: <isonum.txt>

##################################################
Symmetry determination details (program POINTLESS)
##################################################

The indexing of the lattice in an integration program (such as Mosflm,
or the indexing and integration stages of XDS and DIALS)
is based on the lattice geometry, with no regard for the symmetry of the
diffraction pattern, which can only be determined after integration. The
default option here scores potential symmetry operators in the diffraction
pattern, and ranks the possible rotational groups (Laue group,
Patterson group or point group, see Footnote), and also
inspects axial reflections for possible systematic absences which may
indicate a likely space group. Be careful, the presence of
pseudo-symmetry may suggest a higher symmetry than the truth.
POINTLESS tries to allow for this possibility, but inspection of the scores for
individual symmetry elements may help to indicate the correct space
group in difficult cases (fig 3). 
Conversely, severe radiation damage may obscure the true symmetry.
Weak high resolution data does not contain
reliable information for the determination of rotational symmetry (the
Laue or point group), so POINTLESS may cut the resolution just for this
purpose. If possible it cuts the data at the point where CC(|frac12|) falls
below 0.6 (calculated using only the Friedel symmetry (ie in P1)), or
I/σ(I)>6 if there are insufficent data to calculate CC(|frac12|).

This first stage of the task has 4 options, with different result
reports:

#. `Determine Laue group and space group`_
#. `Match index to reference data`_
#. `Choose a known or previous solution`_
#. `No symmetry testing just copy and combine input files`_

`Details of results of scaling and merging are shown on another
page <./scaling_and_merging.html>`__

`Back to main page <./aimless_pipe.html>`__

.. _Determine Laue group and space group:

1. Automatic Determination of Laue group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: aimless_pipe_report_3.png
   :alt: Figure 3: symmetry element scores

   Figure 3: symmetry element scores

*Scoring individual symmetry elements.* Each possible symmetry operator
in the lattice is scored separately, by a pairwise correlation
coefficient (CC) between E\ :sup:`2` for observations related by that
operator, and also by an R-factor, R\ :sub:`meas`. The CC is used to
estimate a probability, allowing for possible small samples by comparing
it with the distribution CCs of equal sized samples of unrelated
observation pairs. This separate scoring is useful in cases of
pseudo-symmetry to indicate the true symmetry, in cases where the
program gets it wrong. In the example in figure 3, the true symmetry is
C222, but because of the accidental combination of cell lengths, b ≈ √3
a, the lattice can be indexed as hexagonal. The scores show that three
orthogonal dyads are present, but the other potential operators of the
hexagonal lattice are absent. 

.. figure:: aimless_pipe_report_4.png
   :alt: Figure 4: Laue group scores

   Figure 4: Laue group scores

*Scoring Laue groups.* The scores from the individual elements are then
combined into a joint score for all possible Laue groups which are
subgroups of the lattice group. A high score for a symmetry element
which is present in the lattice group but not in the Laue group will
count against that group (the CC-, Zc- and R- columns in the table).
Note that in Figure 4 the pseudo-hexagonal lattice can accommodate three
possible Cmmm settings, 60° apart: the one chosen randomly
by the original indexing (Reindex operator [h,k,l]) was wrong here.

.. figure:: aimless_pipe_report_5.png
   :alt: Figure 5: Systematic absences

   Figure 5: Systematic absences

*Systematic absence scores.* Within a chosen Laue group, space groups
(all chiral groups apart from the pairs I222 and
I2\ :sub:`1`\ 2\ :sub:`1`\ 2\ :sub:`1`, and I23 and I2\ :sub:`1`\ 3) may
be distinguished by the presence or absence of screw axes along the cell
edges. A screw axis leads to systematic absences along an axis of the
reciprocal lattice eg a 2\ :sub:`1` screw along b in space group
P2\ :sub:`1` makes axial reflections 0k0 present only when k is even.
Detecting systematic absences may be unreliable because the axial
reflections may be few in number (as along the a axis in figure 5), or
missing from the dataset if they lie along the spindle rotation axis in
the data collection (in the blind region), or may be misleading if there
are approximate non-crystallographic screw axes (tNCS), but in many cases they
can suggest the space group, to be confirmed later: the space group
remains a hypothesis until the structure is satisfactorily solved. In
POINTLESS, a Fourier analysis is used to detect periodicity, on I'/σ:
the intensity I' used here is adjusted by subtraction of a small
fraction (default 0.02) of the intensity of the neighbouring reflection
along the axis, to allow for possible overlap of a nearby strong
reflection (figure 5).

.. figure:: aimless_pipe_report_6.png
   :alt: Figure 6: Space group choice

   Figure 6: Space group choice

| *Choice of space group or point group.* Possible space groups are
  ranked according to their total probability = Laue group probability |times|
  Systematic absence probability. If there is a unique solution with the
  highest total probability, this will be chosen as a "space group"
  solution. If some of the potential systematic absence data is missing,
  then more than one space group has the same score, and the space
  group without translations
  will be chosen, eg P 2, P 2 2 2, P 4 2 2, P 6 2 2 etc.
  Enantiomorphic space groups will
  have the same score, so the first one is chosen as a space group
  solution (see figure 6), unless the other one matches that in the
  input file. As well as the scores, a "confidence" score
  is calculated for both the Laue group and the space group, defined as
  √[Score |times| (Score - NextBestScore)]: 
  these values are printed in a
  summary table as the Best Solution. If you prefer a different solution
  to that chosen by the program, you can rerun the program with the
  Choose a previous solution option.

.. _Match index to reference data:

2. Match index to reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: aimless_pipe_report_7.png
   :alt: Figure 7: alternative indexing scores

   Figure 7: alternative indexing scores

In some point-groups there are more than one (typically two or four)
valid but non-equivalent indexing possibilities. For your first crystal,
you may choose any of these, but subsequent crystals must match the
first. This problem generally arises in cases where the point group
symmetry is less than the lattice symmetry: the alternative indexing
schemes are related by the symmetry operators present in the lattice but
not in the point group. These are the same point groups which may lead
to merohedral twinning, eg point group P3 has four possible indexing
schemes in the lattice point-group P622 (`further
explanation <https://www.ccp4.ac.uk/html/reindexing.html>`__). Within
multiple datasets from the same crystal (eg MAD), you can avoid this
problem by only autoindexing one set, & using the same indexing matrix
for the others, but different crystals must be explicitly checked for
consistent indexing. Reindexing ambiguities may also arise in
lower-symmetry point-groups in case of accidental coincidences or
relationships between cell dimensions (eg a ≈ b in orthorhombic or the
pseudo-hexagonal C222 example in figure 3). 

| A ranked table of scores for each possible indexing scheme is printed.
  The example in figure 7 is unusually complicated: the true space group
  is rhombohedral, R3, but the lattice is pseudo-cubic; in the
  rhombohedral axis system (with a = b = c, α = β = γ) the angles are
  close to 90°. This means that in addition to the usual ambiguity in R3
  with alternative indexing [k,h,-l] (in the hexagonal setting aka H3
  or R3:H) ,
  there are four possible directions for true 3-fold axis, along the
  body-diagonals of the cubic lattice, leading to a total of eight
  possible indexing schemes.

.. _Choose a known or previous solution:

3. Choose known or previous solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| This option allows an explicit choice of indexing or symmetry. The
  Choose solution from search options perform the usual search, but
  selects the given solution even if it does not have the highest score.
  The Specify Laue group name or Specify space group name use the given
  groups without doing the searches, ie the program just changes the
  group, with an option to reindex.

.. _No symmetry testing just copy and combine input files:

4. No symmetry testing just copy and combine input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| This option just sorts one or more input files, usually because there
  are more than one, and you know the space group.

Footnote
^^^^^^^^

Note that strictly a Laue group is the point group
symmetry of the diffraction pattern, ie the crystal point group plus a
centre of inversion at the origin of the reciprocal lattice: here we
include the lattice centring type (P, I, F, C, R) in the "Laue group"
definition, equivalent to the Patterson space group.
