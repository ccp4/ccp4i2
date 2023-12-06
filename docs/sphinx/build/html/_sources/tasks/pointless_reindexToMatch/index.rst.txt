############################
Reindex or change spacegroup
############################

This task provides a number of utility functions for merged reflection
data files, primarily to change space group or reindex the data
relative to different axes, otherwise know as a change of basis. 
It will reindex a set of observations, and optionally a set
of Free-R flags, to match reference data as an external file, or
from explicit operators. There is also an option just to analyse the symmetry
of merged data, to check for under-merging or to check eg the cell conventions
in P1

Sometimes when comparing similar structures, the spacegroup or
indexing convention may be different even when the crystal packing is
similar. In this case the symmetry and the indexing convention changed
to match the existing structure.

Similarly during the course of a structure solution it may become
obvious that the spacegroup has been incorrectly determined. For
example, a symmetry which is present in the data may have been missed
leading to the selection of a spacegroup with lower symmetry, or a
higher spacegroup may have been selected due to the presence of
twinning or pseudo-symmetry. If point group has changed, then the data
processing 
should be repeated to obtain a correct list of unique reflections. If the space
group needs to be changed with the same point group (eg P 2 :sub:`1` 2 :sub:`1` 2 :sub:`1` to P 2 :sub:`1` 2 :sub:`1` 2)
then this utility can be used. Illegal changes will fail.


| 

.. figure:: pointless_reindexToMatch_1.png

    Figure 1: Input pane

|

Input
=====
	
One or two data files are processed, a reflection data object, and optionally
a matching Free-R set **(1.1)**

The pull-down **(1.2)** specifies what is to be done: options 1-4 need a
reference file  **(1.3)** , see below.

Options:

#. Observed data reference: an observed data object
#. Calculated data reference: calculated data
#. 2FoFc coefficient reference: map coefficients, eg from Refmac
#. Coordinate reference: a coordinate file from which intensities will be calculated
#. User specified spacegroup and reindex, see below
#. Remove lattice absences, see below
#. Just analyse data symmetry: possibly useful after an explicit reindex job, to check the result - no output file is produced

Explicit space group and reindex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: pointless_reindexToMatch_2.png

    Figure 2: User specified spacegroup and reindex

As stated in the figure 2 screenshot, you can give either a new space group
(in the same point group), or a reindex operator, or both. Pointless will try
to generate the other, but it doesn't always work in complicated cases.
If the reindexing operator leads to non-integral *hkl* indices, these will be
removed, eg [h/2, k, l] will halve the cell in the a direction (**a drastic
step, be careful that you know what you are doing!**): such Reindex operators
will change the space group - check that you think it correct, and if
necessary run the "Just analyse data symmetry" option on the output

Examples:

- P 3 :sub:`1` to P 3 :sub:`2`
- P 2 :sub:`1` 2 :sub:`1` 2 :sub:`1` to P 2 :sub:`1` 2 :sub:`1` 2
- Reindex P 2 2 :sub:`1` 2 :sub:`1` to P 2 :sub:`1` 2 :sub:`1` 2 using the reindex operator [k,l,h] (a cyclic  permutation to preserve the hand)
- C 2 to I 2
- [k, h, -l] in P 3 etc
- [h/2+k/2, -h/2+k/2, l] to remove reflections which would be C-centred: rotates the axes by 45°

Remove lattice absences
^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: pointless_reindexToMatch_3.png

    Figure 3: Removed centred lattice absences

.. WARNING:: THIS IS NOT RECOMMENDED UNLESS YOU KNOW WHAT YOU ARE DOING.

Occasionally a dataset has been wrongly integrated with a doubled cell, such
that eg half the spots are absent, or you would like to see what happens in a
half-cell when there is tNCS. Spots corresponding to centred lattices can be
removed with this option, eg h+k+l odd for an I lattice. This is easier than
working out the equivalent reindexing operator. The space group will change:
check that you think the new space group is correct

.. note:: an equivalent option is available in the Data reduction pipeline (under Additional options), with the same caveats. Maybe you should go back and look closely at your images
