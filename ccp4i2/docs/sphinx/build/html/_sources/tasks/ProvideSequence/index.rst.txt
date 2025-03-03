############################
Import sequence or alignment
############################

   The sequence import task allows you to import a sequence of sequence
   alignment object into the database. While it is possible to browse
   for a sequence in most task interfaces, it is more convenient to have
   the sequence already available in the database. One chain sequence is
   stored in one data object: complexes of multiple chains are then
   described by a combination of sequence objects.

Input
=====

   |image1|
   First, the type of object to be created must be selected. This may be
   either a sequence, or a sequence alignment. For experimental phasing
   projects it is generally sufficient to provide just the sequence of
   the molecule or complex of interest. For molecular replacement
   projects, it may be useful to import a sequence alignment with
   homologous structures, if one has been determined from an external
   source **(1)**

   Next, the sequence information is entered. This may be in the form of
   an existing file (or even a database object), in which case the file
   or object are selected using one of the options in section **(2)**.

   Otherwise, a sequence or sequence alignment may be pasted into the
   text box **(3)**, using either fasta (.fasta) format for sequences or
   CLUSTALW format for alignments (.aln). Other sequence formats may be
   read correctly, however the results should be checked.

Results
=======

      |image2|
      The sequence entered and input and output files are reported.

.. |image1| image:: ProvideSequence_task_1.png
.. |image2| image:: ProvideSequence_report_1.png
