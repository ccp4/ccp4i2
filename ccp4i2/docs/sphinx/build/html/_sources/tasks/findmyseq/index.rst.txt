######
FindMySequence
######

   FindMySequence is a neural-network-based approach for the rapid identification of unknown proteins in X-ray crystallography and cryo-EM. The implementation in i2 will find the most plausible sequence in a database (sequence search mode) for the user.

Input
=====

   |image1|

   The main input to the task is a reflection file, which may be any of
   the standard formats (mtz, mmcif). A corresponding map is also required, as well
   as a protein model. 

Output
======

   The (best) output sequence is recorded in the ccp4i2 database, and can easily be
   exported for external use (as in other interfaces, by right clicking on
   the data icon and selecting export).

Info
====

   Further information about FindMySequence can be found at the main web-site:-

   -  `FindMySequence Home Page <https://gitlab.com/gchojnowski/findmysequence>`__

.. |image1| image:: FMS_task_1.png


**References**

`G. Chojnowski, A.J. Simpkin, D.A. Leonardo, W. Seifert-Davila, D.E. Vivas-Ruiz, R.M. Keegan & D.J. Rigden. "FindMySequence: a neural-network-based approach for identification of unknown proteins in X-ray crystallography and cryo-EM". IUCr.J (2022). 9(I), 86-97 <https://journals.iucr.org/m/issues/2022/01/00/pw5018/>`_

