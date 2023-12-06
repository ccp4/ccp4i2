############################################################
Convert between Hendrickson-Lattman coefficients and phi/FOM
############################################################

   This task will convert phase probabilities between the two common
   forms: Hendrickson-Lattman coefficients, and phase and
   figure-of-merit. *Note that you should never do this, except to
   generate phases for export. CCP4i2 will convert phases internally as
   and when required. Converting to phase and figure-of-merit loses
   information.*

Input
=====

   |image1|
   
   The input to the task is a set of phases. If the phases are in the
   form of Hendrickson-Lattman coefficients, they will be converted to
   phase and figure-of-merit. If they are in the form of phase and
   figure-of-merit, they will be converted to Hendrickson-Lattman
   coefficients.

.. |image1| image:: chltofom_task_1.png
