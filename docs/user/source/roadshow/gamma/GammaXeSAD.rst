Programs
--------

| The CRANK pipeline:
| (ShelxC, ShelxD, SFTools, Ecalc, REFMAC5, Peakmax, Parrot, Buccaneer)

Scenario
--------

Adapter proteins of the AP2 complex select membrane protein cargo for
capture into vesicles. The gamma adapter binds to peptide motifs on,
amongst others, the EGF receptor, resulting in receptor internalization
and recycling. Recognition of the cytoplasmic tail of EGFR is mediated
by the ear domain, which was crystallized and into which xenon atoms
have been introduced by pressurization of the crystals prior to crash
cooling. A SAD dataset was collected on a home source from these
crystals.

Data location
-------------

Data for this tutorial can be downloaded into a new CCP4i2 project using
the Utilities menu: Utilities->Download test data‑>gamma. Files are then
available in a subdirectory of your project’s top directory

Task 1: Run CRANK to solve the structure
----------------------------------------

|Crank pipeline|

Task location
-------------

All you need to provide is the SAD data, a sequence file, and the
identity of the anomalous scatterer.

|Crank pipeline location|

Task 2: Run the CRANK pipeline part way ...
-------------------------------------------

. |Run the Crank pipeline partway|

Task 3: ... and pick up with other tools
----------------------------------------

|Other tools - density modification, model building|

The autobuild task crashes because the FreeR set that comes out of the
CRANK pipeline is not quite complete. To generate a truly complete Free
R set, you will need to use the CCP4i2 FreeR set (or load one in that
wou have created another way !)

|Other tools Free-R|

.. |Crank pipeline| image:: Figure1.png
   :height: 450px
.. |Crank pipeline location| image:: Figure2.png
   :height: 100px
.. |Run the Crank pipeline partway| image:: Figure3.png
   :height: 100px
.. |Other tools - density modification, model building| image:: Figure4.png
   :height: 100px
.. |Other tools Free-R| image:: Figure5.png
   :height: 100px
