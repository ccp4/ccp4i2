CCP4i2 Tutorial: GerE (PDB:1FSE)
================================

This structure can be solved quickly by either experimental phasing or
molecular replacement and is a good one to try as an introduction to
CCP4.

The data are from the crystal structure of GerE, a transcription
activator from Bacillus subtilis, which was solved by MAD phasing using
the Se signal (V.M.A. Ducros, R.J. Lewis, C.S. Verma, E.J. Dodson, G.
Leonard, J.P. Turkenburg, G.N. Murshudov, A.J. Wilkinson and J.A.
Brannigan, J. Mol. Biol. (2001) 306 759-771).

With modern software, the structure can also be easily solved by SAD
phasing from just the Se peak data using any of the software available.
The structure has 6-fold non-crystallographic symmetry (NCS) with 3
distinct dimers in the assymetric unit, however there are only 2 Se
atoms per molecule, so the NCS is difficult to determine from the heavy
atom coordinates alone. A typical strategy involves density modification
without NCS, then automated model building to obtain a preliminary
model. This model may be used to provide NCS information for a better
density modification calculation.

SIRAS phasing may also be attempted using the native data.

Data files
----------

The data files required for this tutorial are in the demo_data/gere folder. You should
see a demo_data entry in the File Browser side bar. If not you can use 
Utilities -> Copy demo data to project -> gere to copy the files into your project.

-  **gere.fasta**: GerE sequence
-  **gere.mtz**: MTZ file containing the anomolous *intensity*
   data measured at the Se peak, and a Free-R flag.
-  **gere_nat.mtz**: MTZ file containing native data.
   While anomalous data are provided, there is no significant anomalous
   signal.
-  **gere_scaled_data.mtz**: MTZ file containing
   scaled and merged *amplitude* data for all the MAD wavelengths.
-  **gere_molA.pdb**: Structure for a single molecule.

Last modified: Fri Feb 20 17:26:38 GMT 2015
