CCP4i2 Tutorial: RNAse (PDB:1SAR)
=================================

This structure can be solved (1) by experimental phasing using SAD or
SIRAS, or (2) by molecular replacement. In addition, there are data for
two complexes allowing automated or manual ligand fitting.

Reference for the structure of ribonuclease from *Streptomyces
aureofaciens* (RNase Sa) and its complex with guanosine 3′monophosphate,
3′GMP: *Sevcik, J., Dodson, E. J., & Dodson, G. G. (1991). Determination
and restrained least-squares refinement of the structures of
ribonuclease Sa and its complex with 3′guanylic acid at 1.8 Å
resolution. Acta Crystallographica Section B: Structural Science, 47(2),
240-253.*

Experimental phasing
--------------------

The structure was first solved using heavy atom derivatives. Data for
three derivatives (Hg, I and Pt), including the anomalous components
(collected using CuK-α radiation) can be loaded from the first mtz file
(rnase25_Nat_Derivs.mtz), which also contains the native enzyme data,
all at 2.5 Å resolution. This MTZ file contains structure factor
amplitudes (not intensities) for both data sets.

The native data are in the columns labelled:

Column FNAT SIGFNAT

And the derivative data are in the columns labelled:

Column FHG(+) SIGFHG(+) FHG(-) SIGFHG(-)

Column FI(+) SIGFI(+) FI(-) SIGFI(-)

Column FPTNCD25(+) SIGFPTNCD25(+) FPTNCD25(-) SIGFPTNCD25(-)

The structure can be solved using the CRANK2 or SHELXCDE pipelines,
using SAD or SIRAS. The Pt is the best derivative, and there are
estimated to be around 5 Pt atoms per monomer. The Hg derivative does
now lead to a simple solution - the structures was originally solved by
MIR.

Molecular Replacement
---------------------

PDB files for two models are provided. Both lead to a straightforward
solution allowing rebuilding with Buccaneer. Other models can be found
using the MrBUMP-CCP4mg Task, starting from the sequence.

Ligand Fitting
--------------

Amplitudes for the ligand complexes and a higher resolution native set
are in the 1.8 Å data file. These can be used for manual or automated
Ligand Fitting. The names, PDB codes for the complexes, 3 letter ligand
code and canonical smiles strings for the ligands are:

3′GMP PDB: 2SAR 3GP

c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)O)O)N=C(NC2=O)N

2′GMP PDB: 3DGY 2GP

c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)OP(=O)(O)O)N=C(NC2=O)N

Data files
----------

The data files required for this tutorial are in the demo_data/rnase folder. You should
see a demo_data entry in the File Browser side bar. If not you can use 
Utilities -> Copy demo data to project -> rnase to copy the files into your project.


-  **rnase.fasta**: RNAse sequence
-  **rnase18_Nat_Complexes.mtz**: MTZ file
   containing native amplitudes to 1.8 Å plus two complexes, the first
   with 3′GMP (3GP) and the second with 2′GMP (2GP).
-  **rnase25_Nat_Derivs.mtz**: MTZ file
   containing native amplitudes to 2.5 Å plus three heavy atom
   derivatives with anomalous data, Hg, PT and I.
-  **rnase_model.pdb**: Refined structure.
-  **1sar.pdb**: Refined apo-native 1.8 Å structure.
-  **3a5e.pdb**: MR model with 95% sequence identify.
-  **1mgw.pdb**: MR model with 69% sequence identify.
-  **3GMP_smiles.txt**: smiles string for 3′GMP
-  **2GMP_smiles.txt**: smiles string for 2′GMP

Last modified: Fri Feb 20 17:26:38 GMT 2015
