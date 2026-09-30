import_serial
=============

This task imports merged data from serial macromolecular crystallography in CCP4. Data from CrystFEL (HKL format) are converted to MTZ format, including required specification of space group, unit cell parameters and wavelength. The information about symmetry can be given manually or extracted from a reference files (CrystFEL stream file, CrystFEL cell_explorer cell file, reference PDB or mmCIF or MTZ file). Data quality statistics are calculated and reported. To get values of CC1/2 and Rsplit, two half-data-set files have to be provided. Import of data from xia2.ssx is also supported, in this case only a single MTZ file is required. High resolution cutoff is recommended to adjust.

Example of usage
----------------

CrystFEL:

* Merged data file (I): data.hkl
* Merged half data set (I): data.hkl1
* Merged half data set (I): data.hkl2
* CrystFEL cell file: cell.cell
* Space group: P 62 2 2
* Wavelength (A): 1.13
* High resolution cutoff (A): 1.6

xia2.ssx:

* Merged data file (I): merged.mtz
* High resolution cutoff (A): 1.6

It is also possible to run this program in command line or CCP4Console. The obtained MTZ files are recommended to import in CCP4 using the import merged task.

.. code ::

   $ ccp4-python -m import_serial --hklin data.hkl --half-dataset data.hkl1 data.hkl2 --cellfile cell.cell --spacegroup P6222 --wavelength 1.13 --dmin 1.6
   $ ccp4-python -m import_serial --hklin merged.mtz --dmin 1.6

Test and example data are available in at https://github.com/MartinMalyMM/import_serial_test_data

List of all options in command line
-----------------------------------

.. code ::

   $ ccp4-python -m import_serial --help
   
   usage: import_serial [-h] --hklin HKLIN [--half-dataset HKL1 HKL2] [--wavelength WAVELENGTH] 
                        [--spacegroup SPACEGROUP] [--cell a b c alpha beta gamma] [--cellfile CELLFILE]
                        [--streamfile STREAMFILE] [--reference REFERENCE] [--dmin D_MIN] [--dmax D_MAX]
                        [--nbins N_BINS] [--project PROJECT] [--crystal CRYST] [--dataset DATASET] 
   
   Calculate statistics of serial MX data from xia2.ssx or CrystFEL and import them to CCP4
   
   optional arguments:
     -h, --help            show this help message and exit
     --hklin HKLIN, --HKLIN HKLIN
                           Specify merged mtz file from xia2.ssx or merged hkl file from CrystFEL
     --half-dataset HKL1 HKL2
                           CrystFEL only: two half-data-set merge files (usually .hkl1 and .hkl2)
     --wavelength WAVELENGTH, -w WAVELENGTH
                           Wavelength (only for data from CrystFEL)
     --spacegroup SPACEGROUP
                           Space group
     --cell a b c alpha beta gamma
                           Unit cell parameters divided by spaces, e.g. 60 50 40 90 90 90
     --cellfile CELLFILE   Cell file from CrystFEL
     --streamfile STREAMFILE
                           Stream file from CrystFEL
     --reference REFERENCE, --ref REFERENCE, --pdb REFERENCE, --cif REFERENCE, --mmcif REFERENCE
                           Reference file (PDB, mmCIF or MTZ) to provide spacegroup and unit cell
     --dmin D_MIN, --highres D_MIN
                           High-resolution cutoff
     --dmax D_MAX, --lowres D_MAX
                           Low-resolution cutoff
     --nbins N_BINS, --nshells N_BINS
                           Number of resolution bins
     --project PROJECT     Project name
     --crystal CRYST       Crystal name
     --dataset DATASET     Dataset name

This program has been developed by Martin Malý, University of Southampton, `martin.maly@soton.ac.uk <mailto:martin.maly@soton.ac.uk>`_
