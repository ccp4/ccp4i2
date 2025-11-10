# CCP4i2 Changelog

## [2.5.0] - 2025-11-10

- New task to check a model against AU contents
- Standardised ModelCraft output file name
- Twinning analysis in Servalcat report
- Preserving ProSMART parameters when cloning Servalcat jobs

## [2.4.3] - 2025-10-24

- Improved checking/rebuilding of database and project-list XML files on startup

## [2.4.2] - 2025-10-22

- Capturing stderr for subjobs in Refmac, Servalcat and Lorestr pipelines
- Reducing size of reports by fetching log files instead of embedding them
- Fix for reflection CIF files containing both merged and unmerged data
- New Pointless options to remove lattice centering reflections
- Fix to Servalcat report type handling
- Fixes to i2run testing

## [2.4.1] - 2025-10-07

- Fix Coot map coluring with multiple models
- Fix for AceDRG atom name matching
- Improved graphs in Servalcat report
- Added checks before xmlnode append
- Searching Python 3.11 paths for CCP4 10
- PDB-REDO text change

## [2.4.0] - 2025-07-15

- Servalcat refinement against unmerged data
- Servalcat option for van der Waals restraint weight
- Updated Iris validation to work with v0.3.3
- More graphs in the Servalcat report
- More Servalcat i2 run tests
- Fix for deprecated numpy.float
- More AceDrg i2run tests
- Support for 5-letter ligand codes in AceDRG make link
- Fix for Coot 0.9 anomalous map colouring

## [2.3.3] - 2025-05-29

- Fix to Coot output files

## [2.3.2] - 2025-05-22

- Avoid re-importing files with the same checksum
- Fix for the Coot RSR Morph task on Windows
- Fixed collections.Iterable import
- Changed some performance testing thresholds

## [2.3.1] - 2025-05-17

- Changed some performance testing thresholds
- Removed unused CCP4I2Runner code
- Removed deprecated collections imports

## [2.3.0] - 2025-05-13

- Added a Coot 1 task
- Added a MetalCoord task
- Added performance testing to the i2run tests
- Fix to AUSPEX command line for Windows
- Fix for ASU contents view with QtGui/QtWidgets changes

## [2.2.6] - 2025-04-24

- Fix for import merged from a CIF file with a non-standard ASU

## [2.2.5] - 2025-04-23

- Fix for missing useLXML arguments in getEtree

## [2.2.4] - 2025-04-20

- Fixed window maximising on Windows
- ModelCraft used by default in the DR/MR/MB pipeline
- Warning in the Refmac pipeline when validation fails
- Updated neutron refinement options
- Refmac/Servcalcat pipelines only providing X-ray data to validation
- Observation type selection from CIF files in the import merged task
- Hiding nanobind leak warnings from Gemmi 0.7
- Changing Buccaneer to ModelCraft in the Phaser EP task
- Fixed MrParse report on Windows
- Increased maximum Coot files saved from 10 to 250
- Added an i2run script for Windows
- Added more i2run tests to replace test101

## [2.2.3] - 2025-02-07

- Fix Ramachandran plot
- Fix for broken ASU contents
- Fix embedded ModelCraft report

## [2.2.2] - 2025-01-31

- Fix to widgets loading a data file from the database
- Fix to export imported files and project files with a job selection
- Comparing models before and after refinement by Servalcat in Iris
- Updates for gemmi 0.7

## [2.2.1] - 2025-01-20

- Updated ModelCraft control parameters
- Refinement weighting text change
- Fix missing Dummy DbApi entry for ARP/wARP
- Disabling cryo-EM SPA refinement

## [2.2.0] - 2025-01-14

- New Servalcat refinement task
- Fixes to Buster arguments
- Fixed demo data path

## [2.1.1] - 2025-01-10

- Add log files, etc. to reports
- Remove reference to 'manual' model building
- Fix to report XML parsing
- Fix to MTZ cell comparison
- Fix running Lidia on Linux/Mac

## [2.1.0] - 2024-12-13
