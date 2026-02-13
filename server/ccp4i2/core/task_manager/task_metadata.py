# Task metadata extracted from CCP4 CachedLookups.json.
# Contains TASKTITLE, DESCRIPTION, shortTitle for UI display.

# CCP4i2 has a GUI/Plugin split:
# CTaskWidget (GUI) classes define TASKMODULE, DESCRIPTION, and reference a CPluginScript via TASKNAME.
# Multiple GUIs can reference the same plugin (one-to-many).
# This file currently captures a flattened view (one entry per plugin).

# FUTURE: Parse GUI classes (*_gui.py, CTask*.py) without importing Qt.
# Extract TASKNAME to link GUI→Plugin.
# Support multiple GUI entries per plugin for different use cases/folders.

TASK_METADATA = {
    "AMPLE": {
        "TASKTITLE": "Molecular Replacement with unconventional models - AMPLE",
        "DESCRIPTION": "This task is for running Molecular Replacement with unconventional models",
        "shortTitle": "AMPLE Molecular Replacement Pipeline",
    },
    "AUSPEX": {
        "TASKTITLE": "Graphical diagnostics by AUSPEX plots",
        "DESCRIPTION": "Use AUSPEX, generate graphical diagnostics for data set",
        "shortTitle": "AUSPEX",
    },
    "AcedrgLink": {
        "TASKTITLE": "AceDRG in link generation mode",
        "DESCRIPTION": "AceDRG in link generation mode",
        "shortTitle": "AceDRG in link generation mode",
    },
    "AlternativeImportXIA2": {
        "TASKTITLE": "Import Xia2 results",
        "DESCRIPTION": "Harvest merged and unmerged files from selected XIA2 data reduction protocols",
        "shortTitle": "Import Xia2 results",
    },
    "LidiaAcedrg": {
        "TASKTITLE": "MAKE LIGAND",
        "DESCRIPTION": "Generate a PDB file and dictionary (acedrg) from MOL file, SMILES, or sketch (lidia)",
        "shortTitle": "MAKE LIGAND",
        "RANK": 1,
    },
    "LidiaAcedrgNew": {
        "TASKTITLE": "Make Ligand - AceDRG",
        "DESCRIPTION": "Generate a PDB file and dictionary (acedrg) from MOL file, SMILES string/file, or sketch (lidia).<br>Optionally match atom names to known structures.",
        "shortTitle": "Make Ligand - AceDRG",
        "RANK": 1,
    },
    "MakeLink": {
        "TASKTITLE": "Make Covalent Link - AceDRG",
        "DESCRIPTION": "Generate a link dictionary to describe a covalent bond between two monomers, allowing modification of the linked monomers",
        "shortTitle": "Make Covalent Link - AceDRG",
    },
    "MakeMonster": {
        "TASKTITLE": "Export monster mtz",
        "DESCRIPTION": "Build an MTZ file from data objects",
        "shortTitle": "Export monster mtz",
    },
    "MakeProjectsAndDoLigandPipeline": {
        "TASKTITLE": "Make projects and do ligand pipeline",
        "DESCRIPTION": "A task to generate projects and apply ligand pipeline for multiple datasets",
        "shortTitle": "Make projects and do ligand pipeline",
    },
    "PrepareDeposit": {
        "TASKTITLE": "Prepare files for deposition",
        "DESCRIPTION": "Prepare files for deposition",
        "shortTitle": "Prepare files for deposition",
    },
    "ProvideAlignment": {
        "TASKTITLE": "Import an alignment",
        "DESCRIPTION": "Enter an alignment from a file or by cut and paste",
        "shortTitle": "Import an alignment",
    },
    "ProvideAsuContents": {
        "TASKTITLE": "Define AU contents",
        "DESCRIPTION": "Define AU contents, estimate molecular weight, Matthews probability from list of sequences",
        "shortTitle": "Define AU contents",
    },
    "ProvideSequence": {
        "TASKTITLE": "Import sequence(s)",
        "DESCRIPTION": "Enter one or more sequences from a sequence file, from a PDB, or by cut and paste",
        "shortTitle": "Import sequence(s)",
    },
    "ProvideTLS": {
        "TASKTITLE": "Import and/or edit TLS set definitions",
        "DESCRIPTION": "Enter TLS information to be used later in the project",
        "shortTitle": "Import TLS",
    },
    "SIMBAD": {
        "TASKTITLE": "Sequence Free Molecular Replacement - SIMBAD",
        "DESCRIPTION": "This task is for running Molecular Replacement without a sequence",
        "shortTitle": "SIMBAD Molecular Replacement Pipeline",
    },
    "ShelxCD": {
        "TASKTITLE": "Find HA sites - SHELXC/D",
        "DESCRIPTION": "Find sites from SAD/MAD/SIR/SIRAS/RIP/RIPAS data",
        "shortTitle": "ShelxCD",
    },
    "SubstituteLigand": {
        "TASKTITLE": "Automated solution of isomorphous ligand complex",
        "DESCRIPTION": "A ligand workflow, starting from merged or unmerged reflections, SMILES, and an isomorphous parent structure",
        "shortTitle": "Isomorphous ligand solution",
        "RANK": 1,
    },
    "SubtractNative": {
        "TASKTITLE": "Subtract a defied fraction of an atom map from an input map",
        "DESCRIPTION": "Reads a set of map coefficients, generates a corresponding map, and modifies it by subtracting some fraction of the FC map corresponding to a set of coordinates that have been provided. The output is a set of modified map coefficients",
        "shortTitle": "Subtract native map",
    },
    "TestObsConversions": {
        "TASKTITLE": "Test CCP4i2 observed data interconversions",
        "DESCRIPTION": "Exercise CCP4i2 Observed data representation conversions",
        "shortTitle": "Test CCP4i2 observed data interconversions",
    },
    "acorn": {
        "TASKTITLE": "ACORN - Phase Refinement with Dynamic Density Modification",
        "DESCRIPTION": "Un-biased improvement of initial phases for high resolution data (1.5 Angstoms and better)",
        "shortTitle": "ACORN",
    },
    "add_fractional_coords": {
        "TASKTITLE": "Add Fractional Coordinates",
        "DESCRIPTION": "Calculate fractional coordinates and output them in mmCIF format",
        "shortTitle": "Add Fractional Coordinates",
    },
    "adding_stats_to_mmcif_i2": {
        "TASKTITLE": "Prepare and validate files for deposition",
        "DESCRIPTION": "Add data reduction statistics and sequence information into coordinates for deposition",
        "shortTitle": "New deposition task",
    },
    "aimless_pipe": {
        "TASKTITLE": "Data reduction - AIMLESS",
        "DESCRIPTION": "Scale and analyse unmerged data and suggest space group (Pointless, Aimless, Ctruncate, FreeRflag)",
        "shortTitle": "Data reduction",
        "RANK": 1,
    },
    "arcimboldo": {
        "TASKTITLE": "Ab initio phasing and chain tracing - ARCIMBOLDO (LITE, BORGES, SHREDDER)",
        "DESCRIPTION": "Structure solution from ideal molecular fragments using PHASER and SHELXE",
        "shortTitle": "Arcimboldo",
    },
    "arp_warp_classic": {
        "TASKTITLE": "ARP/wARP",
        "DESCRIPTION": "Build model (ARP/wARP classic)",
        "shortTitle": "ARP/wARP",
    },
    "buster": {
        "TASKTITLE": "Refinement - BUSTER",
        "DESCRIPTION": "Refine (BUSTER, Global Phasing Limited) with optional solvent/water update",
        "shortTitle": "BUSTER",
    },
    "ccp4mg_edit_model": {
        "TASKTITLE": "Interactive model preparation - CCP4mg and MrBUMP",
        "DESCRIPTION": "Identify MR models with MrBUMP, display and select with CCP4mg",
        "shortTitle": "CCP4mg MrBUMP",
    },
    "ccp4mg_edit_nomrbump": {
        "TASKTITLE": "Interactive selection of MR model components - CCP4mg",
        "DESCRIPTION": "Use CCP4mg to select components of a search model and output to i2 for MR",
        "shortTitle": "CCP4mg atom select",
    },
    "ccp4mg_general": {
        "TASKTITLE": "Molecular graphics visualization and figure creation - CCP4MG",
        "DESCRIPTION": "Interactive molecular graphics: visualization, figure preparation, analysis.",
        "shortTitle": "CCP4MG",
    },
    "chainsaw": {
        "TASKTITLE": "Truncate search model - CHAINSAW",
        "DESCRIPTION": "Truncate and renumber model prior to molecular replacement",
        "shortTitle": "CHAINSAW",
    },
    "chltofom": {
        "TASKTITLE": "Convert HL phase to and from Phi/Fom",
        "DESCRIPTION": "Interconvert phases between Hendrickson Lattman and Phi/FOM representation",
        "shortTitle": "Convert HL phase to and from Phi/Fom",
    },
    "cif2mtz": {
        "TASKTITLE": "Import merged mmCIF X-ray data",
        "DESCRIPTION": "Import merged mmCIF X-ray data e.g. from the PDB database",
        "shortTitle": "Import merged mmCIF X-ray data",
    },
    "clustalw": {
        "TASKTITLE": "Align sequences - CLUSTALW",
        "DESCRIPTION": "Align sequences using clustalw",
        "shortTitle": "CLUSTALW",
    },
    "cmapcoeff": {
        "TASKTITLE": "Calculate unusual map coefficients",
        "DESCRIPTION": "Compute map coefficients from set of observations and phases (cmapcoeff)",
        "shortTitle": "Calculate map coefficients",
    },
    "comit": {
        "TASKTITLE": "Calculate omit map to reduce model bias after molecular replacement or for validation",
        "DESCRIPTION": "Calculate omit map to reduce model bias",
        "shortTitle": "Calculate omit map",
    },
    "coordinate_selector": {
        "TASKTITLE": "Import a coordinate set - optional selection of subset",
        "DESCRIPTION": "Select (potentially complicated) subset from a coordinate set",
        "shortTitle": "Import coordinates",
    },
    "coot1": {
        "TASKTITLE": "Coot 1",
        "DESCRIPTION": "Interactive model building with Coot 1",
        "shortTitle": "Coot 1",
    },
    "coot_find_waters": {
        "TASKTITLE": "Find waters - COOT",
        "DESCRIPTION": "Find and filter waters based on electron density and contacts (non-interactive Coot)",
        "shortTitle": "Find waters - COOT",
    },
    "coot_rebuild": {
        "TASKTITLE": "Model building - COOT 0.9",
        "DESCRIPTION": "Interactive building (Coot 0.9)",
        "shortTitle": "COOT",
    },
    "coot_rsr_morph": {
        "TASKTITLE": "Real space refinement morphing with coot",
        "DESCRIPTION": "Real space refinement morphing with coot",
        "shortTitle": "Real space refinement morphing with coot",
    },
    "coot_script_lines": {
        "TASKTITLE": "Scripted model building - COOT",
        "DESCRIPTION": "Use scripts to fit sidechains, perform stepped refinement, fill and fit... (non-interactive Coot)",
        "shortTitle": "Scripted COOT",
    },
    "cpatterson": {
        "TASKTITLE": "Calculate Patterson map",
        "DESCRIPTION": "Calculate Patterson map (cpatterson)",
        "shortTitle": "Calculate Patterson map",
    },
    "crank2": {
        "TASKTITLE": "Automated structure solution - CRANK2 phasing and building",
        "DESCRIPTION": "CRANK2 experimental phasing pipeline",
        "shortTitle": "CRANK2",
        "RANK": 1,
    },
    "crank2_comb_phdmmb": {
        "TASKTITLE": "Model building",
        "DESCRIPTION": "Crank2 model building",
        "shortTitle": "Model building",
        "RANK": 1,
    },
    "crank2_dmfull": {
        "TASKTITLE": "Density modification",
        "DESCRIPTION": "Crank2 density modification",
        "shortTitle": "Density modification",
        "RANK": 1,
    },
    "crank2_faest": {
        "TASKTITLE": "FA estimation",
        "DESCRIPTION": "Crank2 FA estimation and other phasing preparations",
        "shortTitle": "FA estimation",
        "RANK": 1,
    },
    "crank2_handdet": {
        "TASKTITLE": "Hand determination",
        "DESCRIPTION": "Crank2 determination of handedness",
        "shortTitle": "Hand determination",
        "RANK": 1,
    },
    "crank2_mbref": {
        "TASKTITLE": "Model building",
        "DESCRIPTION": "Crank2 model building",
        "shortTitle": "Model building",
        "RANK": 1,
    },
    "crank2_phas": {
        "TASKTITLE": "Phasing",
        "DESCRIPTION": "Crank2 phasing",
        "shortTitle": "Phasing",
        "RANK": 1,
    },
    "crank2_phdmmb": {
        "TASKTITLE": "Density modification, poly-Ala tracing",
        "DESCRIPTION": "SHELXE density modification and poly-Ala tracing",
        "shortTitle": "Density modification, poly-Ala tracing",
        "RANK": 1,
    },
    "crank2_ref": {
        "TASKTITLE": "Refinement",
        "DESCRIPTION": "Crank2 refinement",
        "shortTitle": "Refinement",
        "RANK": 1,
    },
    "crank2_refatompick": {
        "TASKTITLE": "Iterative atom picking",
        "DESCRIPTION": "Crank2 substructure improvement",
        "shortTitle": "Iterative atom picking",
        "RANK": 1,
    },
    "crank2_substrdet": {
        "TASKTITLE": "Substructure detection",
        "DESCRIPTION": "Crank2 substructure detection",
        "shortTitle": "Substructure detection",
        "RANK": 1,
    },
    "csymmatch": {
        "TASKTITLE": "Match model to reference structure",
        "DESCRIPTION": "Match symmetry and origin of output model to reference structure (Csymmatch)",
        "shortTitle": "MATCH MODEL",
    },
    "ctruncate": {
        "TASKTITLE": "Convert intensities to amplitudes",
        "DESCRIPTION": "Convert reflection intensities to structure factors (ctruncate)",
        "shortTitle": "Convert intensities to amplitudes",
    },
    "density_calculator": {
        "TASKTITLE": "Density Calculator",
        "DESCRIPTION": "Calculate a map for a structure",
        "shortTitle": "Density Calculator",
    },
    "dials_image": {
        "TASKTITLE": "DIALS Image Viewer",
        "DESCRIPTION": "DIALS Image Viewer",
        "shortTitle": "DIALS Image Viewer",
    },
    "dials_rlattice": {
        "TASKTITLE": "DIALS Reciprocal Lattice Viewer",
        "DESCRIPTION": "DIALS Reciprocal Lattice Viewer",
        "shortTitle": "DIALS rLattice Viewer",
    },
    "dr_mr_modelbuild_pipeline": {
        "TASKTITLE": "Data Reduction, MR, Model build pipeline",
        "DESCRIPTION": "Data Reduction, MR, Model build pipeline",
        "shortTitle": "Data Rection, MR, model building",
        "RANK": 1,
    },
    "dui": {
        "TASKTITLE": "Integrate images with DIALS",
        "DESCRIPTION": "Launch DUI and capture output",
        "shortTitle": "DIALS User Interface",
    },
    "dummy": {
        "TASKTITLE": "Task interface demo code",
        "DESCRIPTION": "Demo code for developers - does not run",
        "shortTitle": "DumMee",
    },
    "editbfac": {
        "TASKTITLE": "Process Predicted Models",
        "DESCRIPTION": "Process Predicted Models - automatically process predicted models",
        "shortTitle": "Process Predicted Models",
    },
    "edstats": {
        "TASKTITLE": "Analyse agreement between model and density - EDSTATS",
        "DESCRIPTION": "Calculates real-space metrics for evaluating the agreement between model and density (Edstats, cfft)",
        "shortTitle": "EDSTATS",
        "RANK": 1,
    },
    "findmyseq": {
        "TASKTITLE": "Find My Sequence",
        "DESCRIPTION": "Find the most likely sequence in a database for a given map and model",
        "shortTitle": "Find My Sequence",
    },
    "freerflag": {
        "TASKTITLE": "Generate a Free R set",
        "DESCRIPTION": "Generate a Free R set for a complete set of reflection indices to a given resolution (FreeRflag)",
        "shortTitle": "Generate a Free R set",
        "RANK": 2,
    },
    "gesamt": {
        "TASKTITLE": "Structural alignment - Gesamt",
        "DESCRIPTION": "Superpose one protein structure on another",
        "shortTitle": "GESAMT",
    },
    "i2Dimple": {
        "TASKTITLE": "DIMPLE - Simple reflections + coordinates to map pipeline",
        "DESCRIPTION": "This task uses the DIMPLE pipeline to generate maps for a new dataset, provided a possible starting model is available.  For simple cases of isomorphous data, the pipeline will use rigid body refinement in REFMAC to 'tweak' the starting model. Where unit cells are incompatible, it will attempt automated molecular replacement.",
        "shortTitle": "DIMPLE",
    },
    "imosflm": {
        "TASKTITLE": "Integrate images with Mosflm",
        "DESCRIPTION": "Launch iMosflm and capture output",
        "shortTitle": "Integrate images with Mosflm",
    },
    "import_merged": {
        "TASKTITLE": "Import merged reflection data",
        "DESCRIPTION": "Import reflection data in any format, report on contents and create CCP4i2 data objects",
        "shortTitle": "Import merged",
    },
    "import_mosflm": {
        "TASKTITLE": "Import iMosflm X-ray data",
        "DESCRIPTION": "Import merged and unmerged X-ray reflections from Mosflm",
        "shortTitle": "Import iMosflm X-ray data",
    },
    "import_serial": {
        "TASKTITLE": "Import Serial",
        "DESCRIPTION": "Import merged data from CrystFEL",
        "shortTitle": "Import Serial Data",
    },
    "import_serial_pipe": {
        "TASKTITLE": "Import Serial Pipeline",
        "DESCRIPTION": "Import merged data from CrystFEL",
        "shortTitle": "Import Serial Data",
    },
    "lorestr_i2": {
        "TASKTITLE": "Low Resolution Refinement Pipeline (LORESTR)",
        "DESCRIPTION": "Automated Low Resolution Structure Refinement Pipeline (LORESTR)",
        "shortTitle": "LORESTR",
    },
    "matthews": {
        "TASKTITLE": "Estimate AU content",
        "DESCRIPTION": "Estimate number of molecules in the asymmetric unit and solvent content (Matthews_coeff)",
        "shortTitle": "Estimate AU content",
        "RANK": 2,
    },
    "mergeMtz": {
        "TASKTITLE": "Merge experimental data objects to MTZ",
        "DESCRIPTION": "Export 'old' style MTZ file",
        "shortTitle": "Merge to MTZ",
    },
    "metalCoord": {
        "TASKTITLE": "MetalCoord",
        "DESCRIPTION": "Generate restraints for metal-containing monomers",
        "shortTitle": "MetalCoord",
    },
    "modelcraft": {
        "TASKTITLE": "Autobuild with ModelCraft, Buccaneer and Nautilus",
        "DESCRIPTION": "Automated model building of protein, nucleic acid and water",
        "shortTitle": "ModelCraft",
    },
    "molrep_den": {
        "TASKTITLE": "Molecular replacement with electron density - MOLREP",
        "DESCRIPTION": "Use electron density as the search model (Molrep)",
        "shortTitle": "MOLREP with density",
    },
    "molrep_mr": {
        "TASKTITLE": "Molecular Replacement and refinement- MOLREP",
        "DESCRIPTION": "Molecular replacement (Molrep)",
        "shortTitle": "MOLREP MR",
    },
    "molrep_pipe": {
        "TASKTITLE": "Molecular Replacement and refinement - MOLREP",
        "DESCRIPTION": "Molecular replacement (Molrep)",
        "shortTitle": "MOLREP",
        "RANK": 1,
    },
    "molrep_selfrot": {
        "TASKTITLE": "Calculate self rotation function",
        "DESCRIPTION": "Evaluate data for anisotropy, optical resolution, pseudo translation and perform self-rotation function (Molrep)",
        "shortTitle": "Calculate self rotation function",
        "RANK": 2,
    },
    "morda_i2": {
        "TASKTITLE": "Automated molecular replacement - MORDA",
        "DESCRIPTION": "Molecular Replacement with Domains and Assemblies",
        "shortTitle": "MORDA",
        "RANK": 1,
    },
    "mosflm": {
        "TASKTITLE": "Integrate images - MOSFLM",
        "DESCRIPTION": "Use a script that you provide",
        "shortTitle": "Integrate images - MOSFLM",
    },
    "mrbump_basic": {
        "TASKTITLE": "Automated structure solution - MrBUMP",
        "DESCRIPTION": "Run a quick MrBUMP job with streamlined settings",
        "shortTitle": "Automated structure solution - MrBUMP",
        "RANK": 1,
    },
    "mrbump_model_prep": {
        "TASKTITLE": "MrBUMP model preparation",
        "DESCRIPTION": "Model preparation with MrBUMP for data reduction/MR/model build pipeline",
        "shortTitle": "MrBUMP model prep",
        "RANK": 1,
    },
    "mrparse": {
        "TASKTITLE": "MrParse",
        "DESCRIPTION": "Search online PDB and EBI-AFDB databases to find and process search models for use in molecular replacement",
        "shortTitle": "mrparse",
    },
    "mtzutils": {
        "TASKTITLE": "Add or delete MTZ columns",
        "shortTitle": "Add or delete MTZ columns",
    },
    "pairef": {
        "TASKTITLE": "Pairef",
        "DESCRIPTION": "Paired Refinement with Pairef",
        "shortTitle": "Pairef",
    },
    "parrot": {
        "TASKTITLE": "Density modification - PARROT",
        "DESCRIPTION": "Modify the electron density (Parrot)",
        "shortTitle": "PARROT",
    },
    "pdb_redo_api": {
        "TASKTITLE": "PDB-REDO Web services",
        "DESCRIPTION": "Refine structures using PDB-REDO Web service",
        "shortTitle": "PDB-REDO",
    },
    "pdbset_ui": {
        "TASKTITLE": "Scripted structure edits - Pdbset",
        "DESCRIPTION": "Structure edits with the pdbset program",
        "shortTitle": "PDBSET",
    },
    "pdbview_edit": {
        "TASKTITLE": "Edit PDB/CIF files by hand",
        "DESCRIPTION": "Edit PDB/CIF files by hand with the PdbView program",
        "shortTitle": "Edit PDB/CIF",
    },
    "phaser_EP": {
        "TASKTITLE": "SAD phasing - PHASER",
        "DESCRIPTION": "Complete a heavy atom model and calculate phases",
        "shortTitle": "SAD phasing - PHASER",
        "RANK": 1,
    },
    "phaser_EP_AUTO": {
        "TASKTITLE": "SAD phasing from heavy atom sites - PHASER",
        "DESCRIPTION": "Complete a heavy atom model and calculate phases",
        "shortTitle": "SAD phasing from heavy atom sites - PHASER",
    },
    "phaser_EP_LLG": {
        "TASKTITLE": "Anomalous map from coordinates - PHASER",
        "DESCRIPTION": "Calculate anomalous LLG map phased by coordinates to highlight anomalous scatterers (Phaser)",
        "shortTitle": "Anomalous map from coordinates - PHASER",
    },
    "phaser_MR_AUTO": {
        "TASKTITLE": "Molecular Replacement - Phaser",
        "DESCRIPTION": "Molecular replacement (Phaser)",
        "shortTitle": "Molecular Replacement - Phaser",
    },
    "phaser_MR_RNP": {
        "TASKTITLE": "Run rigid body refinement - PHASER",
        "DESCRIPTION": "Rigid body refinement usng PHASER in MR_RNP mode",
        "shortTitle": "Run rigid body refinement - PHASER",
    },
    "phaser_ensembler": {
        "TASKTITLE": "Build an ensemble for PHASER",
        "DESCRIPTION": "Compile assorted related structures into an ensemble for use in PHASER",
        "shortTitle": "Build an ensemble for PHASER",
    },
    "phaser_phil": {
        "TASKTITLE": "Phaser auto generated GUI",
        "DESCRIPTION": "Phaser auto generated GUI",
        "shortTitle": "phaser_phil",
    },
    "phaser_pipeline": {
        "TASKTITLE": "Expert Mode Molecular Replacement - PHASER",
        "DESCRIPTION": "Advanced MR options followed by refinement and rebuilding (Phaser, Refmac5, Coot)",
        "shortTitle": "Expert MR - PHASER",
        "RANK": 1,
    },
    "phaser_rnp_pipeline": {
        "TASKTITLE": "Rigid body refinement - PHASER",
        "DESCRIPTION": "Define rigid bodies for refinement (Phaser), fill partial residues (Coot) and refine (Refmac)",
        "shortTitle": "Rigid body PHASER",
    },
    "phaser_simple": {
        "TASKTITLE": "Basic Molecular Replacement - PHASER",
        "DESCRIPTION": "Simple MR with optional refinement and rebuilding (Phaser)",
        "shortTitle": "Basic MR - PHASER",
        "RANK": 1,
    },
    "phaser_singleMR": {
        "TASKTITLE": "Single Atom Molecular Replacement",
        "DESCRIPTION": "Perform Single Atom MR using Phaser",
        "shortTitle": "Single Atom MR",
    },
    "pisapipe": {
        "TASKTITLE": "Structure analysis with Pisa",
        "DESCRIPTION": "Analyse tertiary structure and interfaces of a protein",
        "shortTitle": "Structure analysis with Pisa",
    },
    "pointless_reindexToMatch": {
        "TASKTITLE": "Reindex reflections or change spacegroup",
        "DESCRIPTION": "Reindex: match to reference data/coordinates; change space group; analyse symmetry; or expand to P1 (Pointless)",
        "shortTitle": "Reindex reflections or change spacegroup",
        "RANK": 2,
    },
    "privateer": {
        "TASKTITLE": "Validation of carbohydrate structures - Privateer",
        "DESCRIPTION": "Validation, re-refinement and graphical analysis of carbohydrate structures",
        "shortTitle": "Privateer",
    },
    "prosmart": {
        "TASKTITLE": "ProSMART - Restraint generation and structural comparison",
        "shortTitle": "ProSMART - Restraint generation and structural comparison",
    },
    "prosmart_refmac": {
        "TASKTITLE": "Refinement - Refmacat/Refmac5",
        "DESCRIPTION": "Refine (Refmacat/Refmac5) with optional restraints (Prosmart, Platonyzer)",
        "shortTitle": "REFMAC5",
        "RANK": 1,
    },
    "pyphaser_mr": {
        "TASKTITLE": "MR using Phaser (pythonic)",
        "shortTitle": "MR using Phaser (pythonic)",
    },
    "qtpisa": {
        "TASKTITLE": "Interface and quaternary structure analysis - PISA",
        "DESCRIPTION": "Interface and assembly analysis (qtpisa)",
        "shortTitle": "PISA",
    },
    "scaleit": {
        "TASKTITLE": "Compare two or more datasets",
        "DESCRIPTION": "Compare two or more datasets by running SCALEIT",
        "shortTitle": "Compare two or more datasets",
    },
    "sculptor": {
        "TASKTITLE": "Truncate search model - SCULPTOR",
        "DESCRIPTION": "Truncate model prior to molecular replacement",
        "shortTitle": "SCULPTOR",
    },
    "servalcat_pipe": {
        "TASKTITLE": "Refinement - Servalcat",
        "DESCRIPTION": "Refinement against diffraction data or cryo-EM SPA map with optional restraints from ProSmart and/or MetalCoord",
        "shortTitle": "Servalcat",
        "RANK": 1,
    },
    "sheetbend": {
        "TASKTITLE": "Fast preliminary refinement of atomic model coordinates or temperature factors, including at low resolution",
        "DESCRIPTION": "Fast preliminary refinement of atomic model using the program sheetbend",
        "shortTitle": "Shift field refinement",
    },
    "shelx": {
        "TASKTITLE": "Automated structure solution - SHELXC/D/E phasing and building",
        "DESCRIPTION": "Experimental phasing pipeline SHELX (run via Crank2)",
        "shortTitle": "SHELX",
        "RANK": 1,
    },
    "shelxeMR": {
        "TASKTITLE": "Model building from Molecular Replacement solution using Shelxe",
        "DESCRIPTION": "Use Shelxe to attempt to improve (or verify) a solution from Molecular Replacement",
        "shortTitle": "SHELXE-MR",
    },
    "slicendice": {
        "TASKTITLE": "SliceNDice - Auto model processing and MR",
        "DESCRIPTION": "Automated processing of predicted or deposited search models and Molecular Replacement",
        "shortTitle": "SliceNDice - Auto model processing and MR",
        "RANK": 1,
    },
    "splitMtz": {
        "TASKTITLE": "Import and Split MTZ into experimental data objects",
        "DESCRIPTION": "Select groups of columns from the MTZ file (csplitmtz)",
        "shortTitle": "Import and Split MTZ",
    },
    "tableone": {
        "TASKTITLE": "Generate Table One",
        "DESCRIPTION": "Generate Table One for publications.",
        "shortTitle": "Table One",
    },
    "validate_protein": {
        "TASKTITLE": "Multimetric model validation - Iris",
        "DESCRIPTION": "Calculates per-residue metrics including B-factors, density fit quality, Ramachandran plots, rotamer outliers, and clashes (Iris-Validation & MolProbity)",
        "shortTitle": "Multimetric validation",
        "RANK": 1,
    },
    "workflow": {
    },
    "xia2_dials": {
        "TASKTITLE": "Automated integration of images with DIALS using xia2",
        "DESCRIPTION": "Select a directory containing images and integrate them",
        "shortTitle": "xia2/dials",
        "RANK": 1,
    },
    "xia2_multiplex": {
        "TASKTITLE": "Automated combination of datasets using xia2.multiplex",
        "DESCRIPTION": "Select previous xia2 runs",
        "shortTitle": "multiplex",
        "RANK": 1,
    },
    "xia2_ssx_reduce": {
        "TASKTITLE": "Reduction of serial datasets using xia2.ssx_reduce",
        "DESCRIPTION": "Select integrated data from xia2.ssx or dials.stills_process",
        "shortTitle": "xia2.ssx_reduce",
        "RANK": 1,
    },
    "xia2_xds": {
        "TASKTITLE": "Automated integration of images with XDS using xia2",
        "DESCRIPTION": "Select a directory containing images and integrate them",
        "shortTitle": "xia2/xds",
        "RANK": 1,
    },
    "zanuda": {
        "TASKTITLE": "Zanuda",
        "DESCRIPTION": "Space group validation",
        "shortTitle": "Zanuda",
    },
}
