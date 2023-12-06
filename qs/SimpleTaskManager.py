import os

MODULE_ORDER = ['data_entry', 'data_processing', 'data_reduction', 'expt_phasing', 'bioinformatics',
                'molecular_replacement', 'density_modification', 'model_building', 'refinement', 'ligands',
                'validation', 'export', 'expt_data_utility', 'model_data_utility', 'developer_tools',
                ]

MODULE_ICONS = { 'data_entry'            : 'data_entry.svg',
                 'data_processing'       : 'data_processing.svg',
                 'data_reduction'        : 'data_reduction.svg',
                 'expt_phasing'          : 'expt_phasing.svg',
                 'bioinformatics'        : 'bioinformatics.svg',
                 'molecular_replacement' : 'molecular_replacement.svg',
                 'density_modification'  : 'density_modification.svg',
                 'model_building'        : 'model_building.svg',
                 'refinement'            : 'refinement.svg',
                 'ligands'               : 'ligands.svg',
                 'validation'            : 'validation.svg',
                 'export'                : 'export.svg',
                 'expt_data_utility'     : 'expt_data_utility.svg',
                 'model_data_utility'    : 'model_data_utility.svg',
                 'developer_tools'       : 'ccp4.svg',
                 'wrappers'              : 'ccp4.svg',
                 'demo'                  : 'ccp4.svg',
                 'test'                  : 'ccp4.svg'}

MODULE_TITLES = {'data_entry'            : 'Import merged data, AU contents, alignments or coordinates',
                 'data_processing'       : 'Integrate X-ray images',
                 'data_reduction'        : 'X-ray data reduction and analysis',
                 'expt_phasing'          : 'Experimental phasing',
                 'bioinformatics'        : 'Bioinformatics including model preparation for Molecular Replacement',
                 'molecular_replacement' : 'Molecular Replacement',
                 'density_modification'  : 'Density modification',
                 'model_building'        : 'Model building and Graphics',
                 'refinement'            : 'Refinement',
                 'ligands'               : 'Ligands',
                 'validation'            : 'Validation and analysis',
                 'export'                : 'Export and Deposition',
                 'expt_data_utility'     : 'Reflection data tools',
                 'model_data_utility'    : 'Coordinate data tools',
                 'developer_tools'       : 'Developer tools',
                 'wrappers'              : 'Wrappers',
                 'demo'                  : 'Demo code for developers only',
                 'test'                  : 'Test code for developers only'}

MODULE_DEFAULTS = {
        'test'        :    ['demo_copycell', 'refmac_wrk', 'crank2_comb_phdmmb', 'crank2_dmfull', 'crank2_faest', 'crank2_handdet', 'crank2_mbref', 'crank2_phas', 'crank2_phdmmb', 'crank2_ref', 'crank2_refatompick', 'crank2_substrdet', 'LidiaAcedrg', 'phaser_MR_RNP', 'pisapipe', 'reindex_minimtz', 'reindex_processed_data', 'uniqueify', 'newProject_fromMerged', 'phaser_mr', 'unique', 'aimless', 'arcimboldo', 'cif2mtz', 'convert2mtz', 'fft', 'freerflag', 'libcheck', 'molrep_mr', 'mosflm', 'mrbump_basic', 'mtzdump', 'mtzheader', 'mtzutils', 'phaser_phil', 'pyphaser_mr', 'scalepack2mtz', 'ShelxCD', 'ShelxCE', 'ShelxCECompareHands'],
        'model_data_utility'        :    ['csymmatch', 'gesamt', 'coordinate_selector', 'qtpisa', 'pdbview_edit'],
        'molecular_replacement'        :    ['mrbump_basic', 'phaser_simple', 'phaser_pipeline', 'molrep_pipe', 'molrep_den', 'csymmatch', 'parrot', 'phaser_rnp_pipeline', 'molrep_mr'],
        'data_entry'        :    ['import_merged', 'ProvideAsuContents', 'ProvideSequence', 'ProvideAlignment', 'AlternativeImportXIA2', 'coordinate_selector', 'splitMtz'],
        'data_processing'        :    ['xia2_dials', 'xia2_xds', 'imosflm'],
        'expt_data_utility'        :    ['pointless_reindexToMatch', 'phaser_EP_LLG', 'cmapcoeff', 'chltofom', 'cphasematch', 'ctruncate', 'splitMtz', 'freerflag', 'matthews'],
        'developer_tools'        :    ['MakeProjectsAndDoLigandPipeline', 'i2Dimple', 'MakeMonster', 'SyncToDjango', 'TestObsConversions'],
        'refinement'        :    ['prosmart_refmac', 'phaser_rnp_pipeline', 'lorestr_i2', 'ProvideTLS'],
        'density_modification'        :    ['acorn', 'parrot'],
        'export'        :    ['PrepareDeposit', 'mergeMtz'],
        'bioinformatics'        :    ['ccp4mg_edit_model', 'ccp4mg_edit_nomrbump', 'chainsaw', 'sculptor', 'phaser_ensembler', 'clustalw'],
        'data_reduction'        :    ['aimless_pipe', 'freerflag', 'matthews', 'molrep_selfrot'],
        'model_building'        :    ['buccaneer_build_refine_mr', 'coot_rebuild', 'coot_script_lines', 'coot_find_waters', 'nautilus_build_refine', 'arp_warp_classic', 'ccp4mg_general', 'shelxeMR'],
        'expt_phasing'        :    ['crank2', 'shelx', 'phaser_EP_AUTO'],
        'validation'        :    ['edstats', 'privateer_validate', 'qtpisa', 'privateer', 'validate_protein'],
        'ligands'        :    ['LidiaAcedrgNew', 'SubstituteLigand'],
                   }

MODULE_DEFAULTS['core'] = []
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['data_entry'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['model_data_utility'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['molecular_replacement'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['data_entry'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['data_processing'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['expt_data_utility'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['refinement'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['density_modification'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['export'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['bioinformatics'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['data_reduction'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['model_building'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['expt_phasing'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['validation'])
MODULE_DEFAULTS['core'].extend(MODULE_DEFAULTS['ligands'])


TASK_TITLES = {
        'arcimboldo'        :    'Ab initio Phasing with ARCIMBOLDO',
        'prosmart'        :    'Restraint generation and structural comparison (Prosmart)',
        'crank2'        :    'Automated structure solution - CRANK2 phasing and building',
        'fft'        :    'Export map',
        'phaser_pipeline'        :    'Expert Mode Molecular Replacement - PHASER',
        'cmapcoeff'        :    'Calculate unusual map coefficients',
        'phaser_rnp_pipeline'        :    'Rigid body refinement - PHASER',
        'crank2_dmfull'        :    'Density modification',
        'coordinate_selector'        :    'Import a coordinate set - optional selection of subset',
        'crank2_comb_phdmmb'        :    'Model building',
        'prosmart_martin'        :    'Simple prosmart for use in prosmart_refmac pipeline',
        'crank2_substrdet'        :    'Substructure detection',
        'phaser_simple'        :    'Basic Molecular Replacement - PHASER',
        'acedrgNew'        :    'acedrgNew',
        'SubstituteLigand'        :    'Automated solution of isomorphous ligand complex',
        'pdbset'        :    'PDBSet',
        'buccaneer_build_refine_mr'        :    'Autobuild protein - BUCCANEER',
        'chainsaw'        :    'Truncate search model - CHAINSAW',
        'PrepareDeposit'        :    'Prepare files for deposition',
        'cif2mtz'        :    'Import merged mmCIF X-ray data',
        'unique'        :    'Create dummy dataset',
        'arp_warp_classic'        :    'ARP/wARP',
        'uniqueify'        :    'Standardise reflection data file',
        'parrot'        :    'Parrot',
        'ZZPluginNameZZ'        :    'ZZLongTitleZZ',
        'mrbump_basic'        :    'MrBUMP Basic',
        'pdb_extract_wrapper'        :    'pdb_extract_wrapper',
        'ProvideAsuContents'        :    'Define AU contents',
        'xia2_run'        :    'Import data processing results from XIA2',
        'coot_rebuild'        :    'Manual model building - COOT',
        'pisa_xml'        :    'Structure analysis with Pisa',
        'cphasematch'        :    'Match and analyse phases to reference set',
        'import_mosflm'        :    'Import iMosflm X-ray data',
        'i2Dimple'        :    'DIMPLE - Simple reflections + coordinates to map pipeline',
        'pisa_list'        :    'Structure analysis with Pisa',
        'crank2_handdet'        :    'Hand determination',
        'phaser_EP_LLG'        :    'Anomalous map from coordinates - PHASER',
        'mtzheader'        :    'Read MTZ header',
        'ccp4mg_edit_model'        :    'Interactive model preparation - CCP4mg and MrBUMP',
        'xia2_xds'        :    'Automated integration of images with XDS using xia2',
        'mosflm'        :    'Integrate images - MOSFLM',
        'MakeProjectsAndDoLigandPipeline'        :    'Make projects and do ligand pipeline',
        'ctruncate'        :    'Convert intensities to amplitudes',
        'gesamt'        :    'Gesamt - structural alignment',
        'phaser_EP_AUTO'        :    'SAD phasing from heavy atom sites - PHASER',
        'acorn'        :    'ACORN - Phase Refinement with Dynamic Density Modification',
        'LidiaAcedrgNew'        :    'Make Ligand - Acedrg',
        'coot_find_waters'        :    'Find waters - COOT',
        'cad_copy_column'        :    'Copy an MTZ column between files',
        'prosmart_refmac'        :    'Refinement - REFMAC5',
        'xia2_aimless'        :    'Aimless in XIA2',
        'SyncToDjango'        :    'Exchange projects with Django Back end',
        'crank2_mbref'        :    'Model building',
        'demo_multi_mtzdump'        :    'Demo multi MTZ dump',
        'scalepack2mtz'        :    'Convert scalepack merged reflection file to MTZ',
        'refmac_martin'        :    'Basic refinement with Refmac5',
        'demo_copycell'        :    'Demo copy cell',
        'MakeMonster'        :    'Export monster mtz',
        'splitMtz'        :    'Import and Split MTZ into experimental data objects',
        'pisa_analyse'        :    'Structure analysis with Pisa',
        'Lidia'        :    'Lidia',
        'refmac'        :    'Refinement (Refmac5)',
        'xia2_pointless'        :    'Pointless in XIA2',
        'pdbview_edit'        :    'Edit PDB/CIF files by hand',
        'TestObsConversions'        :    'Test CCP4i2 observed data interconversions',
        'molrep_mr'        :    'Molecular replacement using Molrep',
        'molrep_den'        :    'Molecular replacement with electron density - MOLREP',
        'phaser_MR_AUTO'        :    'Molecular Replacement - Phaser',
        'pisapipe'        :    'Structure analysis with Pisa',
        'matthews'        :    'Estimate AU content',
        'crank2_phdmmb'        :    'Density modification, poly-Ala tracing',
        'lorestr_i2'        :    'Low Resolution Refinement Pipeline (LORESTR)',
        'qtpisa'        :    'Interface and quaternary structure analysis - PISA',
        'ccp4mg_edit_nomrbump'        :    'Interactive selection of MR model components - CCP4mg',
        'ProvideAlignment'        :    'Import an alignment',
        'pointless_reindexToMatch'        :    'Reindex or  change spacegroup',
        'reindex_minimtz'        :    'Reindex any miniMTZ',
        'ccp4mg_general'        :    'Molecular graphics visualization and figure creation - CCP4MG',
        'aimless_pipe'        :    'Data reduction - AIMLESS',
        'import_merged'        :    'Import merged reflection data',
        'ShelxCECompareHands'        :    'ShelxE Compare hands pipeline',
        'shelxeMR'        :    'ShelxeMR',
        'baverage'        :    'Average B over main and side chain atoms',
        'validate_protein'        :    'Analyse model geometry',
        'AlternativeImportXIA2'        :    'Import Xia2 results',
        'chltofom'        :    'Convert HL phase to and from Phi/Fom',
        'shelx'        :    'Automated structure solution - SHELXC/D/E phasing and building',
        'newProject_fromMerged'        :    'Import project starting data from merged data',
        'sculptor'        :    'Truncate search model - SCULPTOR',
        'xia2_ctruncate'        :    'CTruncate in XIA2',
        'clustalw'        :    'Align sequences - CLUSTALW',
        'pyphaser_mr'        :    'MR using Phaser (pythonic)',
        'molrep_pipe'        :    'Molecular Replacement and refinement - MOLREP',
        'imosflm'        :    'Integrate images - iMosflm',
        'refmac_wrk'        :    'Refmac Replica',
        'csymmatch'        :    'Match model to reference structure',
        'xia2_integration'        :    'Integration in XIA2',
        'ShelxCD'        :    'Find HA sites - SHELXC/D',
        'ShelxCE'        :    'Experimental phasing from heavy atom sites with ShelxE',
        'convert2mtz'        :    'Convert merged reflection file to MTZ',
        'privateer'        :    'Validation of carbohydrate structures',
        'privateer_validate'        :    'Validation of carbohydrate structures',
        'ProvideTLS'        :    'Import and/or edit TLS set definitions',
        'phaser_mr'        :    'MR using Phaser',
        'coot_script_lines'        :    'Scripted model building - COOT',
        'freerflag'        :    'Add a freeR flag',
        'acedrg'        :    'acedrg',
        'mtzutils'        :    'Add or delete MTZ columns',
        'crank2_faest'        :    'FA estimation',
        'guipreferences'        :    'CCP4i2 Preferences',
        'phaser_MR_RNP'        :    'Run rigid body refinement - PHASER',
        'mtzdump'        :    'MTZDump',
        'nautilus_build_refine'        :    'Autobuild RNA/DNA (Nautilus pipeline)',
        'xia2_dials'        :    'Automated integration of images with DIALS using xia2',
        'LidiaAcedrg'        :    'MAKE LIGAND',
        'libcheck'        :    'Manage dictionary libraries',
        'ProvideSequence'        :    'Import a sequence',
        'molrep_selfrot'        :    'Calculate self rotation function',
        'crank2_refatompick'        :    'Iterative atom picking',
        'reindex_processed_data'        :    'Reindex observed and FreeR data',
        'buccaneer_mr'        :    'Buccaneer_mr',
        'dummy'        :    'Task interface demo code',
        'aimless'        :    'Scale and merge dataset (AIMLESS)',
        'crank2_ref'        :    'Refinement',
        'mergeMtz'        :    'Merge experimental data objects to MTZ',
        'phaser_ensembler'        :    'Build an ensemble for PHASER',
        'phaser_phil'        :    'Phaser auto generated GUI',
        'edstats'        :    'Measure agreement between model and density',
        'ZZPipelineNameZZ'        :    'ZZLongTitleZZ',
        'crank2_phas'        :    'Phasing',
        'xia2_multiplex'     :    'Combine multiple datasets with xia2.multiplex'
}

TASK_ICONS = {
        'arcimboldo'        :    'arcimboldo.svg',
        'crank2'        :    'crank2.svg',
        'phaser_pipeline'        :    'phaser_pipeline.svg',
        'cmapcoeff'        :    'cmapcoeff.svg',
        'phaser_rnp_pipeline'        :    'phaser_rnp_pipeline.svg',
        'splitMtz'        :    'splitMtz.svg',
        'crank2_comb_phdmmb'        :    'crank2_comb_phdmmb.svg',
        'crank2_substrdet'        :    'crank2_substrdet.svg',
        'phaser_simple'        :    'phaser_simple.svg',
        'acedrgNew'        :    'acedrgNew.svg',
        'SubstituteLigand'        :    'SubstituteLigand.svg',
        'chainsaw'        :    'chainsaw.svg',
        'PrepareDeposit'        :    'PrepareDeposit.svg',
        'cif2mtz'        :    'cif2mtz.svg',
        'crank2_dmfull'        :    'crank2_dmfull.svg',
        'arp_warp_classic'        :    'ccp4.svg',
        'molrep_mr'        :    'ccp4.svg',
        'parrot'        :    'parrot.svg',
        'mrbump_basic'        :    'mrbump_basic.png',
        'phaser_MR_AUTO'        :    'phaser_MR_AUTO.svg',
        'ProvideAsuContents'        :    'ProvideAsuContents.svg',
        'coot_rebuild'        :    'coot_rebuild.svg',
        'cphasematch'        :    'cphasematch.svg',
        'import_mosflm'        :    'import_mosflm.svg',
        'i2Dimple'        :    'i2Dimple.svg',
        'crank2_handdet'        :    'crank2_handdet.svg',
        'phaser_EP_LLG'        :    'phaser_EP_LLG.svg',
        'ccp4mg_edit_model'        :    'ccp4mg_edit_model.svg',
        'xia2_xds'        :    'xia2_xds.svg',
        'mosflm'        :    'mosflm.svg',
        'MakeProjectsAndDoLigandPipeline'        :    'MakeProjectsAndDoLigandPipeline.svg',
        'ctruncate'        :    'ctruncate.svg',
        'gesamt'        :    'gesamt.svg',
        'phaser_EP_AUTO'        :    'phaser_EP_AUTO.svg',
        'acorn'        :    'acorn.svg',
        'matthews'        :    'matthews.svg',
        'coot_find_waters'        :    'coot_find_waters.svg',
        'prosmart_refmac'        :    'prosmart_refmac.svg',
        'SyncToDjango'        :    'SyncToDjango.svg',
        'crank2_mbref'        :    'crank2_mbref.svg',
        'molrep_selfrot'        :    'molrep_selfrot.svg',
        'demo_copycell'        :    'demo_copycell.svg',
        'MakeMonster'        :    'MakeMonster.svg',
        'Lidia'        :    'Lidia.svg',
        'ccp4mg_edit_nomrbump'        :    'ccp4mg_edit_nomrbump.svg',
        'TestObsConversions'        :    'TestObsConversions.svg',
        'pdb_extract_wrapper'        :    'pdb_extract_wrapper.svg',
        'pisapipe'        :    'pisapipe.svg',
        'LidiaAcedrgNew'        :    'LidiaAcedrgNew.svg',
        'crank2_phdmmb'        :    'crank2_phdmmb.svg',
        'qtpisa'        :    'qtpisa.svg',
        'pdbview_edit'        :    'pdbview_edit.svg',
        'ProvideAlignment'        :    'ProvideAlignment.svg',
        'pointless_reindexToMatch'        :    'pointless_reindexToMatch.svg',
        'ccp4mg_general'        :    'ccp4mg_general.svg',
        'aimless_pipe'        :    'aimless_pipe.svg',
        'lorestr_i2'        :    'ccp4.svg',
        'molrep_den'        :    'molrep_den.svg',
        'privateer_validate'        :    'privateer.svg',
        'ShelxCECompareHands'        :    'ShelxCECompareHands.svg',
        'shelxeMR'        :    'shelxeMR.svg',
        'validate_protein'        :    'validate_protein.svg',
        'AlternativeImportXIA2'        :    'AlternativeImportXIA2.svg',
        'chltofom'        :    'chltofom.svg',
        'shelx'        :    'shelx.svg',
        'coordinate_selector'        :    'coordinate_selector.svg',
        'sculptor'        :    'sculptor.svg',
        'buccaneer_build_refine_mr'        :    'buccaneer_build_refine_mr.svg',
        'clustalw'        :    'clustalw.svg',
        'molrep_pipe'        :    'molrep_pipe.svg',
        'imosflm'        :    'imosflm.svg',
        'csymmatch'        :    'csymmatch.svg',
        'import_merged'        :    'import_merged.svg',
        'ShelxCD'        :    'ShelxCD.svg',
        'ShelxCE'        :    'ShelxCE.svg',
        'ProvideSequence'        :    'ProvideSequence.svg',
        'privateer'        :    'privateer.svg',
        'ProvideTLS'        :    'ProvideTLS.svg',
        'coot_script_lines'        :    'coot_script_lines.svg',
        'freerflag'        :    'freerflag.svg',
        'acedrg'        :    'acedrg.svg',
        'xia2_dials'        :    'xia2_dials.svg',
        'crank2_faest'        :    'crank2_faest.svg',
        'phaser_MR_RNP'        :    'phaser_MR_RNP.svg',
        'nautilus_build_refine'        :    'nautilus_build_refine.svg',
        'LidiaAcedrg'        :    'LidiaAcedrg.svg',
        'crank2_refatompick'        :    'crank2_refatompick.svg',
        'reindex_processed_data'        :    'reindex_processed_data.svg',
        'dummy'        :    'dummy.svg',
        'crank2_ref'        :    'crank2_ref.svg',
        'mergeMtz'        :    'mergeMtz.svg',
        'phaser_ensembler'        :    'phaser_ensembler.svg',
        'phaser_phil'        :    'phaser_phil.svg',
        'edstats'        :    'edstats.svg',
        'crank2_phas'        :    'crank2_phas.svg',
        'xia2_multiplex'        :    'xia2_multiplex.svg',
}

TASK_KEYWORDS = {
        'crank2'                 :    ['experimental phasing','ep'],
        'privateer'              :    ['validation', 'sugars', 'carbohydrates'],
        'ProvideSequence'        :    ['import sequence','sequence import']
}

TASK_DESC = {
        'arcimboldo'        :    'To be done',
        'crank2'        :    'CRANK2 experimental phasing pipeline',
        'phaser_pipeline'        :    'Advanced MR options followed by refinement and rebuilding (Phaser, Refmac5, Coot)',
        'cmapcoeff'        :    'Compute map coefficients from set of observations and phases (cmapcoeff)',
        'phaser_rnp_pipeline'        :    'Define rigid bodies for refinement (Phaser), fill partial residues (Coot) and refine (Refmac)',
        'splitMtz'        :    'Select groups of columns from the MTZ file (csplitmtz)',
        'crank2_comb_phdmmb'        :    'Crank2 model building',
        'ProvideAlignment'        :    'Enter an alignment from a file or by cut and paste',
        'phaser_simple'        :    'Simple MR with optional refinement and rebuilding (Phaser)',
        'acedrgNew'        :    'Create a ligand dictionary with Acedrg',
        'SubstituteLigand'        :    'A ligand workflow, starting from merged or unmerged reflections, SMILES, and an isomorphous parent structure',
        'chainsaw'        :    'Truncate and renumber model prior to molecular replacement',
        'PrepareDeposit'        :    'Prepare files for deposition',
        'cif2mtz'        :    'Import merged mmCIF X-ray data e.g. from the PDB database',
        'crank2_dmfull'        :    'Crank2 density modification',
        'arp_warp_classic'        :    'Build model (ARP/wARP classic)',
        'molrep_mr'        :    'Molecular replacement (Molrep)',
        'molrep_den'        :    'Use electron density as the search model (Molrep)',
        'ZZPluginNameZZ'        :    'ZZDescriptionZZ',
        'parrot'        :    'Modify the electron density (Parrot)',
        'mrbump_basic'        :    'Run a quick MrBUMP job with streamlined settings',
        'pdb_extract_wrapper'        :    'Run pdb_extract',
        'ProvideAsuContents'        :    'Provide list of chain sequences- import sequence or from coordinate file',
        'coot_rebuild'        :    'Interactive building (Coot)',
        'cphasematch'        :    'Compare phases from different sources (with option change of origin/hand)',
        'import_mosflm'        :    'Import merged and unmerged X-ray reflections from Mosflm',
        'i2Dimple'        :    '''This task uses the DIMPLE pipeline to generate maps for a new dataset, provided a possible starting model is available.  For simple cases of isomorphous data, the pipeline will use rigid body refinement in REFMAC to 'tweak' the starting model. Where unit cells are incompatible, it will attempt automated molecular replacement.''',
        'crank2_handdet'        :    'Crank2 determination of handedness',
        'phaser_EP_LLG'        :    'Calculate anomalous LLG map phased by coordinates to highlight anomalous scatterers (Phaser)',
        'ccp4mg_edit_model'        :    'Identify MR models with MrBUMP, display and select with CCP4mg',
        'xia2_xds'        :    'Select a directory containing images and integrate them',
        'mosflm'        :    'Use a script that you provide',
        'MakeProjectsAndDoLigandPipeline'        :    'A task to generate projects and apply ligand pipeline for multiple datasets',
        'ctruncate'        :    'Convert reflection intensities to structure factors (ctruncate)',
        'gesamt'        :    'Superpose one protein structure on another',
        'phaser_EP_AUTO'        :    'Complete a heavy atom model and calculate phases',
        'acorn'        :    'Un-biased improvement of initial phases for high resolution data (1.5 Angstoms and better)',
        'matthews'        :    'Estimate number of molecules in the asymmetric unit and solvent content (Matthews_coeff)',
        'coot_find_waters'        :    'Find and filter waters based on electron density and contacts (non-interactive Coot)',
        'ProvideTLS'        :    'Enter TLS information to be used later in the project',
        'SyncToDjango'        :    'Synchronization utility for CCP4i2 projects in locations which have a django-based project archive set up',
        'crank2_mbref'        :    'Crank2 model building',
        'molrep_selfrot'        :    'Evaluate data for anisotropy, optical resolution, pseudo translation and perform self-rotation function (Molrep)',
        'demo_copycell'        :    'Demo pipeline - copy cell from MTZ file to PDB file',
        'MakeMonster'        :    'Build an MTZ file from data objects',
        'Lidia'        :    'Sketch a ligand',
        'pdbview_edit'        :    'Edit PDB/CIF files by hand with the PdbView program',
        'TestObsConversions'        :    'Exercise CCP4i2 Observed data representation conversions',
        'phaser_MR_AUTO'        :    'Molecular replacement (Phaser)',
        'pisapipe'        :    'Analyse tertiary structure and interfaces of a protein',
        'LidiaAcedrgNew'        :    'Generate a PDB file and dictionary (acedrg) from MOL file, SMILES string/file, or sketch (lidia).<br>Optionally match atom names to known structures.',
        'crank2_phdmmb'        :    'SHELXE density modification and poly-Ala tracing',
        'qtpisa'        :    'Interface and assembly analysis (qtpisa)',
        'ccp4mg_edit_nomrbump'        :    'Use CCP4mg to select components of a search model and output to i2 for MR',
        'crank2_substrdet'        :    'Crank2 substructure detection',
        'pointless_reindexToMatch'        :    'Reindex to match reference data/coordinates and/or change space group of reflections or Free R set (Pointless)',
        'ccp4mg_general'        :    'Interactive molecular graphics: visualization, figure preparation, analysis.',
        'aimless_pipe'        :    'Scale and analyse unmerged data and suggest space group (Pointless, Aimless, Ctruncate, FreeRflag)',
        'import_merged'        :    'Import reflection data in any format, report on contents and create CCP4i2 data objects',
        'ShelxCECompareHands'        :    'Phasing, density modification and autobuilding - both hands of the HA structure',
        'shelxeMR'        :    'Use Shelxe to attempt to improve (or verify) a solution from Molecular Replacement',
        'validate_protein'        :    'Calculate mean B-factors, Ramachandran plots and other metrics to aid in validation (Clipper)',
        'AlternativeImportXIA2'        :    'Harvest merged and unmerged files from selected XIA2 data reduction protocols',
        'chltofom'        :    'Interconvert phases between Hendrickson Lattman and Phi/FOM representation',
        'shelx'        :    'Experimental phasing pipeline SHELX (run via Crank2)',
        'coordinate_selector'        :    'Select (potentially complicated) subset from a coordinate set',
        'sculptor'        :    'Truncate model prior to molecular replacement',
        'buccaneer_build_refine_mr'        :    'Iterations of model building (Buccaneer) and refinement (Refmac5, Prosmart and Coot)',
        'clustalw'        :    'Align sequences using clustalw',
        'molrep_pipe'        :    'Molecular replacement (Molrep)',
        'imosflm'        :    'Launch iMosflm and capture output',
        'csymmatch'        :    'Match symmetry and origin of output model to reference structure (Csymmatch)',
        'lorestr_i2'        :    'Automated Low Resolution Structure Refinement Pipeline (LORESTR)',
        'ShelxCD'        :    'Find sites from SAD/MAD/SIR/SIRAS/RIP/RIPAS data',
        'ShelxCE'        :    'Phasing, density modification and autobuilding',
        'privateer'        :    'Checks stereochemistry, conformation and fit to density. Describes glycans (Privateer)',
        'privateer_validate'        :    'Checks stereochemistry, conformation and fit to density. Describes glycans (Privateer)',
        'prosmart_refmac'        :    'Refine (Refmac5) with optional restraints (Prosmart, Platonyzer)',
        'coot_script_lines'        :    'Use scripts to fit sidechains, perform stepped refinement, fill and fit... (non-interactive Coot)',
        'freerflag'        :    'Generate a Free R set for a complete set of reflection indices to a given resolution (FreeRflag)',
        'acedrg'        :    'Sketch a ligand',
        'xia2_dials'        :    'Select a directory containing images and integrate them',
        'crank2_faest'        :    'Crank2 FA estimation and other phasing preparations',
        'phaser_MR_RNP'        :    'Rigid body refinement usng PHASER in MR_RNP mode',
        'nautilus_build_refine'        :    'Iterations of model building (Nautilus) and refinement (Refmac5)',
        'LidiaAcedrg'        :    'Generate a PDB file and dictionary (acedrg) from MOL file, SMILES, or sketch (lidia)',
        'crank2_refatompick'        :    'Crank2 substructure improvement',
        'reindex_processed_data'        :    'Reindex processed data to correspond to a master dataset (Reindex)',
        'dummy'        :    'Demo code for developers - does not run',
        'ZZPipelineNameZZ'        :    'ZZDescriptionZZ',
        'crank2_ref'        :    'Crank2 refinement',
        'mergeMtz'        :    '''Export 'old' style MTZ file''',
        'phaser_ensembler'        :    'Compile assorted related structures into an ensemble for use in PHASER',
        'phaser_phil'        :    'Phaser auto generated GUI',
        'edstats'        :    'Calculates real-space metrics for evaluating the agreement between model and density (Edstats)',
        'ProvideSequence'        :    'Enter a sequence from a sequence file, from a PDB, or by cut and paste',
        'crank2_phas'        :    'Crank2 phasing',
        'xia2_multiplex'     :    'Combine data sets from previous xia2 runs'
}

def get_def_xml_for_task(taskName):
    if taskName == "ccp4mg_edit_model":
        xmlfn = os.path.join(os.path.dirname(__file__),"..","wrappers","ccp4mg_edit_model","script","ccp4mg_edit_model.def.xml")
    elif taskName == "aimless_pipe":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","aimless_pipe","script","aimless_pipe.def.xml")
    elif taskName == "buccaneer_build_refine_mr":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","buccaneer_build_refine_mr","script","buccaneer_build_refine_mr.def.xml")
    elif taskName == "crank2":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","crank2","script","crank2.def.xml")
    elif taskName == "crank2_createfree":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","crank2","wrappers","crank2_createfree","script","crank2_createfree.def.xml")
    elif taskName == "import_merged":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_merged","script","import_merged.def.xml")
    elif taskName == "import_xia2":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","script","import_xia2.def.xml")
    elif taskName == "xia2_aimless":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","wrappers","xia2_aimless","script","xia2_aimless.def.xml")
    elif taskName == "xia2_ctruncate":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","wrappers","xia2_ctruncate","script","xia2_ctruncate.def.xml")
    elif taskName == "xia2_integration":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","wrappers","xia2_integration","script","xia2_integration.def.xml")
    elif taskName == "xia2_pointless":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","wrappers","xia2_pointless","script","xia2_pointless.def.xml")
    elif taskName == "xia2_run":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","import_xia2","wrappers","xia2_run","script","xia2_run.def.xml")
    elif taskName == "LidiaAcedrg":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","LidiaAcedrg","script","LidiaAcedrg.def.xml")
    elif taskName == "LidiaAcedrgNew":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","LidiaAcedrgNew","script","LidiaAcedrgNew.def.xml")
    elif taskName == "MakeProjectsAndDoLigandPipeline":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","MakeProjectsAndDoLigandPipeline","script","MakeProjectsAndDoLigandPipeline.def.xml")
    elif taskName == "molrep_pipe":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","molrep_pipe","script","molrep_pipe.def.xml")
    elif taskName == "nautilus_build_refine":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","nautilus_build_refine","script","nautilus_build_refine.def.xml")
    elif taskName == "phaser_pipeline":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","script","phaser_pipeline.def.xml")
    elif taskName == "phaser_EP_AUTO":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","phaser_EP_AUTO","script","phaser_EP_AUTO.def.xml")
    elif taskName == "phaser_EP_LLG":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","phaser_EP_LLG","script","phaser_EP_LLG.def.xml")
    elif taskName == "phaser_MR":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","phaser_MR","script","phaser_MR.def.xml")
    elif taskName == "phaser_MR_AUTO":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","phaser_MR_AUTO","script","phaser_MR_AUTO.def.xml")
    elif taskName == "phaser_MR_RNP":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","phaser_MR_RNP","script","phaser_MR_RNP.def.xml")
    elif taskName == "pointless_reindexToMatch":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_pipeline","wrappers","pointless_reindexToMatch","script","pointless_reindexToMatch.def.xml")
    elif taskName == "phaser_rnp_pipeline":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_rnp_pipeline","script","phaser_rnp_pipeline.def.xml")
    elif taskName == "phaser_simple":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","phaser_simple","script","phaser_simple.def.xml")
    elif taskName == "pisapipe":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","pisapipe","script","pisapipe.def.xml")
    elif taskName == "pisa_analyse":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","pisapipe","wrappers","pisa_analyse","script","pisa_analyse.def.xml")
    elif taskName == "pisa":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","pisapipe","wrappers","pisa_list","script","pisa.def.xml")
    elif taskName == "pisa_list":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","pisapipe","wrappers","pisa_list","script","pisa_list.def.xml")
    elif taskName == "pisa_xml":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","pisapipe","wrappers","pisa_xml","script","pisa_xml.def.xml")
    elif taskName == "PrepareDeposit":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","PrepareDeposit","script","PrepareDeposit.def.xml")
    elif taskName == "hklin2cif":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","PrepareDeposit","wrappers","hklin2cif","script","hklin2cif.def.xml")
    elif taskName == "pdb_extract_wrapper":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","PrepareDeposit","wrappers","pdb_extract_wrapper","script","pdb_extract_wrapper.def.xml")
    elif taskName == "prosmart_refmac":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","script","prosmart_refmac.def.xml")
    elif taskName == "coot_find_waters":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","coot_find_waters","script","coot_find_waters.def.xml")
    elif taskName == "coot_fit_residues":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","coot_fit_residues","script","coot_fit_residues.def.xml")
    elif taskName == "coot_script_lines":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","coot_script_lines","script","coot_script_lines.def.xml")
    elif taskName == "coot_stepped_refine":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","coot_stepped_refine","script","coot_stepped_refine.def.xml")
    elif taskName == "reindex_minimtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","reindex_minimtz","script","reindex_minimtz.def.xml")
    elif taskName == "reindex_processed_data":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","prosmart_refmac","wrappers","reindex_processed_data","script","reindex_processed_data.def.xml")
    elif taskName == "shelx":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","shelx","script","shelx.def.xml")
    elif taskName == "SubstituteLigand":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","SubstituteLigand","script","SubstituteLigand.def.xml")
    elif taskName == "uniqueify":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","pipelines","uniqueify","script","uniqueify.def.xml")
    elif taskName == "acedrg":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","acedrg","script","acedrg.def.xml")
    elif taskName == "acedrgNew":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","acedrgNew","script","acedrgNew.def.xml")
    elif taskName == "acorn":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","acorn","script","acorn.def.xml")
    elif taskName == "aimless":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","aimless","script","aimless.def.xml")
    elif taskName == "AlternativeImportXIA2":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","AlternativeImportXIA2","script","AlternativeImportXIA2.def.xml")
    elif taskName == "arcimboldo":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","arcimboldo","script","arcimboldo.def.xml")
    elif taskName == "arp_warp_classic":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","arp_warp_classic","script","arp_warp_classic.def.xml")
    elif taskName == "baverage":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","baverage","script","baverage.def.xml")
    elif taskName == "boilerplate":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","boilerplate","script","boilerplate.def.xml")
    elif taskName == "buccaneer_mr":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","buccaneer_mr","script","buccaneer_mr.def.xml")
    elif taskName == "cad_copy_column":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","cad_copy_column","script","cad_copy_column.def.xml")
    elif taskName == "ccp4mg_edit_model":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ccp4mg_edit_model","script","ccp4mg_edit_model.def.xml")
    elif taskName == "ccp4mg_edit_nomrbump":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ccp4mg_edit_nomrbump","script","ccp4mg_edit_nomrbump.def.xml")
    elif taskName == "ccp4mg_general":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ccp4mg_general","script","ccp4mg_general.def.xml")
    elif taskName == "chainsaw":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","chainsaw","script","chainsaw.def.xml")
    elif taskName == "chltofom":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","chltofom","script","chltofom.def.xml")
    elif taskName == "cif2mtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","cif2mtz","script","cif2mtz.def.xml")
    elif taskName == "clustalw":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","clustalw","script","clustalw.def.xml")
    elif taskName == "cmapcoeff":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","cmapcoeff","script","cmapcoeff.def.xml")
    elif taskName == "convert2mtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","convert2mtz","script","convert2mtz.def.xml")
    elif taskName == "coordinate_selector":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","coordinate_selector","script","coordinate_selector.def.xml")
    elif taskName == "coot_rebuild":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","coot_rebuild","script","coot_rebuild.def.xml")
    elif taskName == "cphasematch":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","cphasematch","script","cphasematch.def.xml")
    elif taskName == "csymmatch":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","csymmatch","script","csymmatch.def.xml")
    elif taskName == "ctruncate":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ctruncate","script","ctruncate.def.xml")
    elif taskName == "edstats":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","edstats","script","edstats.def.xml")
    elif taskName == "fft":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","fft","script","fft.def.xml")
    elif taskName == "freerflag":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","freerflag","script","freerflag.def.xml")
    elif taskName == "gesamt":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","gesamt","script","gesamt.def.xml")
    elif taskName == "i2Dimple":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","i2Dimple","script","i2Dimple.def.xml")
    elif taskName == "imosflm":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","imosflm","script","imosflm.def.xml")
    elif taskName == "import_mosflm":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","import_mosflm","script","import_mosflm.def.xml")
    elif taskName == "libcheck":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","libcheck","script","libcheck.def.xml")
    elif taskName == "Lidia":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","Lidia","script","Lidia.def.xml")
    elif taskName == "lorestr_i2":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","lorestr_i2","script","lorestr_i2.def.xml")
    elif taskName == "lorestr_i2_autogenerated":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","lorestr_i2","script","lorestr_i2_autogenerated.def.xml")
    elif taskName == "MakeMonster":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","MakeMonster","script","MakeMonster.def.xml")
    elif taskName == "matthews":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","matthews","script","matthews.def.xml")
    elif taskName == "mergeMtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mergeMtz","script","mergeMtz.def.xml")
    elif taskName == "molrep_den":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","molrep_mr","script","molrep_den.def.xml")
    elif taskName == "molrep_mr":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","molrep_mr","script","molrep_mr.def.xml")
    elif taskName == "molrep_selfrot":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","molrep_selfrot","script","molrep_selfrot.def.xml")
    elif taskName == "morda_i2":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","morda_i2","script","morda_i2.def.xml")
    elif taskName == "mosflm":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mosflm","script","mosflm.def.xml")
    elif taskName == "mrbump_basic":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mrbump_basic","script","mrbump_basic.def.xml")
    elif taskName == "mtzdump":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mtzdump","script","mtzdump.def.xml")
    elif taskName == "mtzheader":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mtzheader","script","mtzheader.def.xml")
    elif taskName == "mtzutils":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","mtzutils","script","mtzutils.def.xml")
    elif taskName == "nautilus":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","nautilus","script","nautilus.def.xml")
    elif taskName == "parrot":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","parrot","script","parrot.def.xml")
    elif taskName == "pdbset":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","pdbset","script","pdbset.def.xml")
    elif taskName == "pdbview_edit":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","pdbview_edit","script","pdbview_edit.def.xml")
    elif taskName == "phaser_ensembler":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","phaser_ensembler","script","phaser_ensembler.def.xml")
    elif taskName == "phaser_mr":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","phaser_mr","script","phaser_mr.def.xml")
    elif taskName == "phaser_phil":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","phaser_phil","script","phaser_phil.def.xml")
    elif taskName == "pointless":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","pointless","script","pointless.def.xml")
    elif taskName == "privateer_validate":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","privateer","script","privateer.def.xml")
    elif taskName == "prosmart":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","prosmart","script","prosmart.def.xml")
    elif taskName == "ProvideAlignment":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ProvideAlignment","script","ProvideAlignment.def.xml")
    elif taskName == "ProvideAsuContents":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ProvideAsuContents","script","ProvideAsuContents.def.xml")
    elif taskName == "ProvideSequence":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ProvideSequence","script","ProvideSequence.def.xml")
    elif taskName == "ProvideTLS":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ProvideTLS","script","ProvideTLS.def.xml")
    elif taskName == "pyphaser_mr":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","pyphaser_mr","script","pyphaser_mr.def.xml")
    elif taskName == "qtpisa":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","qtpisa","script","qtpisa.def.xml")
    elif taskName == "refmac":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","refmac_i2","script","refmac.def.xml")
    elif taskName == "scalepack2mtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","scalepack2mtz","script","scalepack2mtz.def.xml")
    elif taskName == "sculptor":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","sculptor","script","sculptor.def.xml")
    elif taskName == "ShelxCD":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ShelxCDE","script","ShelxCD.def.xml")
    elif taskName == "ShelxCDEBase":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ShelxCDE","script","ShelxCDEBase.def.xml")
    elif taskName == "ShelxCE":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ShelxCDE","script","ShelxCE.def.xml")
    elif taskName == "ShelxCECompareHands":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","ShelxCDE","script","ShelxCECompareHands.def.xml")
    elif taskName == "shelxeMR":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","shelxeMR","script","shelxeMR.def.xml")
    elif taskName == "splitMtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","splitMtz","script","splitMtz.def.xml")
    elif taskName == "SyncToDjango":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","SyncToDjango","script","SyncToDjango.def.xml")
    elif taskName == "TestObsConversions":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","TestObsConversions","script","TestObsConversions.def.xml")
    elif taskName == "unique":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","unique","script","unique.def.xml")
    elif taskName == "validate_protein":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","validate_protein","script","validate_protein.def.xml")
    elif taskName == "x2mtz":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","x2mtz","script","x2mtz.def.xml")
    elif taskName == "xia2_dials":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","xia2_dials","script","xia2_dials.def.xml")
    elif taskName == "xia2_xds":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","xia2_xds","script","xia2_xds.def.xml")
    elif taskName == "xia2_multiplex":
        xmlfn= os.path.join(os.path.dirname(__file__),"..","wrappers","xia2_multiplex","script","xia2_multiplex.def.xml")
    else:
        return None
    return xmlfn
