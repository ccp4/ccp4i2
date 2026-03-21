# CCP4i2 Test Coverage Tracker

Last updated: 2026-03-21

Coverage of GUI-visible tasks (from `client/renderer/lib/task-registry.ts`).
Tasks appearing in multiple categories are listed once, under their primary
category.  "Developer tools" and "Uncategorized" tasks are excluded.

Guide for writing new tests:
- [docs/writing-i2run-tests.md](writing-i2run-tests.md)
- [docs/writing-api-tests.md](writing-api-tests.md)

## Legend

| Symbol | Meaning                          |
|--------|----------------------------------|
| Y      | Test exists and passes           |
| D      | Test exists but disabled/skipped |
| N/A    | Not testable (interactive GUI)   |
| -      | No test                          |

### mmCIF column (coordinate input compatibility)

| Symbol | Meaning                                                        |
|--------|----------------------------------------------------------------|
| Y      | Wrapper is mmCIF-safe (preserves format or handles correctly)  |
| P      | PDB-only — auto-converts mmCIF→PDB (safe for ≤3-char codes)   |
| !      | PDB-only — broken for mmCIF (will fail or corrupt long codes)  |
| ·      | No coordinate input (MTZ/sequence only)                        |

## Summary

| Metric             | Count |
|--------------------|-------|
| GUI-visible tasks  | 81    |
| Interactive (N/A)  | 10    |
| i2run tested       | 65    |
| i2run disabled     | 7     |
| API tested         | 54    |
| API disabled       | 0     |
| Deprecated         | 1     |
| Not installed      | 1     |
| No tests at all    | 12    |
| Full coverage      | 45    |


## Data Entry

| Task                 | i2run | API | mmCIF | i2run test file                | API test file                |
|----------------------|-------|-----|-------|--------------------------------|------------------------------|
| import_merged        | Y     | Y   | ·     | test_import_merged.py          | test_data_reduction_api.py   |
| ProvideAsuContents   | Y     | Y   | ·     | test_asu_contents.py           | test_utilities_api.py        |
| ProvideSequence      | Y     | Y   | ·     | test_provide_sequence.py       | test_utilities_api.py        |
| ProvideAlignment     | Y     | Y   | ·     | test_provide_alignment.py      | test_utilities_api.py        |
| coordinate_selector  | Y     | Y   | Y     | test_coordinate_selector.py    | test_utilities_api.py        |
| AlternativeImportXIA2| D     | -   | ·     | test_alternative_import_xia2.py|                              |
| import_serial_pipe   | D     | -   | ·     | test_import_serial_pipe.py     |                              |

## Data Processing

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| xia2_dials           | Y     | -   | ·     | test_xia2.py               |                              |
| xia2_xds             | Y     | -   | ·     | test_xia2.py               |                              |
| xia2_multiplex       | D     | -   | ·     | test_xia2_multiplex.py     |                              |
| imosflm              | N/A   | N/A | ·     |                            |                              |
| dials_image          | N/A   | N/A | ·     |                            |                              |
| dials_rlattice       | N/A   | N/A | ·     |                            |                              |
| dui                  | N/A   | N/A | ·     |                            |                              |

## Data Reduction

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| aimless_pipe         | Y     | Y   | Y     | test_aimless.py            | test_data_reduction_api.py   |
| freerflag            | Y     | Y   | ·     | test_freerflag.py          | test_data_reduction_api.py   |
| matthews             | Y     | -   | ·     | test_matthews.py           |                              |
| molrep_selfrot       | Y     | Y   | P     | test_molrep_selfrot.py     | test_mr_pipelines_api.py     |
| AUSPEX               | Y     | -   | ·     | test_auspex.py             |                              |
| xia2_ssx_reduce      | D     | -   | ·     | test_xia2_ssx_reduce.py    |                              |

## Big Pipelines

| Task                      | i2run | API | mmCIF | i2run test file              | API test file                   |
|---------------------------|-------|-----|-------|------------------------------|---------------------------------|
| SubstituteLigand          | Y     | Y   | P     | test_substitute_ligand.py    | test_substitute_ligand_api.py   |
| dr_mr_modelbuild_pipeline | D     | -   | P     | test_drmrmb.py (disabled)    |                                 |

## AlphaFold Utilities

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                     |
|----------------------|-------|-----|-------|----------------------------|-----------------------------------|
| mrparse              | Y     | Y   | ·     | test_mrparse.py            | test_utilities_api.py             |
| editbfac             | Y     | Y   | P     | test_editbfac.py           | test_utilities_api.py             |
| arcimboldo           | Y     | -   | ·     | test_arcimboldo.py         |                                   |
| slicendice           | N/A   | N/A | P     |                            |                                   |

## Experimental Phasing

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| crank2               | Y     | Y   | Y     | test_crank2.py             | test_ep_pipelines_api.py     |
| shelx                | Y     | Y   | ·     | test_shelx.py              | test_ep_pipelines_api.py     |
| phaser_EP_AUTO       | Y     | Y   | P     | test_phaser_ep_auto.py     | test_ep_pipelines_api.py     |
| phaser_EP            | Y     | Y   | P     | test_phaser_ep.py          | test_ep_pipelines_api.py     |

## Bioinformatics

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| chainsaw             | Y     | -   | P     | test_chainsaw.py           |                              |
| sculptor             | Y     | -   | Y     | test_sculptor.py           |                              |
| phaser_ensembler     | Y     | -   | P     | test_phaser_ensembler.py   |                              |
| clustalw             | Y     | Y   | ·     | test_clustalw.py           | test_utilities_api.py        |
| findmyseq            | Y     | -   | P     | test_findmyseq.py         |                              |

## Molecular Replacement

| Task                 | i2run | API | mmCIF | i2run test file             | API test file                    |
|----------------------|-------|-----|-------|-----------------------------|----------------------------------|
| mrbump_basic         | Y     | Y   | P     | test_mrbump.py              | test_mr_pipelines_api.py         |
| phaser_simple        | Y     | Y   | P     | test_phaser_simple.py       | test_mr_pipelines_api.py         |
| phaser_pipeline      | Y     | Y   | P     | test_phaser_expert.py       | test_mr_pipelines_api.py         |
| molrep_pipe          | Y     | Y   | P     | test_molrep.py              | test_mr_pipelines_api.py         |
| molrep_den           | Y     | -   | P     | test_molrep_den.py          |                                  |
| parrot               | Y     | Y   | Y     | test_parrot.py              | test_model_building_api.py       |
| phaser_rnp_pipeline  | Y     | Y   | P     | test_phaser_rnp_pipeline.py | test_mr_pipelines_api.py         |
| AMPLE                | Y     | -   | P     | test_ample.py               |                                  |
| SIMBAD               | Y     | -   | P     | test_simbad.py              |                                  |
| morda_i2             | -     | -   | P     |                             |                                  |
| phaser_singleMR      | Y     | -   | P     | test_single_atom_mr.py      |                                  |
| comit                | Y     | Y   | ·     | test_comit.py               | test_utilities_api.py            |
| i2Dimple             | Y     | Y   | P     | test_dimple.py              | test_mr_pipelines_api.py         |

## Density Modification

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| acorn                | Y     | -   | P     | test_acorn.py              |                              |
| parrot               | Y     | Y   | Y     | test_parrot.py             | test_model_building_api.py   |

## Model Building

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| modelcraft           | Y     | Y   | Y     | test_modelcraft.py         | test_model_building_api.py   |
| coot1                | N/A   | N/A | Y     |                            |                              |
| coot_script_lines    | D     | -   | Y     | test_coot_script_lines.py  |                              |
| coot_find_waters     | Y     | Y   | Y     | test_find_waters.py        | test_utilities_api.py        |
| arp_warp_classic     | Y     | -   | P     | test_arpwarp.py            |                              |
| shelxeMR             | Y     | Y   | P     | test_shelxe_mr.py          | test_model_building_api.py   |

## Refinement

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                     |
|----------------------|-------|-----|-------|----------------------------|-----------------------------------|
| servalcat_pipe       | Y     | Y   | Y     | test_servalcat.py          | test_refinement_api.py            |
| prosmart_refmac      | Y     | Y   | P     | test_refmac.py             | test_refinement_api.py            |
| buster               | Y     | -   | P     | test_buster.py             |                                   |
| metalCoord           | Y     | Y   | Y     | test_metalcoord.py         | test_refinement_api.py            |
| ProvideTLS           | Y     | -   | P     | test_provide_tls.py        |                                   |
| coot_rsr_morph       | Y     | Y   | Y     | test_rsr_morph.py          | test_utilities_api.py             |
| pdb_redo_api         | D     | -   | Y     | test_pdb_redo_api.py       |                                   |
| sheetbend            | Y     | Y   | Y     | test_sheetbend.py          | test_refinement_api.py            |
| SubtractNative       | Y     | Y   | ·     | test_subtract_native.py    | test_utilities_api.py             |
| lorestr_i2           | Y     | Y   | Y     | test_lorestr.py            | test_refinement_api.py            |
| pairef               | Y     | -   | P     | test_pairef.py             |                                   |
| zanuda               | Y     | Y   | P     | test_zanuda.py             | test_refinement_api.py            |

## Ligands

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| LidiaAcedrgNew       | Y     | Y   | ·     | test_acedrg.py             | test_utilities_api.py        |
| MakeLink             | Y     | Y   | ·     | test_make_link.py          | test_utilities_api.py        |

## Validation

| Task                 | i2run | API | mmCIF | i2run test file            | API test file                |
|----------------------|-------|-----|-------|----------------------------|------------------------------|
| validate_protein     | Y     | Y   | Y     | test_validate.py           | test_utilities_api.py        |
| edstats              | Y     | Y   | P     | test_edstats_wrapper.py    | test_structure_analysis_api.py |
| privateer            | Y     | Y   | Y     | test_privateer.py          | test_structure_analysis_api.py |
| modelASUCheck        | Y     | -   | P     | test_model_asu_check.py    |                              |

## Export

| Task                     | i2run | API | mmCIF | i2run test file            | API test file                |
|--------------------------|-------|-----|-------|----------------------------|------------------------------|
| adding_stats_to_mmcif_i2 | D     | -   | ·     | test_adding_stats.py       |                              |
| mergeMtz                 | Y     | -   | ·     | test_merge_mtz.py          |                              |

## Reflection Data Tools

| Task                      | i2run | API | mmCIF | i2run test file            | API test file                     |
|---------------------------|-------|-----|-------|----------------------------|-----------------------------------|
| pointless_reindexToMatch  | Y     | Y   | ·     | test_pointless_reindex.py  | test_data_reduction_api.py        |
| phaser_EP_LLG             | Y     | -   | ·     | test_phaser_ep_llg.py      |                                   |
| cmapcoeff                 | Y     | -   | ·     | test_cmapcoeff.py          |                                   |
| chltofom                  | Y     | Y   | ·     | test_chltofom.py           | test_data_reduction_api.py        |
| cphasematch               | Y     | -   | ·     | test_cphasematch.py        |                                   |
| ctruncate                 | D     | Y   | ·     | test_ctruncate.py (skipped)| test_data_reduction_api.py        |
| splitMtz                  | Y     | Y   | ·     | test_split_mtz.py          | test_data_reduction_api.py        |
| scaleit                   | Y     | Y   | ·     | test_scaleit.py            | test_data_reduction_api.py        |
| cpatterson                | Y     | -   | ·     | test_cpatterson.py         |                                   |
| density_calculator        | Y     | -   | Y     | test_density_calculator.py |                                   |

## Coordinate Data Tools

| Task                 | i2run | API | mmCIF | i2run test file               | API test file                |
|----------------------|-------|-----|-------|-------------------------------|------------------------------|
| csymmatch            | Y     | Y   | Y     | test_csymmatch.py             | test_utilities_api.py        |
| gesamt               | Y     | Y   | Y     | test_gesamt.py                | test_structure_analysis_api.py |
| pdbset_ui            | Y     | -   | P     | test_pdbset_ui.py             |                              |
| add_fractional_coords| Y     | -   | Y     | test_add_fractional_coords.py |                              |


## Excluded from Testing

### Interactive GUI Tasks (Not Testable)

These tasks launch interactive graphical programs (CCP4MG, Coot, QtPISA, etc.)
and cannot be tested via i2run or the API.

| Task                 | Category               |
|----------------------|------------------------|
| ccp4mg_edit_model    | AlphaFold Utilities    |
| ccp4mg_edit_nomrbump | Bioinformatics         |
| ccp4mg_general       | Model Building         |
| coot_rebuild         | Model Building         |
| coot1                | Model Building         |
| pdbview_edit         | Coordinate Data Tools  |
| qtpisa               | Validation             |
| imosflm              | Data Processing        |
| dials_image          | Data Processing        |
| dials_rlattice       | Data Processing        |
| dui                  | Data Processing        |
| slicendice           | AlphaFold Utilities    |

### Not Installed

| Task                 | Category               | Notes                              |
|----------------------|------------------------|------------------------------------|
| morda_i2             | Molecular Replacement  | MoRDa not distributed with CCP4   |

### Deprecated Tasks

| Task                 | Replacement              | Notes                              |
|----------------------|--------------------------|------------------------------------|
| PrepareDeposit       | adding_stats_to_mmcif_i2 | Moved to Uncategorized in task-chooser with "do not use" warning |


## Known issues

| Task | Suite | Issue | Notes |
|------|-------|-------|-------|
| ctruncate | i2run | Skipped | KeywordExtractor bug: `get_merged_metadata` is None for ctruncate's container |
| adding_stats_to_mmcif_i2 | i2run | Skipped | Requires prior refinement job in same project (complex setup) |
| pdb_redo_api | i2run | Skipped | Requires external PDB-REDO web service |

## mmCIF Compatibility Summary

| Status | Count | Tasks |
|--------|-------|-------|
| Y (mmCIF-safe) | 19 | coordinate_selector, crank2, parrot, coot_find_waters, coot_script_lines, servalcat_pipe, metalCoord, coot_rsr_morph, pdb_redo_api, sheetbend, lorestr_i2, validate_protein, privateer, sculptor, csymmatch, density_calculator, add_fractional_coords, gesamt, modelcraft |
| P (auto-converts) | 19 | molrep_selfrot, SubstituteLigand, dr_mr_modelbuild_pipeline, editbfac, phaser_EP_AUTO, phaser_EP, phaser_ensembler, mrbump_basic, phaser_simple, phaser_pipeline, molrep_pipe, molrep_den, phaser_rnp_pipeline, AMPLE, SIMBAD, i2Dimple, arp_warp_classic, prosmart_refmac, zanuda, edstats, chainsaw, shelxeMR |
| ! (broken) | 0 | (none) |
| · (no coords) | 26 | (MTZ/sequence/ligand-only tasks) |

## Remaining gaps (no i2run test possible)

| Task                      | Category               | Reason                                    |
|---------------------------|------------------------|-------------------------------------------|
| morda_i2                  | Molecular Replacement  | MoRDa not distributed with CCP4           |
