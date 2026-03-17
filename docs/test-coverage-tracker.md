# CCP4i2 Test Coverage Tracker

Last updated: 2026-03-17

Coverage of GUI-visible tasks (from `client/renderer/lib/task-registry.ts`).
Tasks appearing in multiple categories are listed once, under their primary
category.  "Developer tools" and "Uncategorized" tasks are excluded.

Guide for writing new tests: [docs/writing-i2run-tests.md](writing-i2run-tests.md)

## Legend

| Symbol | Meaning                          |
|--------|----------------------------------|
| Y      | Test exists and passes           |
| D      | Test exists but disabled/skipped |
| -      | No test                          |

## Summary

| Metric             | Count |
|--------------------|-------|
| GUI-visible tasks  | 63    |
| i2run tested       | 35    |
| i2run disabled     | 2     |
| API tested         | 25    |
| No tests at all    | 24    |
| Full coverage      | 22    |


## Data Entry

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| import_merged        | Y     | Y   | test_import_merged.py      | test_data_reduction_api.py   |
| ProvideAsuContents   | Y     | Y   | test_asu_contents.py       | test_utilities_api.py        |
| ProvideSequence      | -     | -   |                            |                              |
| ProvideAlignment     | -     | -   |                            |                              |
| coordinate_selector  | Y     | Y   | test_coordinate_selector.py| test_utilities_api.py        |
| AlternativeImportXIA2| -     | -   |                            |                              |

## Data Processing

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| xia2_dials           | Y     | -   | test_xia2.py               |                              |
| xia2_xds             | Y     | -   | test_xia2.py               |                              |

## Data Reduction

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| aimless_pipe         | Y     | Y   | test_aimless.py            | test_data_reduction_api.py   |
| freerflag            | Y     | Y   | test_freerflag.py          | test_data_reduction_api.py   |
| matthews             | -     | -   |                            |                              |
| molrep_selfrot       | -     | -   |                            |                              |
| xia2_ssx_reduce      | -     | -   |                            |                              |

## Big Pipelines

| Task                      | i2run | API | i2run test file              | API test file                   |
|---------------------------|-------|-----|------------------------------|---------------------------------|
| SubstituteLigand          | Y     | Y   | test_substitute_ligand.py    | test_substitute_ligand_api.py   |
| dr_mr_modelbuild_pipeline | D     | -   | test_drmrmb.py (disabled)    |                                 |

## AlphaFold Utilities

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| ccp4mg_edit_model    | -     | -   |                            |                              |
| mrparse              | Y     | Y   | test_mrparse.py            | test_utilities_api.py        |
| editbfac             | Y     | -   | test_editbfac.py           |                              |
| arcimboldo           | Y     | -   | test_arcimboldo.py         |                              |

## Experimental Phasing

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| crank2               | Y     | Y   | test_crank2.py             | test_ep_pipelines_api.py     |
| shelx                | Y     | Y   | test_shelx.py              | test_ep_pipelines_api.py     |
| phaser_EP_AUTO       | -     | -   |                            |                              |
| phaser_EP            | Y     | Y   | test_phaser_ep.py          | test_ep_pipelines_api.py     |

## Bioinformatics

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| ccp4mg_edit_nomrbump | -     | -   |                            |                              |
| chainsaw             | -     | -   |                            |                              |
| sculptor             | -     | -   |                            |                              |
| phaser_ensembler     | -     | -   |                            |                              |
| clustalw             | -     | -   |                            |                              |

## Molecular Replacement

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| mrbump_basic         | Y     | Y   | test_mrbump.py             | test_mr_pipelines_api.py     |
| phaser_simple        | Y     | Y   | test_phaser_simple.py      | test_mr_pipelines_api.py     |
| phaser_pipeline      | Y     | -   | test_phaser_expert.py      |                              |
| molrep_pipe          | Y     | Y   | test_molrep.py             | test_mr_pipelines_api.py     |
| molrep_den           | -     | -   |                            |                              |
| parrot               | Y     | Y   | test_parrot.py             | test_model_building_api.py   |
| phaser_rnp_pipeline  | Y     | -   | test_phaser_rnp_pipeline.py|                              |
| AMPLE                | Y     | -   | test_ample.py              |                              |
| SIMBAD               | -     | -   |                            |                              |
| morda_i2             | -     | -   |                            |                              |
| comit                | Y     | -   | test_comit.py              |                              |
| i2Dimple             | Y     | Y   | test_dimple.py             | test_mr_pipelines_api.py     |

## Model Building

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| modelcraft           | Y     | Y   | test_modelcraft.py         | test_model_building_api.py   |
| coot_rebuild         | -     | -   |                            |                              |
| coot_script_lines    | -     | -   |                            |                              |
| coot_find_waters     | Y     | Y   | test_find_waters.py        | test_utilities_api.py        |
| arp_warp_classic     | Y     | -   | test_arpwarp.py            |                              |
| shelxeMR             | Y     | Y   | test_shelxe_mr.py          | test_model_building_api.py   |
| ccp4mg_general       | -     | -   |                            |                              |

## Refinement

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| servalcat_pipe       | Y     | Y   | test_servalcat.py          | test_refinement_api.py       |
| prosmart_refmac      | Y     | Y   | test_refmac.py             | test_refinement_api.py       |
| metalCoord           | Y     | -   | test_metalcoord.py         |                              |
| coot_rsr_morph       | Y     | Y   | test_rsr_morph.py          | test_utilities_api.py        |
| pdb_redo_api         | -     | -   |                            |                              |
| sheetbend            | Y     | Y   | test_sheetbend.py          | test_refinement_api.py       |
| SubtractNative       | Y     | -   | test_subtract_native.py    |                              |
| lorestr_i2           | D     | -   | test_lorestr.py (disabled) |                              |
| zanuda               | -     | -   |                            |                              |

## Ligands

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| LidiaAcedrgNew       | Y     | Y   | test_acedrg.py             | test_utilities_api.py        |
| MakeLink             | Y     | -   | test_make_link.py          |                              |

## Validation

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| validate_protein     | Y     | Y   | test_validate.py           | test_utilities_api.py        |
| edstats              | -     | -   |                            |                              |
| privateer            | -     | -   |                            |                              |
| qtpisa               | -     | -   |                            |                              |

## Export

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| PrepareDeposit       | -     | -   |                            |                              |
| adding_stats_to_mmcif_i2 | -  | -   |                            |                              |
| mergeMtz             | -     | -   |                            |                              |

## Reflection Data Tools

| Task                      | i2run | API | i2run test file            | API test file                |
|---------------------------|-------|-----|----------------------------|------------------------------|
| pointless_reindexToMatch  | -     | -   |                            |                              |
| phaser_EP_LLG             | -     | -   |                            |                              |
| cmapcoeff                 | -     | -   |                            |                              |
| chltofom                  | -     | -   |                            |                              |
| cphasematch               | -     | -   |                            |                              |
| ctruncate                 | -     | -   |                            |                              |
| splitMtz                  | Y     | -   | test_split_mtz.py          |                              |
| scaleit                   | Y     | -   | test_scaleit.py            |                              |
| cpatterson                | -     | -   |                            |                              |
| density_calculator        | -     | -   |                            |                              |

## Coordinate Data Tools

| Task                 | i2run | API | i2run test file            | API test file                |
|----------------------|-------|-----|----------------------------|------------------------------|
| csymmatch            | Y     | Y   | test_csymmatch.py          | test_utilities_api.py        |
| gesamt               | -     | -   |                            |                              |
| pdbview_edit         | -     | -   |                            |                              |
| add_fractional_coords| -     | -   |                            |                              |


## Tasks with no test coverage (priority list)

These 24 tasks have neither i2run nor API tests:

| Priority | Task                      | Category               | Notes                            |
|----------|---------------------------|------------------------|----------------------------------|
| HIGH     | chainsaw                  | Bioinformatics         | Widely used MR prep tool         |
| HIGH     | phaser_EP_AUTO            | Experimental Phasing   | Key EP pipeline                  |
| HIGH     | zanuda                    | Refinement             |                                  |
| HIGH     | pdb_redo_api              | Refinement             | External API dependency          |
| HIGH     | edstats                   | Validation             |                                  |
| HIGH     | privateer                 | Validation             |                                  |
| MED      | matthews                  | Data Reduction         | Simple calculation               |
| MED      | molrep_selfrot            | Data Reduction         |                                  |
| MED      | molrep_den                | Molecular Replacement  |                                  |
| MED      | coot_rebuild              | Model Building         | Requires Coot binary             |
| MED      | coot_script_lines         | Model Building         | Requires Coot binary             |
| MED      | qtpisa                    | Validation             |                                  |
| MED      | pointless_reindexToMatch  | Refl. Data Tools       | Tested indirectly via phaser_rnp |
| MED      | ctruncate                 | Refl. Data Tools       |                                  |
| MED      | gesamt                    | Coord. Data Tools      |                                  |
| LOW      | ProvideSequence           | Data Entry             | Simple data entry                |
| LOW      | ProvideAlignment          | Data Entry             | Simple data entry                |
| LOW      | AlternativeImportXIA2     | Data Entry             |                                  |
| LOW      | xia2_ssx_reduce           | Data Reduction         | Needs SSX data                   |
| LOW      | ccp4mg_edit_model         | AlphaFold / Bioinf.    | Requires CCP4MG                  |
| LOW      | ccp4mg_edit_nomrbump      | Bioinformatics         | Requires CCP4MG                  |
| LOW      | sculptor                  | Bioinformatics         |                                  |
| LOW      | phaser_ensembler          | Bioinformatics         |                                  |
| LOW      | clustalw                  | Bioinformatics         | External binary                  |
| LOW      | SIMBAD                    | Molecular Replacement  | Large search — slow              |
| LOW      | morda_i2                  | Molecular Replacement  |                                  |
| LOW      | ccp4mg_general            | Model Building         | Requires CCP4MG                  |
| LOW      | PrepareDeposit            | Export                 |                                  |
| LOW      | adding_stats_to_mmcif_i2  | Export                 |                                  |
| LOW      | mergeMtz                  | Export                 |                                  |
| LOW      | phaser_EP_LLG             | Refl. Data Tools       |                                  |
| LOW      | cmapcoeff                 | Refl. Data Tools       |                                  |
| LOW      | chltofom                  | Refl. Data Tools       |                                  |
| LOW      | cphasematch               | Refl. Data Tools       |                                  |
| LOW      | cpatterson                | Refl. Data Tools       |                                  |
| LOW      | density_calculator        | Refl. Data Tools       |                                  |
| LOW      | pdbview_edit              | Coord. Data Tools      |                                  |
| LOW      | add_fractional_coords     | Coord. Data Tools      |                                  |
