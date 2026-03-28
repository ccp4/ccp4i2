/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
// Maps task category module names to their representative icon filename in /svgicons/
export const MODULE_ICON_NAMES: Record<string, string> = {
  data_entry: "import_merged",
  data_processing: "xia2_dials",
  data_reduction: "aimless_pipe",
  bigpipes: "dr_mr_modelbuild_pipeline",
  alpha_fold: "mrparse",
  expt_phasing: "crank2",
  bioinformatics: "chainsaw",
  molecular_replacement: "phaser_simple",
  density_modification: "parrot",
  model_building: "coot_rebuild",
  refinement: "refmac",
  ligands: "acedrg",
  validation: "validate_protein",
  export: "export",
  expt_data_utility: "pointless",
  model_data_utility: "gesamt",
  developer_tools: "ccp4i2",
  wrappers: "ccp4i2",
  test: "ccp4i2",
  demo: "ccp4i2",
};

export const moduleCategoryIconSrc = (moduleName: string): string =>
  `/svgicons/${MODULE_ICON_NAMES[moduleName] ?? "ccp4i2"}.svg`;
