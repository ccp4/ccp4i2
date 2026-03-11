/**
 * Client-side task category registry.
 *
 * This file is the single source of truth for how tasks are organised into
 * categories in the task-chooser UI.  The server provides a flat lookup of
 * task metadata (title, description, etc.) via the `task_tree` API; this
 * registry layers the categorisation on top.
 *
 * ## Adding a new task to an existing category
 *   Add its task name to the relevant `tasks` array below.
 *   Put it in `preferred` if it should appear near the top.
 *
 * ## Creating a new category
 *   Add a new entry to TASK_CATEGORIES with a unique `id`.
 *   Add a corresponding icon mapping in `lib/icons.ts`.
 *
 * ## Related files
 *   - `lib/icons.ts` — category icon mapping
 *   - `task-interfaces/task-container.tsx` — task name → custom React component
 */

// ── Category definition ────────────────────────────────────────────

export interface TaskCategory {
  /** Unique slug used as key and for icon lookup. */
  id: string;
  /** Human-readable description shown in the UI. */
  description: string;
  /**
   * Tasks that should appear first, in exactly this order.
   * Any task listed here that doesn't exist in the server lookup is skipped.
   */
  preferred: string[];
  /**
   * All tasks in this category (including preferred ones).
   * Non-preferred tasks are sorted alphabetically by title at display time.
   */
  tasks: string[];
}

// ── Registry ────────────────────────────────────────────────────────

export const TASK_CATEGORIES: TaskCategory[] = [
  {
    id: "data_entry",
    description: "Import merged data, AU contents, alignments or coordinates",
    preferred: ["import_merged", "ProvideAsuContents", "ProvideSequence"],
    tasks: [
      "import_merged",
      "ProvideAsuContents",
      "ProvideSequence",
      "ProvideAlignment",
      "AlternativeImportXIA2",
      "cif2mtz",
      "coordinate_selector",
      "import_mosflm",
      "import_serial",
      "import_serial_pipe",
      "import_xia2",
      "splitMtz",
    ],
  },
  {
    id: "data_processing",
    description: "Integrate X-ray images",
    preferred: ["xia2_dials", "xia2_xds"],
    tasks: [
      "xia2_dials",
      "xia2_xds",
      "imosflm",
      "mosflm",
      "dials_image",
      "dials_rlattice",
      "dui",
    ],
  },
  {
    id: "data_reduction",
    description: "X-ray data reduction and analysis",
    preferred: ["aimless_pipe"],
    tasks: [
      "aimless_pipe",
      "AUSPEX",
      "molrep_selfrot",
      "xia2_multiplex",
      "xia2_ssx_reduce",
    ],
  },
  {
    id: "alpha_fold",
    description: "AlphaFold and RoseTTAFold Utilities",
    preferred: [],
    tasks: ["editbfac", "slicendice"],
  },
  {
    id: "expt_phasing",
    description: "Experimental phasing",
    preferred: ["crank2"],
    tasks: ["crank2", "shelx", "ShelxCD", "phaser_EP_AUTO", "phaser_EP"],
  },
  {
    id: "bioinformatics",
    description:
      "Bioinformatics including model preparation for Molecular Replacement",
    preferred: ["ccp4mg_edit_model", "chainsaw", "mrparse"],
    tasks: [
      "ccp4mg_edit_model",
      "ccp4mg_edit_nomrbump",
      "chainsaw",
      "clustalw",
      "findmyseq",
      "matthews",
      "mrparse",
      "phaser_ensembler",
      "sculptor",
    ],
  },
  {
    id: "molecular_replacement",
    description: "Molecular Replacement",
    preferred: [
      "mrbump_basic",
      "phasertng_picard",
      "phaser_simple",
      "phaser_pipeline",
    ],
    tasks: [
      "mrbump_basic",
      "phasertng_picard",
      "phasertng_riker",
      "phaser_simple",
      "phaser_pipeline",
      "phaser_phil",
      "phaser_rnp_pipeline",
      "phaser_singleMR",
      "pyphaser_mr",
      "molrep_pipe",
      "molrep_den",
      "csymmatch",
      "AMPLE",
      "SIMBAD",
      "arcimboldo",
      "comit",
      "i2Dimple",
      "morda_i2",
    ],
  },
  {
    id: "density_modification",
    description: "Density modification",
    preferred: ["acorn", "parrot"],
    tasks: ["acorn", "parrot"],
  },
  {
    id: "model_building",
    description: "Model building and Graphics",
    preferred: ["modelcraft", "coot_rebuild"],
    tasks: [
      "modelcraft",
      "coot_rebuild",
      "coot1",
      "coot_script_lines",
      "coot_find_waters",
      "coot_find_ligand",
      "arp_warp_classic",
      "shelxeMR",
      "dr_mr_modelbuild_pipeline",
      "ccp4mg_general",
    ],
  },
  {
    id: "refinement",
    description: "Refinement",
    preferred: ["servalcat_pipe", "prosmart_refmac"],
    tasks: [
      "servalcat_pipe",
      "prosmart_refmac",
      "refmac",
      "prosmart",
      "ProvideTLS",
      "SubtractNative",
      "buster",
      "coot_rsr_morph",
      "lorestr_i2",
      "pairef",
      "pdb_redo_api",
      "sheetbend",
      "zanuda",
    ],
  },
  {
    id: "ligands",
    description: "Ligands",
    preferred: ["LidiaAcedrgNew", "MakeLink", "SubstituteLigand"],
    tasks: ["LidiaAcedrgNew", "MakeLink", "SubstituteLigand", "AcedrgLink"],
  },
  {
    id: "validation",
    description: "Validation and analysis",
    preferred: ["validate_protein"],
    tasks: [
      "validate_protein",
      "edstats",
      "metalCoord",
      "modelASUCheck",
      "pisapipe",
      "privateer",
      "qtpisa",
    ],
  },
  {
    id: "export",
    description: "Export and Deposition",
    preferred: ["PrepareDeposit"],
    tasks: [
      "PrepareDeposit",
      "adding_stats_to_mmcif_i2",
      "mergeMtz",
      "tableone",
    ],
  },
  {
    id: "expt_data_utility",
    description: "Reflection data tools",
    preferred: ["pointless_reindexToMatch"],
    tasks: [
      "pointless_reindexToMatch",
      "phaser_EP_LLG",
      "chltofom",
      "cmapcoeff",
      "cpatterson",
      "cphasematch",
      "ctruncate",
      "density_calculator",
      "freerflag",
      "mtzutils",
      "scaleit",
    ],
  },
  {
    id: "model_data_utility",
    description: "Coordinate data tools",
    preferred: ["gesamt"],
    tasks: [
      "gesamt",
      "add_fractional_coords",
      "areaimol",
      "pdbset_ui",
      "pdbview_edit",
    ],
  },
  {
    id: "developer_tools",
    description: "Developer tools",
    preferred: [],
    tasks: [
      "MakeMonster",
      "MakeProjectsAndDoLigandPipeline",
      "TestObsConversions",
      "mrparse_simple",
    ],
  },
];

// ── Tree builder ────────────────────────────────────────────────────

type TaskLookup = Record<
  string,
  { TASKTITLE?: string; DESCRIPTION?: string; shortTitle?: string; isAutogenerated?: boolean }
>;

export type TaskTreeEntry = [
  categoryId: string,
  description: string,
  taskNames: string[],
];

/**
 * Build the task tree by merging the client-side category registry with
 * the server-provided task lookup.  Tasks listed in a category but absent
 * from the server lookup are silently dropped (e.g. plugins not installed).
 * Tasks in the lookup but not in any category appear under "Uncategorized".
 */
export function buildTaskTree(lookup: TaskLookup): TaskTreeEntry[] {
  const validKeys = new Set(Object.keys(lookup));
  const assigned = new Set<string>();

  const tree: TaskTreeEntry[] = TASK_CATEGORIES.map((cat) => {
    const preferredSet = new Set(cat.preferred);

    // Preferred tasks first (in specified order), then the rest alphabetically
    const preferred = cat.preferred.filter((t) => validKeys.has(t));
    const rest = cat.tasks
      .filter((t) => !preferredSet.has(t) && validKeys.has(t))
      .sort((a, b) => {
        const titleA = (lookup[a]?.TASKTITLE ?? a).toLowerCase();
        const titleB = (lookup[b]?.TASKTITLE ?? b).toLowerCase();
        return titleA.localeCompare(titleB);
      });

    const ordered = [...preferred, ...rest];
    ordered.forEach((t) => assigned.add(t));

    return [cat.id, cat.description, ordered] as TaskTreeEntry;
  }).filter(([, , tasks]) => tasks.length > 0);

  // Uncategorized: tasks the server knows about but the registry doesn't list
  const uncategorized = [...validKeys]
    .filter((t) => !assigned.has(t))
    .sort((a, b) => {
      const titleA = (lookup[a]?.TASKTITLE ?? a).toLowerCase();
      const titleB = (lookup[b]?.TASKTITLE ?? b).toLowerCase();
      return titleA.localeCompare(titleB);
    });

  if (uncategorized.length > 0) {
    tree.push(["uncategorized", "Uncategorized", uncategorized]);
  }

  return tree;
}
