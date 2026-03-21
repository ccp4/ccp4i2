/**
 * Client-side task category registry.
 *
 * This file is the **single source of truth** for which tasks are exposed in
 * the task-chooser UI, which folder each task lives in, and the preferred
 * ordering within each folder.
 *
 * The structure mirrors the Qt GUI's task tree as defined by MODULE_ORDER,
 * MODULE_TITLES and per-task TASKMODULE attributes on the main branch
 * (see task_menu_structure.json for the reference extraction).
 *
 * The server provides a flat lookup of task metadata (title, description,
 * etc.) via the `task_tree` API; this registry layers the categorisation on
 * top.  Tasks listed here but absent from the server lookup are silently
 * dropped (e.g. plugins not installed).
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
 *   - `task_menu_structure.json` — reference extraction from Qt main branch
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
//
// Aligned with Qt GUI main-branch task tree.  Tasks marked [django-only]
// in comments are exposed here but have no Qt GUI equivalent.

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
      "coordinate_selector",
      "AlternativeImportXIA2",
      "import_serial_pipe",
    ],
  },
  {
    id: "data_processing",
    description: "Integrate X-ray images",
    preferred: ["xia2_dials", "xia2_xds"],
    tasks: [
      "xia2_dials",
      "xia2_xds",
      "xia2_multiplex",
      "imosflm",
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
      "freerflag",
      "matthews",
      "molrep_selfrot",
      "AUSPEX",
      "xia2_ssx_reduce",
    ],
  },
  {
    id: "bigpipes",
    description: "Data to complete structure solution including ligand fitting",
    preferred: ["SubstituteLigand", "dr_mr_modelbuild_pipeline"],
    tasks: [
      "SubstituteLigand",
      "dr_mr_modelbuild_pipeline",
    ],
  },
  {
    id: "alpha_fold",
    description: "AlphaFold and RoseTTAFold Utilities",
    preferred: ["ccp4mg_edit_model", "mrparse", "editbfac"],
    tasks: [
      "ccp4mg_edit_model",
      "mrparse",
      "editbfac",
      "arcimboldo",
      "slicendice",
    ],
  },
  {
    id: "expt_phasing",
    description: "Experimental phasing",
    preferred: ["crank2"],
    tasks: [
      "crank2",
      "shelx",
      "phaser_EP_AUTO",
      "phaser_EP",
    ],
  },
  {
    id: "bioinformatics",
    description:
      "Bioinformatics including model preparation for Molecular Replacement",
    preferred: ["ccp4mg_edit_model", "chainsaw"],
    tasks: [
      "ccp4mg_edit_model",
      "ccp4mg_edit_nomrbump",
      "chainsaw",
      "sculptor",
      "phaser_ensembler",
      "clustalw",
      "findmyseq",
    ],
  },
  {
    id: "molecular_replacement",
    description: "Molecular Replacement",
    preferred: [
      "mrbump_basic",
      "phaser_simple",
      "phaser_pipeline",
    ],
    tasks: [
      "mrbump_basic",
      "phaser_simple",
      "phaser_pipeline",
      "molrep_pipe",
      "molrep_den",
      "csymmatch",
      "parrot",
      "phaser_rnp_pipeline",
      "AMPLE",
      "SIMBAD",
      "morda_i2",
      "phaser_singleMR",
      "comit",
      "i2Dimple",
      "arcimboldo",
    ],
  },
  {
    id: "density_modification",
    description: "Density modification",
    preferred: ["acorn"],
    tasks: [
      "acorn",
      "parrot",
    ],
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
      "buster",
      "phaser_rnp_pipeline",
      "metalCoord",
      "ProvideTLS",
      "coot_rsr_morph",
      "pdb_redo_api",
      "sheetbend",
      "SubtractNative",
      "lorestr_i2",
      "pairef",
      "zanuda",
    ],
  },
  {
    id: "ligands",
    description: "Ligands",
    preferred: ["LidiaAcedrgNew", "SubstituteLigand", "MakeLink"],
    tasks: [
      "LidiaAcedrgNew",
      "SubstituteLigand",
      "MakeLink",
    ],
  },
  {
    id: "validation",
    description: "Validation and analysis",
    preferred: ["validate_protein"],
    tasks: [
      "validate_protein",
      "edstats",
      "privateer",
      "modelASUCheck",
      "qtpisa",
    ],
  },
  {
    id: "export",
    description: "Export and Deposition",
    preferred: ["adding_stats_to_mmcif_i2"],
    tasks: [
      "adding_stats_to_mmcif_i2",
      "mergeMtz",
    ],
  },
  {
    id: "expt_data_utility",
    description: "Reflection data tools",
    preferred: ["pointless_reindexToMatch"],
    tasks: [
      "pointless_reindexToMatch",
      "phaser_EP_LLG",
      "cmapcoeff",
      "chltofom",
      "cphasematch",
      "ctruncate",
      "splitMtz",
      "scaleit",
      "cpatterson",
      "density_calculator",
    ],
  },
  {
    id: "model_data_utility",
    description: "Coordinate data tools",
    preferred: ["csymmatch", "gesamt"],
    tasks: [
      "csymmatch",
      "gesamt",
      "coordinate_selector",
      "qtpisa",
      "pdbview_edit",
      "pdbset_ui",
      "editbfac",
      "add_fractional_coords",
    ],
  },
  {
    id: "developer_tools",
    description: "Developer tools",
    preferred: [],
    tasks: [
      "MakeMonster",
      "TestObsConversions",
      "MakeProjectsAndDoLigandPipeline",
    ],
  },
];

// ── Deprecated tasks ─────────────────────────────────────────────────
//
// Tasks listed here will show a prominent warning in the task-chooser card.
// They remain available (via Uncategorized) but are discouraged.

export const DEPRECATED_TASKS: Record<string, string> = {
  PrepareDeposit: "Superseded by adding_stats_to_mmcif_i2 — do not use",
};

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
