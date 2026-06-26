/**
 * Moorhen Scene types.
 *
 * A Scene is a portable, human-editable description of how to render one
 * or more structures in Moorhen: which files to load, what domains to
 * recognise, what representations and colour rules to apply, and where the
 * camera should sit.
 *
 * Scenes are designed to be re-applied across different PDB files of the
 * same protein. Residue ranges in domain definitions are resolved against
 * the actual residues present in each loaded structure at apply-time, with
 * the resolver policy controlled by the `resolver` block.
 *
 * Two related artefacts:
 *   - `*.scene.yaml`  — this format; the only file a human writes/edits.
 *   - `*.session.json` — Moorhen-native backupSession, regenerable from
 *                       (scene, structure). Disposable cache.
 *
 * See `moorhen-scene.md` for the full grammar and worked examples.
 */

// --------------------------------------------------------------------------
// Top-level Scene
// --------------------------------------------------------------------------

/** Current schema version. Bump on breaking changes; readers may refuse
 *  to load scenes with a higher major version than they support. */
export const SCENE_SCHEMA_VERSION = 1 as const;

export interface MoorhenScene {
  /** Human-readable identifier for the scene. Free text. */
  scene: string;

  /** Schema version. Must equal SCENE_SCHEMA_VERSION for v1 readers. */
  version: number;

  /** Provenance: where this scene was authored. Never consulted by the
   *  resolver — only shown to humans and surfaced in bug reports. */
  authoredIn?: SceneProvenance;

  /** Files to load, named so that elements can reference them. */
  files?: SceneFileRef[];

  /** Superpositions to apply after fetching but before rendering. Each
   *  entry aligns one file (`move`) onto another (`onto`). Applied in
   *  declared order. */
  superpose?: SceneSuperpose[];

  /** Dictionary file names (from the `files:` block, with kind:
   *  dictionary) loaded globally — visible to every coordinate molecule
   *  in the scene. Use for cofactors, common buffer components, etc.
   *  Loaded before any coordinate molecules so coords containing those
   *  monomers parse correctly. */
  globalDictionaries?: string[];

  /** Reusable domain definitions, referenced by `colour: by-domain`
   *  inside elements. Hoisted to the top level so a multi-structure
   *  scene doesn't duplicate them per element. */
  domains?: SceneDomain[];

  /** Per-file rendering instructions: representations and colour rules. */
  elements?: SceneElement[];

  /** Map instances to render. Each references either an MTZ file from the
   *  `files:` block (kind: "mtz") with a column-label spec, or a real-space
   *  CCP4 map file (kind: "map", incl. masks) which needs no columns. Carries
   *  optional contour/colour/style; masks default to a translucent solid
   *  surface. */
  maps?: SceneMap[];

  /** Name of the map (from `maps:`) that should be Moorhen's active
   *  map after apply. Moorhen refines against the active map, so this
   *  is normally the best/calculated map (not a difference map). */
  activeMap?: string;

  /** Camera, clip, fog, background. Portable subset of Moorhen's
   *  viewDataSession; lighting/SSAO/shadow params intentionally omitted. */
  view?: SceneView;

  /** Apply-time policy. */
  resolver?: SceneResolverOptions;
}

// --------------------------------------------------------------------------
// Provenance
// --------------------------------------------------------------------------

export interface SceneProvenance {
  projectId?: string;       // UUID of the authoring project
  projectName?: string;     // human-readable name (advisory)
  createdAt?: string;       // ISO-8601 timestamp
  createdBy?: string;       // email or username
  ccp4i2Version?: string;   // for debugging
}

// --------------------------------------------------------------------------
// File references
// --------------------------------------------------------------------------

/**
 * A named file reference. Exactly one of {pdb, url, fileId+project,
 * job+param+project, path} should be set. The resolver first looks for an
 * already-loaded molecule that matches; failing that (and when the ref
 * carries enough info), the resolver fetches the coords and registers a
 * new molecule before binding.
 */
export interface SceneFileRef {
  /** Local name used by elements (e.g. "protein", "apo", "ref"). Unique
   *  within the files block. */
  name: string;

  /** What kind of file this ref points at. Defaults to "coordinates".
   *  "dictionary" refs are CIF monomer dicts (refmac/coot format) that
   *  the resolver loads into Coot's dictionary store rather than as
   *  separate molecules; they're then scoped to specific molecules via
   *  the per-element `dictionaries:` list. "mtz" refs are reflection
   *  data — referenced by entries in the `maps:` block which carry the
   *  column-label spec and contour/colour settings. "map" refs are
   *  real-space CCP4 map files (application/CCP4-map), including masks —
   *  also referenced by `maps:` entries, but loaded directly (no columns)
   *  via loadToCootFromMapData. */
  kind?: "coordinates" | "dictionary" | "mtz" | "map";

  /** PDB ID (4-letter or extended). Fetched via the PDBe proxy on apply
   *  if not already loaded. The most portable, share-friendly form for
   *  deposited structures. Only meaningful for `kind: coordinates`. */
  pdb?: string;

  /** Inline CIF text — only valid on `kind: dictionary` refs. Used by
   *  the lifter for dicts loaded from job outputs (no stable URL). The
   *  resolver hands the text straight to Coot's `read_dictionary_string`
   *  without any network round-trip. Multi-block dicts work in one shot
   *  because Coot parses every `data_comp_*` block in the input. */
  cifText?: string;

  /** Asset path inside a scene bundle (`.scene.zip`). Resolved against
   *  the in-memory asset map the editor populates when a .zip is opened.
   *  The portable shape for sharing scenes with their data attached —
   *  works for both coord and dictionary refs, and round-trips across
   *  machines without needing local filesystem access. */
  bundle?: string;

  /** Fully qualified URL — most portable across machines for non-PDB
   *  structures (refined coords on a shared server, e.g.). */
  url?: string;

  /** Absolute path on the local filesystem. Not portable. */
  path?: string;

  // -- ccp4i2 project-internal references --------------------------------

  /** UUID of the project containing the file. Required for `fileId` or
   *  `job`+`param` references. */
  projectId?: string;

  /** Project name (advisory; for human inspection only). */
  projectName?: string;

  /** Stable file id within the project. */
  fileId?: number;

  /** Job number within the project. Pair with `param`. */
  job?: number;

  /** Job parameter name, e.g. "XYZOUT". Pair with `job`. */
  param?: string;
}

// --------------------------------------------------------------------------
// Domains
// --------------------------------------------------------------------------

export interface SceneDomain {
  /** Free-text name used in `colour: by-domain` and in the resolver log. */
  name: string;

  /** CID selection — the preferred, general form. Any valid Coot CID:
   *
   *   - `"//F"`        — the whole of chain F.
   *   - `"//F/32-64"`  — residues 32-64 of chain F.
   *   - a wildcard chain (`//<star>/32-64`) — that range across every chain.
   *   - `"//A/(ALA,GLY)"`, `"//A/55/CA[C]"` — residue names, atoms, … things
   *     the chain+range form can't express.
   *
   *  When the CID is the `//chain/start-end` shape the resolver still clamps the
   *  range to present residues and warns (parity with `chain`+`range`); any other
   *  CID is passed straight to Coot. Use this in preference to `chain`+`range`. */
  selection?: string;

  /** @deprecated Legacy chain selector — use `selection`. Kept for a short
   *  migration window. Forms: `"A"` (single chain), `"*"` (every chain present),
   *  `["A","B","C"]` (explicit list). The resolver fans out across the resolved
   *  chains and clamps the range per-chain. */
  chain?: string | string[];

  /** @deprecated Legacy inclusive residue range "start-end" — use `selection`.
   *  Omitted ⇒ the whole chain. Clamped to present residues by the resolver. */
  range?: string;

  /** Hex colour, e.g. "#4b8bbe". */
  color: string;
}

// --------------------------------------------------------------------------
// Superpose
// --------------------------------------------------------------------------

/**
 * Align one loaded structure onto another. Run after fetch, before
 * representations/camera, so the saved view's quaternion is meaningful
 * against the aligned coords. Mutates the moving structure's display
 * transform in place; does not write to disk.
 *
 * Two methods, mirroring Moorhen's own API:
 *
 *   - `ssm`: secondary-structure matching. Cheap default; needs one
 *     chain id per side. Use for homologues or different conformations
 *     of the same chain.
 *
 *   - `lsq`: least-squares on explicit residue ranges. Use when you
 *     know the correspondences (e.g. matching specific binding-site
 *     residues across distantly related structures).
 */
export type SceneSuperpose = SceneSuperposeSsm | SceneSuperposeLsq;

export interface SceneSuperposeSsm {
  method: "ssm";
  /** File name (from the `files:` block) being transformed. */
  move: string;
  /** File name being aligned to (reference; unchanged). */
  onto: string;
  /** Chain id in the moving structure to use for the alignment. */
  movChain: string;
  /** Chain id in the reference structure to use for the alignment. */
  refChain: string;
}

export interface SceneSuperposeLsq {
  method: "lsq";
  move: string;
  onto: string;
  /** Per-range correspondences. Each entry pairs a residue range in the
   *  reference with one in the moving structure. Omit when using the
   *  `chain`+`range` shorthand below. */
  matches?: SceneLsqMatch[];
  /** Shorthand chain id, applied to both reference and moving structures.
   *  Use with `range:` when the simple "same chain, same residue numbers
   *  on both sides" case applies — saves writing a single-entry `matches`
   *  block. Mutually exclusive with `matches`. */
  chain?: string;
  /** Shorthand residue range "start-end", applied to both sides. Paired
   *  with `chain:`. */
  range?: string;
  /** Which atoms to fit:
   *   - `"all"`   — all atoms in the residue ranges
   *   - `"main"`  — main-chain atoms only (default; usually what you want)
   *   - `"ca"`    — Cα atoms only (fastest, most tolerant of differences) */
  matchType?: "all" | "main" | "ca";
}

export interface SceneLsqMatch {
  refChain: string;
  /** Inclusive residue range "start-end" in the reference. */
  refRange: string;
  movChain: string;
  /** Inclusive residue range "start-end" in the moving structure. */
  movRange: string;
}

// --------------------------------------------------------------------------
// Elements
// --------------------------------------------------------------------------

export interface SceneElement {
  /** Name of a file from the top-level `files:` block. */
  file: string;

  /** Dictionary file names (with `kind: dictionary`) to associate with
   *  this molecule specifically. Loaded globally first (so the coord
   *  parses) then re-associated with this molecule's molNo. This is the
   *  pattern that lets two molecules both call their ligand "LIG" with
   *  different chemistry — each gets its own scoped dictionary. */
  dictionaries?: string[];

  /** Representations to draw on this file. */
  representations?: SceneRepresentation[];
}

export interface SceneRepresentation {
  /** Moorhen RepresentationStyle, e.g. "CRs" (ribbons), "CBs" (sticks),
   *  "MolecularSurface", "ligands". See moorhen.RepresentationStyles
   *  for the full set of 27 styles. */
  style: string;

  /** CID selection, e.g. "//A" (chain A), "//STAR/LIG" (all ligands
   *  named LIG, where STAR is the wildcard). Default is the all-atoms
   *  wildcard CID. */
  selection?: string;

  /** Colour specification — see SceneColour for the variants. */
  colour?: SceneColour;

  /** Opacity in [0, 1] (1 = fully opaque). Maps to Moorhen's per-representation
   *  `nonCustomOpacity`; applies to surfaces (MolecularSurface, VdWSurface,
   *  gaussian, MetaBalls) as well as ribbons/sticks. Omitted ⇒ opaque. This is
   *  the opacity used when colour is scene-driven (rules / hex / named schemes);
   *  a hand-picked colour-picker colour would instead carry alpha in its RGBA. */
  alpha?: number;
}

/**
 * Colour specification. Shapes:
 *
 *   1. Hex literal:           colour: "#4b8bbe"
 *   2. Named scheme:          colour: by-domain | b-factor | af2-plddt | ...
 *   3. Per-selection list:    colour: [ { selection: "//A", colour: "#a08766" }, ... ]
 *   4. Raw escape hatch:      colour: { raw: { ruleType, args } }
 *
 * The per-selection list is the general "colour these CIDs these colours" form:
 * one entry per selection, each compiled to a single colour rule. A whole chain
 * (`//A`) and a residue range (`//A/121-130`) are the same shape at different
 * granularities — so this is what coot's default per-chain colouring lifts to
 * (one entry per chain, the colour coot assigned), instead of a single
 * misleading hex. `by-domain` is the named shortcut that compiles the top-level
 * `domains:` block into this same form. Named schemes map 1:1 to Moorhen's
 * multi-rule ruleTypes. The raw form preserves anything we can't lift.
 */
export type SceneColour =
  | string
  | SceneNamedColour
  | SceneColourSelection[]
  | SceneRawColour;

/** One entry of a per-selection colour list: a CID and the hex it gets. */
export interface SceneColourSelection {
  /** CID, e.g. "//A" (whole chain) or "//A/121-130" (residue range). */
  selection: string;
  /** Hex colour (#rrggbb or #rrggbbaa). */
  colour: string;
}

export type SceneNamedColour =
  | "by-domain"
  | "b-factor"
  | "b-factor-norm"
  | "af2-plddt"
  | "secondary-structure"
  | "jones-rainbow"
  | "mol-symm";

export interface SceneRawColour {
  raw: {
    ruleType: string;
    args: (string | number)[];
    isMultiColourRule?: boolean;
    applyColourToNonCarbonAtoms?: boolean;
  };
}

// --------------------------------------------------------------------------
// Maps
// --------------------------------------------------------------------------

/**
 * Column-label spec for reading an MTZ. Mirrors moorhen.selectedMtzColumns
 * but kept as a scene-local type so the YAML grammar isn't tied to the
 * installed Moorhen version. F + PHI are the minimum for a coefficient
 * map; Fobs/SigFobs/FreeR appear when the caller wants Moorhen to do its
 * own coefficient calculation (calcStructFact).
 */
export interface SceneMapColumns {
  /** Structure-factor label (e.g. "FWT", "F"). */
  F?: string;
  /** Phase label (e.g. "PHWT", "PHI"). */
  PHI?: string;
  /** Observed structure factor (when calcStructFact = true). */
  Fobs?: string;
  /** Sigma of observed F (when calcStructFact = true). */
  SigFobs?: string;
  /** Free-R flag column (when calcStructFact = true). */
  FreeR?: string;
  /** Weight column (FOM). */
  useWeight?: boolean;
  /** True iff Moorhen should compute coefficients on the fly
   *  rather than reading FWT/PHWT directly. */
  calcStructFact?: boolean;
}

/**
 * One Moorhen map. References an MTZ file from the `files:` block and
 * carries the column-label spec plus optional render state.
 *
 * Render-state fields are all optional — any omitted, Moorhen's
 * defaults apply at load time. The lifter only emits a field when its
 * captured value differs from those defaults, so the YAML stays small.
 */
export interface SceneMap {
  /** Local name, unique within the maps block. Used by `activeMap:` at
   *  the scene root and as a join key for the lifter/promoter. */
  name: string;

  /** Name of an entry in `files:` (kind: "mtz" or "map"). */
  file: string;

  /** Column-label spec for reading the referenced MTZ. Required for
   *  `kind: "mtz"` files; omitted for `kind: "map"` (real-space CCP4 maps
   *  / masks are read directly, no columns). */
  columns?: SceneMapColumns;

  /** Mark this entry as a mask (a real-space region map). Advisory: when
   *  set and style/contour/alpha are not otherwise specified, the resolver
   *  applies the translucent solid-surface mask defaults. Only meaningful
   *  for `kind: "map"` files. */
  isMask?: boolean;

  /** Difference-map flag. Maps with this set render both positive and
   *  negative contour lobes using `positiveColour` / `negativeColour`. */
  isDifference?: boolean;

  /** Contour level (rmsd-relative for Moorhen). */
  contourLevel?: number;

  /** Map radius (Å) — sphere around the camera origin to contour. */
  radius?: number;

  /** Opacity (0–1). */
  alpha?: number;

  /** Render style. "lit-lines" is Moorhen's prettier default for
   *  workstation-class GPUs. */
  style?: "lines" | "solid" | "lit-lines";

  /** Hex colour (#rrggbb or #rrggbbaa). For non-difference maps only;
   *  difference maps use `positiveColour` + `negativeColour`. */
  colour?: string;

  /** Hex colour for the positive contour lobe of a difference map. */
  positiveColour?: string;

  /** Hex colour for the negative contour lobe of a difference map. */
  negativeColour?: string;

  /** Whether the map is visible. Defaults to true. */
  visible?: boolean;
}

// --------------------------------------------------------------------------
// View
// --------------------------------------------------------------------------

export interface SceneView {
  /** Camera centre as an explicit Cartesian point. */
  origin?: [number, number, number];
  /** Camera centre derived from a selection's centroid — the resolver computes
   *  it at apply-time, so "centred on chain A" needs no coordinates. Takes
   *  precedence over `origin` when both are present. */
  centre?: SceneCentre;
  quat?: [number, number, number, number];
  zoom?: number;
  clipStart?: number;
  clipEnd?: number;
  fogStart?: number;
  fogEnd?: number;
  /** Clip/fog intent. Coot derives clip and fog from zoom and a shared pair of
   *  field depths, and recomputes them on zoom unless told not to. This is the
   *  stable, intent-level control over that. See SceneClip. When present it
   *  drives clip/fog (and the lock); the raw clipStart/End/fogStart/End above
   *  stay as an escape hatch. */
  clip?: SceneClip;
  /** Set the clip/fog DEPTH window to a selection's bounding sphere — a z-depth
   *  control only. Computed at apply-time from the selection's atoms. It does NOT
   *  move the camera: the slab brackets the current origin in depth, so to "show
   *  just chain A's region" pair it with `centre` on the same selection. `slab`
   *  and `centre` are independent and both settable. Takes precedence over
   *  `clip` when present (both control the same clip planes). */
  slab?: SceneSlab;
  background?: string;        // hex
}

/**
 * Clip/fog intent:
 *
 *   - `"auto"`              — let coot recompute clip+fog from zoom (its default).
 *   - `"lock"`             — freeze the current clip+fog so zoom won't change them.
 *   - `{ front, back }`    — set coot's field depths (zoom-independent depth of
 *                            field, in front of / behind the centre) and lock.
 *                            Drives clip AND fog together, the way coot does.
 *
 * Any explicit `clipStart/End/fogStart/End` are also locked on apply, so a
 * scene's clip sticks instead of being recomputed away on the next zoom.
 */
export type SceneClip = "auto" | "lock" | SceneClipFieldDepth;

export interface SceneClipFieldDepth {
  /** Depth of field in front of the view centre (coot default 8). */
  front: number;
  /** Depth of field behind the view centre (coot default 21). */
  back: number;
}

/** Selection whose centroid the camera centres on (see SceneView.centre). */
export interface SceneCentre {
  /** Name of a file from the top-level `files:` block. Optional when exactly one
   *  molecule is loaded — the resolver then defaults to it. */
  file?: string;
  /** CID selection within that file; omitted ⇒ the whole molecule. */
  selection?: string;
}

/**
 * Z-depth clip window for a selection (see SceneView.slab). The resolver walks
 * the selection's atoms for a bounding radius R and sets a symmetric clip/fog
 * depth of R + pad about the current origin. It does NOT centre — pair it with
 * `centre` on the same selection to frame it. Orientation-independent (a sphere);
 * once `orient` exists it can tighten to the depth along the view axis.
 */
export interface SceneSlab {
  /** Name of a file from the top-level `files:` block. Optional when exactly one
   *  molecule is loaded — the resolver then defaults to it. */
  file?: string;
  /** CID selection within that file; omitted ⇒ the whole molecule. */
  selection?: string;
  /** Extra Ångström added to the radius on each side (default 0). */
  pad?: number;
}

// --------------------------------------------------------------------------
// Resolver options
// --------------------------------------------------------------------------

export interface SceneResolverOptions {
  /**
   * Policy when a domain's residue range references residues that aren't
   * present in the loaded structure.
   *
   *   "clamp-and-log"  — clamp endpoints inward to the nearest present
   *                      residue, split across internal gaps, write a
   *                      sidecar log entry. Render silently. (Default.)
   *   "strict"         — fail to apply the scene; surface an error.
   */
  onMissingResidues?: "clamp-and-log" | "strict";
}

// --------------------------------------------------------------------------
// Type guards
// --------------------------------------------------------------------------

const NAMED_COLOURS: ReadonlySet<string> = new Set([
  "by-domain",
  "b-factor",
  "b-factor-norm",
  "af2-plddt",
  "secondary-structure",
  "jones-rainbow",
  "mol-symm",
]);

export function isSceneHexColour(c: SceneColour): c is string {
  return typeof c === "string" && c.startsWith("#");
}

export function isSceneNamedColour(c: SceneColour): c is SceneNamedColour {
  return typeof c === "string" && NAMED_COLOURS.has(c);
}

export function isSceneRawColour(c: SceneColour): c is SceneRawColour {
  return (
    typeof c === "object" &&
    c !== null &&
    "raw" in c &&
    typeof (c as SceneRawColour).raw === "object"
  );
}
