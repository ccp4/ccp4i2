/**
 * Moorhen Scene — Zod schema, portable CORE.
 *
 * Single source of truth for the scene format's structure. TS types are
 * derived (`z.infer`), runtime validation IS this schema (`.parse`), and the
 * published JSON Schema is generated from it (see ./index.ts `buildJsonSchema`).
 *
 * This module is the **upstreamable core**: only Moorhen-meaningful, portable
 * concepts. ccp4i2-specific, deployment-coupled references (relativeUrl,
 * fileId, job+param, projectId) live in ./dialect.ts, which extends the file
 * ref defined here. See MOORHEN_SCENES_SCHEMA_V1_DESIGN.md §3, §5.
 *
 * Descriptions are kept terse on purpose — they flow into the JSON Schema and
 * the generated authoring brief, where tokens cost (design doc §7a).
 */
import { z } from "zod";

import { SCENE_SCHEMA_VERSION } from "../../types/moorhen-scene";

// --- primitives -----------------------------------------------------------

/** #rrggbb or #rrggbbaa. */
export const HexColour = z
  .string()
  .regex(/^#(?:[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$/, "must be #rrggbb or #rrggbbaa")
  .describe("hex colour #rrggbb or #rrggbbaa");

/** Inclusive residue range "start-end", e.g. "32-64". */
export const ResidueRange = z
  .string()
  .regex(/^-?\d+-{1}-?\d+$|^-?\d+$/, 'residue range "start-end"');

// --- colour ---------------------------------------------------------------

export const NamedColour = z.enum([
  "by-domain",
  "b-factor",
  "b-factor-norm",
  "af2-plddt",
  "secondary-structure",
  "jones-rainbow",
  "mol-symm",
]);

/** One entry of a per-selection colour list: a CID and the hex it gets. */
export const ColourSelection = z
  .object({
    selection: z.string().describe('CID, e.g. "//A" or "//A/121-130"'),
    colour: HexColour,
  })
  .strict();

/** Escape hatch: a raw Moorhen colour rule we can't otherwise express. */
export const RawColour = z
  .object({
    raw: z
      .object({
        ruleType: z.string(),
        args: z.array(z.union([z.string(), z.number()])),
        isMultiColourRule: z.boolean().optional(),
        applyColourToNonCarbonAtoms: z.boolean().optional(),
      })
      .strict(),
  })
  .strict();

/**
 * Colour spec. Shapes: hex literal | named scheme | per-selection list |
 * raw rule. See design doc / moorhen-scene.md.
 */
export const Colour = z.union([
  HexColour,
  NamedColour,
  z.array(ColourSelection),
  RawColour,
]);

// --- file refs (CORE: portable only) --------------------------------------

export const FileKind = z.enum(["coordinates", "dictionary", "mtz", "map"]);

/**
 * Base file-ref fields shared by core and dialect. Kept as a plain object
 * (not refined) so the dialect can `.extend()` it; the "exactly one source"
 * rule is applied per-profile in the refined exports below.
 */
export const CoreFileRefShape = z.object({
  name: z.string().describe("local name referenced by elements/maps"),
  kind: FileKind.optional().describe('default "coordinates"'),
  pdb: z.string().optional().describe("PDB id; fetched via proxy on apply"),
  url: z.string().url().optional().describe("absolute URL (portable)"),
  bundle: z.string().optional().describe("asset path inside a .scene.zip"),
  cifText: z.string().optional().describe("inline CIF (dictionary refs only)"),
});

/** Portable file ref: exactly one of pdb | url | bundle | cifText. */
export const CoreFileRef = CoreFileRefShape.strict().superRefine((ref, ctx) => {
  const sources = [ref.pdb, ref.url, ref.bundle, ref.cifText].filter(
    (v) => v != null && v !== "",
  );
  if (sources.length !== 1) {
    ctx.addIssue({
      code: z.ZodIssueCode.custom,
      message: "set exactly one of: pdb, url, bundle, cifText",
    });
  }
  if (ref.cifText && ref.kind !== "dictionary") {
    ctx.addIssue({
      code: z.ZodIssueCode.custom,
      path: ["cifText"],
      message: "cifText is only valid on kind: dictionary",
    });
  }
  if (ref.pdb && ref.kind === "dictionary") {
    ctx.addIssue({
      code: z.ZodIssueCode.custom,
      path: ["pdb"],
      message: "pdb is for coordinates, not dictionaries",
    });
  }
});

// --- domains --------------------------------------------------------------

export const Domain = z
  .object({
    name: z.string().describe("used by colour: by-domain and the resolver log"),
    selection: z.string().describe("CID selection, e.g. //F or //F/32-64"),
    color: HexColour,
  })
  .strict();

// --- superpose ------------------------------------------------------------

const SuperposeSsm = z
  .object({
    method: z.literal("ssm"),
    move: z.string().describe("file being transformed"),
    onto: z.string().describe("reference file (unchanged)"),
    movChain: z.string(),
    refChain: z.string(),
  })
  .strict();

const LsqMatch = z
  .object({
    refChain: z.string(),
    refRange: ResidueRange,
    movChain: z.string(),
    movRange: ResidueRange,
  })
  .strict();

const SuperposeLsq = z
  .object({
    method: z.literal("lsq"),
    move: z.string(),
    onto: z.string(),
    matches: z.array(LsqMatch).optional(),
    chain: z.string().optional().describe("shorthand chain, both sides"),
    range: ResidueRange.optional().describe("shorthand range, both sides"),
    matchType: z.enum(["all", "main", "ca"]).optional().describe('default "main"'),
  })
  .strict()
  .superRefine((s, ctx) => {
    if (s.matches && (s.chain || s.range)) {
      ctx.addIssue({
        code: z.ZodIssueCode.custom,
        message: "use either matches or chain+range, not both",
      });
    }
  });

export const Superpose = z.discriminatedUnion("method", [
  SuperposeSsm,
  SuperposeLsq,
]);

// --- elements / representations -------------------------------------------

/**
 * Honoured geometry — physical dimensions (Å) a conforming renderer MUST
 * reproduce; ignoring them draws the wrong-sized object (design doc §4).
 * Maps to Moorhen's per-representation m2tParameters. Quality-only knobs
 * (tessellation/sampling) are intentionally excluded — those are render
 * quality, not depicted geometry.
 */
export const Geometry = z
  .object({
    bondRadius: z.number().positive().optional().describe("bond cylinder radius (Å)"),
    ballRadius: z.number().positive().optional().describe("ball-and-stick atom ball radius (Å)"),
    vdwScale: z.number().positive().optional().describe("VdW sphere radius multiplier (×element radius)"),
    probeRadius: z.number().positive().optional().describe("molecular-surface solvent probe radius (Å)"),
    ribbonCoilThickness: z.number().positive().optional().describe("ribbon coil thickness (Å)"),
    ribbonHelixWidth: z.number().positive().optional().describe("ribbon helix width (Å)"),
    ribbonStrandWidth: z.number().positive().optional().describe("ribbon strand width (Å)"),
    ribbonArrowWidth: z.number().positive().optional().describe("ribbon arrow width (Å)"),
    ribbonDNARNAWidth: z.number().positive().optional().describe("nucleotide ribbon width (Å)"),
  })
  .strict();

/**
 * Moorhen representation styles. Mirrors Moorhen's `RepresentationStyles`
 * union (MoorhenMoleculeRepresentation.d.ts). An enum (not a bare string) so
 * the JSON Schema drives editor autocompletion + conformance for the style
 * field. Keep in sync with Moorhen on version bumps (a CI check could assert
 * equality against the installed type).
 */
export const RepresentationStyle = z.enum([
  "VdwSpheres", "ligands", "CAs", "CBs", "CDs", "gaussian", "allHBonds",
  "rama", "rotamer", "CRs", "MolecularSurface", "DishyBases", "VdWSurface",
  "Calpha", "unitCell", "hover", "environment", "ligand_environment",
  "contact_dots", "chemical_features", "ligand_validation", "glycoBlocks",
  "restraints", "residueSelection", "MetaBalls", "adaptativeBonds",
  "StickBases", "residue_environment", "transformation",
]);

export const Representation = z
  .object({
    style: RepresentationStyle.describe('Moorhen RepresentationStyle, e.g. "CRs", "CBs"'),
    selection: z.string().optional().describe("CID; default all-atoms"),
    colour: Colour.optional(),
    alpha: z.number().min(0).max(1).optional().describe("opacity 0..1 (honoured: governs visibility)"),
    geometry: Geometry.optional(),
  })
  .strict();

export const Element = z
  .object({
    file: z.string().describe("name of a files[] entry"),
    dictionaries: z.array(z.string()).optional(),
    colour: Colour.optional().describe(
      "molecule-scoped colour: the default for every representation of this file; a representation's own `colour` overrides it",
    ),
    representations: z.array(Representation).optional(),
  })
  .strict();

// --- maps -----------------------------------------------------------------

export const MapColumns = z
  .object({
    F: z.string().optional(),
    PHI: z.string().optional(),
    Fobs: z.string().optional(),
    SigFobs: z.string().optional(),
    FreeR: z.string().optional(),
    useWeight: z.boolean().optional(),
    calcStructFact: z.boolean().optional(),
  })
  .strict();

export const SceneMapSchema = z
  .object({
    name: z.string(),
    file: z.string().describe("name of a files[] entry (kind mtz or map)"),
    columns: MapColumns.optional().describe("required for mtz, omit for map"),
    isMask: z.boolean().optional(),
    isDifference: z.boolean().optional(),
    contourLevel: z.number().optional().describe("rmsd-relative"),
    radius: z.number().optional().describe("contour radius (Å)"),
    alpha: z.number().min(0).max(1).optional(),
    style: z.enum(["lines", "solid", "lit-lines"]).optional(),
    colour: HexColour.optional().describe("non-difference maps only"),
    positiveColour: HexColour.optional(),
    negativeColour: HexColour.optional(),
    visible: z.boolean().optional(),
  })
  .strict();

// --- view -----------------------------------------------------------------

const Centre = z
  .object({
    file: z.string().optional(),
    selection: z.string().optional(),
  })
  .strict();

const ClipFieldDepth = z
  .object({ front: z.number(), back: z.number() })
  .strict();

const Clip = z.union([z.literal("auto"), z.literal("lock"), ClipFieldDepth]);

const Slab = z
  .object({
    file: z.string().optional(),
    selection: z.string().optional(),
    pad: z.number().optional().describe("extra Å each side"),
  })
  .strict();

export const View = z
  .object({
    origin: z.tuple([z.number(), z.number(), z.number()]).optional(),
    centre: Centre.optional().describe("centroid of a selection; beats origin"),
    quat: z.tuple([z.number(), z.number(), z.number(), z.number()]).optional(),
    zoom: z.number().optional(),
    clipStart: z.number().optional(),
    clipEnd: z.number().optional(),
    fogStart: z.number().optional(),
    fogEnd: z.number().optional(),
    clip: Clip.optional(),
    slab: Slab.optional().describe("z-depth window for a selection; beats clip"),
    background: HexColour.optional(),
  })
  .strict();

// --- provenance / resolver ------------------------------------------------

const Provenance = z
  .object({
    projectId: z.string().optional(),
    projectName: z.string().optional(),
    createdAt: z.string().optional(),
    createdBy: z.string().optional(),
    ccp4i2Version: z.string().optional(),
  })
  .strict();

const ResolverOptions = z
  .object({
    onMissingResidues: z.enum(["clamp-and-log", "strict"]).optional(),
  })
  .strict();

// --- hints (advisory render layer) ----------------------------------------

/**
 * Scene lighting. Mirrors Moorhen's single scene-global light (glRefSlice):
 * a directional light plus ambient/diffuse/specular colours and a specular
 * power. `direction` is the "substituted" class — a renderer that can't honour
 * it falls back to its own default rather than omitting light (design doc §4a).
 */
export const Lighting = z
  .object({
    direction: z
      .tuple([z.number(), z.number(), z.number()])
      .optional()
      .describe("principal directional light vector → Moorhen lightPosition"),
    ambient: HexColour.optional().describe("ambient light colour"),
    diffuse: HexColour.optional().describe("diffuse light colour"),
    specular: HexColour.optional().describe("specular light colour"),
    shininess: z.number().optional().describe("specular power → Moorhen specularPower"),
  })
  .strict();

/** Perceptual post-processing. Additive: each has a meaningful "off"; a
 *  renderer that omits one renders correct-but-plainer (design doc §4a). */
export const Effects = z
  .object({
    ssao: z
      .object({
        enabled: z.boolean().optional(),
        radius: z.number().optional(),
        bias: z.number().optional(),
      })
      .strict()
      .optional()
      .describe("screen-space ambient occlusion"),
    edgeDetect: z
      .object({
        enabled: z.boolean().optional(),
        depthThreshold: z.number().optional(),
        normalThreshold: z.number().optional(),
        depthScale: z.number().optional(),
        normalScale: z.number().optional(),
      })
      .strict()
      .optional(),
    shadows: z.boolean().optional(),
    depthBlur: z
      .object({ radius: z.number().optional(), depth: z.number().optional() })
      .strict()
      .optional()
      .describe("depth-of-field blur"),
    perspective: z.boolean().optional().describe("perspective vs orthographic"),
  })
  .strict();

/** Advisory render hints — a conforming renderer MAY ignore any of these and
 *  still produce a correct image (design doc §3, §4). */
export const Hints = z
  .object({
    lighting: Lighting.optional(),
    effects: Effects.optional(),
  })
  .strict();

// --- top-level scene builder ----------------------------------------------

/**
 * Build the full scene schema over a given file-ref schema, so the portable
 * core and the ccp4i2 dialect share one structure and differ only in which
 * file references they admit. Cross-reference integrity (names that must
 * resolve against files[]/maps[]) is enforced here in superRefine — it is the
 * one class of rule no JSON Schema can express.
 */
export function buildScene<T extends z.ZodTypeAny>(fileRef: T) {
  return z
    .object({
      scene: z.string().describe("human-readable scene name"),
      version: z.number().describe(`must equal ${SCENE_SCHEMA_VERSION}`),
      authoredIn: Provenance.optional(),
      files: z.array(fileRef).optional(),
      superpose: z.array(Superpose).optional(),
      globalDictionaries: z.array(z.string()).optional(),
      domains: z.array(Domain).optional(),
      elements: z.array(Element).optional(),
      maps: z.array(SceneMapSchema).optional(),
      activeMap: z.string().optional(),
      view: View.optional(),
      hints: Hints.optional(),
      resolver: ResolverOptions.optional(),
    })
    .strict()
    .superRefine((s, ctx) => {
      if (s.version !== SCENE_SCHEMA_VERSION) {
        ctx.addIssue({
          code: z.ZodIssueCode.custom,
          path: ["version"],
          message: `unsupported version ${s.version}; this reader handles ${SCENE_SCHEMA_VERSION}`,
        });
      }
      const files = (s.files ?? []) as { name: string; kind?: string }[];
      const fileNames = new Set(files.map((f) => f.name));
      const dictNames = new Set(
        files.filter((f) => f.kind === "dictionary").map((f) => f.name),
      );
      const mapNames = new Set((s.maps ?? []).map((m) => m.name));

      const ref = (
        name: string | undefined,
        set: Set<string>,
        path: (string | number)[],
        what: string,
      ) => {
        if (name != null && !set.has(name)) {
          ctx.addIssue({
            code: z.ZodIssueCode.custom,
            path,
            message: `unknown ${what} "${name}"`,
          });
        }
      };

      (s.elements ?? []).forEach((e, i) => {
        ref(e.file, fileNames, ["elements", i, "file"], "file");
        (e.dictionaries ?? []).forEach((d, j) =>
          ref(d, dictNames, ["elements", i, "dictionaries", j], "dictionary"),
        );
      });
      (s.maps ?? []).forEach((m, i) =>
        ref(m.file, fileNames, ["maps", i, "file"], "file"),
      );
      (s.superpose ?? []).forEach((sp, i) => {
        ref(sp.move, fileNames, ["superpose", i, "move"], "file");
        ref(sp.onto, fileNames, ["superpose", i, "onto"], "file");
      });
      (s.globalDictionaries ?? []).forEach((d, i) =>
        ref(d, dictNames, ["globalDictionaries", i], "dictionary"),
      );
      ref(s.activeMap, mapNames, ["activeMap"], "map");
      if (s.view?.centre?.file)
        ref(s.view.centre.file, fileNames, ["view", "centre", "file"], "file");
      if (s.view?.slab?.file)
        ref(s.view.slab.file, fileNames, ["view", "slab", "file"], "file");
    });
}

/** The portable, upstreamable scene schema. */
export const CoreScene = buildScene(CoreFileRef);
export type CoreScene = z.infer<typeof CoreScene>;
