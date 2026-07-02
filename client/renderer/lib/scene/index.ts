/**
 * Moorhen Scene — public API over the Zod schema.
 *
 * parse / validate / serialise, plus the JSON-Schema generator that emits the
 * published contracts. Two profiles (design doc §2):
 *
 *   - permissive (default)  → ccp4i2 dialect: deployment refs allowed.
 *   - strict-portable       → core schema: only portable refs; the same check
 *                             that gates export/publish/share.
 *
 * This is additive: the legacy hand validator in ../moorhen-scene.ts is still
 * live and unchanged. Wiring the lifter/resolver onto this schema is a later
 * step (design doc §3 order).
 */
import * as YAML from "yaml";
import { z } from "zod";

import type { MoorhenScene } from "../../types/moorhen-scene";
import { CoreScene } from "./core";
import { Ccp4i2Scene } from "./dialect";

export { CoreScene } from "./core";
export { Ccp4i2Scene } from "./dialect";

export interface SceneValidationError {
  path: string;
  message: string;
}

export interface ParseOptions {
  /** strict-portable profile: reject deployment-coupled refs. */
  portable?: boolean;
}

/** Strip a Markdown ```yaml fence (LLM output); plain YAML is unaffected. */
function stripCodeFence(text: string): string {
  const m = text.match(/```[ \t]*[\w+-]*[ \t]*\r?\n([\s\S]*?)\r?\n?```/);
  return m ? m[1] : text;
}

function schemaFor(opts?: ParseOptions) {
  return opts?.portable ? CoreScene : Ccp4i2Scene;
}

/** Render a Zod path as files[1].fileId style. */
function fmtPath(path: PropertyKey[]): string {
  return path
    .map((p) => (typeof p === "number" ? `[${p}]` : `.${String(p)}`))
    .join("")
    .replace(/^\./, "");
}

function toErrors(err: z.ZodError): SceneValidationError[] {
  return err.issues.map((i) => ({ path: fmtPath(i.path), message: i.message }));
}

/** Validate a parsed JS value. Returns scene-or-null plus errors.
 *  Typed as MoorhenScene for the app boundary; the Zod-inferred shape is a
 *  structural subset (no deprecated domain chain/range), so the cast is sound. */
export function validateScene(
  raw: unknown,
  opts?: ParseOptions,
): { scene: MoorhenScene | null; errors: SceneValidationError[] } {
  const result = schemaFor(opts).safeParse(raw);
  if (result.success) return { scene: result.data as MoorhenScene, errors: [] };
  return { scene: null, errors: toErrors(result.error) };
}

export class SceneParseError extends Error {
  constructor(public readonly errors: SceneValidationError[]) {
    super(
      `Invalid Moorhen scene: ${errors.length} error(s)\n` +
        errors.map((e) => `  - ${e.path}: ${e.message}`).join("\n"),
    );
    this.name = "SceneParseError";
  }
}

/** Parse a YAML scene document; throws SceneParseError on invalid input. */
export function parseScene(yamlText: string, opts?: ParseOptions): MoorhenScene {
  let raw: unknown;
  try {
    raw = YAML.parse(stripCodeFence(yamlText));
  } catch (e) {
    throw new SceneParseError([
      { path: "", message: `YAML parse error: ${(e as Error).message}` },
    ]);
  }
  const { scene, errors } = validateScene(raw, opts);
  if (errors.length > 0) throw new SceneParseError(errors);
  return scene as MoorhenScene;
}

/** Serialise a scene to YAML. */
export function serialiseScene(scene: unknown): string {
  return YAML.stringify(scene);
}

/** Recursively drop null-valued object keys (keeps array elements). */
function stripNulls(v: unknown): unknown {
  if (Array.isArray(v)) return v.map(stripNulls);
  if (v && typeof v === "object") {
    const out: Record<string, unknown> = {};
    for (const [k, val] of Object.entries(v as Record<string, unknown>)) {
      if (val === null) continue;
      out[k] = stripNulls(val);
    }
    return out;
  }
  return v;
}

/**
 * Prepare raw LLM output for the editor. Strips a ```yaml/```json fence, then —
 * because strict OpenAI Structured Outputs emits explicit `null`s for absent
 * optionals and Zod treats *optional* ≠ *null* — drops those nulls, and
 * re-serialises to tidy YAML (JSON ⊂ YAML, so JSON output normalises too). Falls
 * back to the raw text if it doesn't parse, so the editor still shows something
 * the validity markers can flag. Parsing/validation stay with parseScene; this
 * is only ingest tid-up.
 */
export function normaliseGeneratedScene(text: string): string {
  try {
    const raw = YAML.parse(stripCodeFence(text));
    if (raw == null || typeof raw !== "object") return text;
    return YAML.stringify(stripNulls(raw));
  } catch {
    return text;
  }
}

// --- JSON Schema generation (the published contracts) ---------------------

/**
 * Build the published JSON Schemas from the Zod source. Returns both the
 * upstreamable core and the ccp4i2 dialect. Cross-reference rules (superRefine)
 * are not representable in JSON Schema and are intentionally absent here — they
 * live only in the runtime validator.
 */
export function buildJsonSchemas(): { core: unknown; ccp4i2: unknown } {
  const opts = { target: "draft-2020-12" as const };
  return {
    core: z.toJSONSchema(CoreScene, opts),
    ccp4i2: z.toJSONSchema(Ccp4i2Scene, opts),
  };
}

/** Canonical serialisation used for both the committed file and the drift test. */
export function serialiseJsonSchema(schema: unknown): string {
  return JSON.stringify(schema, null, 2) + "\n";
}

// --- Strict OpenAI Structured Outputs profile -----------------------------

/**
 * OpenAI's `json_schema` strict mode rejects most JSON-Schema validation
 * vocabulary and imposes its own shape rules. This is the closed set we keep on
 * any node; everything else (pattern, format, numeric/length/array bounds,
 * default, $schema, …) is dropped. The runtime Zod validator still enforces
 * those — the model's job here is to emit the right *shape*; parseScene
 * re-checks the constraints (and the cross-references JSON Schema can't express)
 * on the way back in.
 */
const STRUCTURED_KEEP = new Set([
  "type",
  "enum",
  "const",
  "properties",
  "required",
  "items",
  "anyOf",
  "additionalProperties",
  "description",
]);

/* eslint-disable @typescript-eslint/no-explicit-any */
type JsonNode = Record<string, any>;

function isPlainObject(v: unknown): v is JsonNode {
  return typeof v === "object" && v !== null && !Array.isArray(v);
}

/** Make an already-strictified node accept `null`, so an optional field can be
 *  present-as-null (strict mode requires every property in `required`). */
function makeNullable(node: JsonNode): void {
  if (Array.isArray(node.anyOf)) {
    if (!node.anyOf.some((b: JsonNode) => b && b.type === "null")) {
      node.anyOf.push({ type: "null" });
    }
    return;
  }
  if (node.const !== undefined) {
    const c = node.const;
    delete node.const;
    node.anyOf = [{ const: c }, { type: "null" }];
    return;
  }
  if (Array.isArray(node.type)) {
    if (!node.type.includes("null")) node.type.push("null");
    if (node.enum && !node.enum.includes(null)) node.enum.push(null);
    return;
  }
  if (typeof node.type === "string") {
    node.type = [node.type, "null"];
    if (node.enum && !node.enum.includes(null)) node.enum.push(null);
    return;
  }
  // Enum-only / shapeless node: wrap so null is representable.
  const inner = { ...node };
  for (const k of Object.keys(node)) delete node[k];
  node.anyOf = [inner, { type: "null" }];
}

/** Recursively rewrite a draft-2020-12 schema node into OpenAI strict form:
 *  drop unsupported keywords, oneOf→anyOf, every object gets
 *  additionalProperties:false + all keys required, optionals become nullable. */
function strictify(input: unknown): JsonNode {
  if (!isPlainObject(input)) return input as JsonNode;

  const node: JsonNode = {};
  for (const [k, v] of Object.entries(input)) {
    if (k === "oneOf") {
      node.anyOf = v;
    } else if (STRUCTURED_KEEP.has(k)) {
      node[k] = v;
    }
    // else: dropped (pattern, format, minimum, …)
  }

  if (Array.isArray(node.anyOf)) node.anyOf = node.anyOf.map(strictify);

  // Azure strict mode rejects tuple arrays (`prefixItems` / array-without-`items`).
  // Lower a fixed-length tuple to a homogeneous array: `items` = the common branch
  // schema, or an `anyOf` of the distinct branch types when heterogeneous. The real
  // tuple length (3-vector origin, 4-vector quat) is re-checked by validateScene —
  // same trade as dropping pattern/format for strict.
  const inputPrefix = (input as JsonNode).prefixItems;
  if (Array.isArray(inputPrefix)) {
    const branches = inputPrefix.map(strictify);
    const uniq: JsonNode[] = [];
    for (const b of branches) {
      if (!uniq.some((u) => JSON.stringify(u) === JSON.stringify(b))) uniq.push(b);
    }
    node.items = uniq.length === 1 ? uniq[0] : { anyOf: uniq };
  } else if (node.items !== undefined) {
    node.items = Array.isArray(node.items) ? node.items.map(strictify) : strictify(node.items);
  }

  if (isPlainObject(node.properties)) {
    const originalRequired = new Set<string>(
      Array.isArray(node.required) ? node.required : [],
    );
    const keys = Object.keys(node.properties);
    const props: JsonNode = {};
    for (const key of keys) {
      const child = strictify(node.properties[key]);
      if (!originalRequired.has(key)) makeNullable(child);
      props[key] = child;
    }
    node.properties = props;
    // Strict mode: ALL properties must be required (optionals are nullable).
    node.required = keys;
    node.additionalProperties = false;
  }

  return node;
}
/* eslint-enable @typescript-eslint/no-explicit-any */

/**
 * Authoring-core prune list: fields removed from the STRICT generation profile
 * (not from the validator). Azure OpenAI caps a strict `json_schema` at **100
 * object properties total**; the full ccp4i2 contract is ~142, so the strict
 * profile is narrowed to the "authoring core" the design doc §7a anticipated —
 * dropping advisory, escape-hatch, provenance and rarely-authored fields. The
 * runtime `validateScene` still accepts every one of them; this only narrows
 * what the constrained model can *emit*. Each entry is justified:
 */
const STRUCTURED_PRUNE = {
  /** Whole top-level blocks the model should not author. */
  topLevel: [
    "authoredIn", // provenance — stamped by the app, never by the model
    "hints", // advisory lighting/effects; the resolver doesn't even apply lighting yet
  ],
  /** Property names removed wherever they occur. Each is unique to one block in
   *  this schema, so name-matching is unambiguous. */
  props: new Set([
    "relativeUrl", // deployment ref the system prompt already forbids
    "bundle", // rare file source (.scene.zip asset)
    "cifText", // rare file source (inline dictionary CIF)
    "ribbonCoilThickness", "ribbonHelixWidth", "ribbonStrandWidth",
    "ribbonArrowWidth", "ribbonDNARNAWidth", // keep bond/ball/vdw/probe radii; drop ribbon dims
    "clipStart", "clipEnd", "fogStart", "fogEnd", // slab/clip cover the common case
    "columns", // MTZ column spec — the resolver derives it from the file
  ]),
} as const;

/* eslint-disable @typescript-eslint/no-explicit-any */
/** A colour-union branch that is the raw-rule escape hatch (`{ raw: {…} }`). */
function isRawColourBranch(node: any): boolean {
  return (
    isPlainObject(node) &&
    isPlainObject(node.properties) &&
    Object.keys(node.properties).length === 1 &&
    "raw" in node.properties
  );
}

/** Deep-clone the schema while removing the authoring-core exclusions: top-level
 *  blocks, named properties, and the raw-colour escape-hatch union branch. Run
 *  BEFORE strictify so it recomputes required/additionalProperties cleanly. */
function pruneForAuthoringCore(input: unknown, atRoot = false): any {
  if (Array.isArray(input)) return input.map((n) => pruneForAuthoringCore(n));
  if (!isPlainObject(input)) return input;

  const out: JsonNode = {};
  for (const [k, v] of Object.entries(input)) {
    if (k === "properties" && isPlainObject(v)) {
      const props: JsonNode = {};
      for (const [pk, pv] of Object.entries(v)) {
        if (STRUCTURED_PRUNE.props.has(pk)) continue;
        if (atRoot && (STRUCTURED_PRUNE.topLevel as readonly string[]).includes(pk)) continue;
        props[pk] = pruneForAuthoringCore(pv);
      }
      out[k] = props;
    } else if (k === "required" && Array.isArray(v)) {
      out[k] = v.filter(
        (r: string) =>
          !STRUCTURED_PRUNE.props.has(r) &&
          !(atRoot && (STRUCTURED_PRUNE.topLevel as readonly string[]).includes(r)),
      );
    } else if (k === "anyOf" && Array.isArray(v)) {
      out[k] = v.filter((b) => !isRawColourBranch(b)).map((n) => pruneForAuthoringCore(n));
    } else {
      out[k] = pruneForAuthoringCore(v);
    }
  }
  return out;
}
/* eslint-enable @typescript-eslint/no-explicit-any */

/**
 * Build the strict OpenAI Structured Outputs profile from the ccp4i2 contract.
 *
 * Derived from the ccp4i2 dialect (not core): the production endpoint sends the
 * project manifest as grounding, so the model must be able to emit
 * job/param/fileId refs. Two transforms, both lossy by design — the dropped
 * information is re-checked by parseScene/validateScene, the sole authority:
 *   1. pruneForAuthoringCore — narrows to the authoring core to fit Azure's
 *      100-property strict cap (see STRUCTURED_PRUNE).
 *   2. strictify — OpenAI strict shape (oneOf→anyOf, all-required+nullable,
 *      additionalProperties:false, unsupported vocabulary dropped).
 */
export function buildStructuredJsonSchema(): unknown {
  const { ccp4i2 } = buildJsonSchemas();
  return strictify(pruneForAuthoringCore(ccp4i2, true));
}
