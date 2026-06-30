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
  "prefixItems",
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
  if (Array.isArray(node.prefixItems)) node.prefixItems = node.prefixItems.map(strictify);
  if (node.items !== undefined) {
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
 * Build the strict OpenAI Structured Outputs profile from the ccp4i2 contract.
 *
 * Derived from the ccp4i2 dialect (not core): the production endpoint sends the
 * project manifest as grounding, so the model must be able to emit
 * job/param/fileId refs. The transform is mechanical and lossy by design —
 * constraints the model can't be trusted to honour structurally (CID patterns,
 * hex format, Å bounds, exactly-one-file-source) are dropped here and re-checked
 * by parseScene/validateScene, which remain the sole authority.
 */
export function buildStructuredJsonSchema(): unknown {
  const { ccp4i2 } = buildJsonSchemas();
  return strictify(ccp4i2);
}
