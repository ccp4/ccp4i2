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
