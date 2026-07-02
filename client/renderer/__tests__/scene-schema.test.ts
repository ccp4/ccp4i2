/**
 * Moorhen Scene Zod schema — drift guard + validity tests.
 *
 * The "published contracts" block regenerates the JSON Schemas from the Zod
 * source and asserts byte-equality with the committed files, so the published
 * moorhen-scene.{core,ccp4i2}.v1.json can never drift from the schema. To
 * (re)generate after an intentional schema change:
 *
 *   UPDATE_SCHEMA=1 npx vitest run renderer/__tests__/scene-schema.test.ts
 */
import { readFileSync, writeFileSync, existsSync, mkdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";

import { describe, it, expect } from "vitest";

import {
  parseScene,
  validateScene,
  buildJsonSchemas,
  buildStructuredJsonSchema,
  serialiseJsonSchema,
  normaliseGeneratedScene,
  SceneParseError,
} from "../lib/scene";
import { sceneMarkers, pathToSegments } from "../lib/scene/yaml-markers";
import { buildSceneBrief, buildSceneMarkdown } from "../lib/scene/brief";
import { buildSceneSystemPrompt } from "../lib/moorhen-scene-prompt";

const SCENE_DIR = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  "../lib/scene",
);
const CONTRACTS = {
  core: path.join(SCENE_DIR, "moorhen-scene.core.v1.json"),
  ccp4i2: path.join(SCENE_DIR, "moorhen-scene.ccp4i2.v1.json"),
} as const;

// Server-side mirror of the artifacts a slim (frontend-free) Django deployment
// consumes — Materia's nlp_scene endpoint reads them from the installed ccp4i2
// package via importlib.resources. client/renderer is canonical (the Zod source
// lives here); this keeps the server copies byte-identical so they can never
// drift. See MOORHEN_SCENES_SCHEMA_V1_DESIGN.md §12.
const SERVER_CONTRACTS_DIR = path.resolve(
  path.dirname(fileURLToPath(import.meta.url)),
  "../../../server/ccp4i2/scene_contracts",
);

/** Assert a server-tree mirror is byte-identical to the canonical content
 *  (regenerating it under UPDATE_SCHEMA). */
function assertServerMirror(filename: string, content: string): void {
  const p = path.join(SERVER_CONTRACTS_DIR, filename);
  if (process.env.UPDATE_SCHEMA) {
    mkdirSync(SERVER_CONTRACTS_DIR, { recursive: true });
    writeFileSync(p, content);
  }
  expect(
    existsSync(p),
    `server mirror ${filename} missing; run with UPDATE_SCHEMA=1`,
  ).toBe(true);
  expect(readFileSync(p, "utf8")).toBe(content);
}

describe("published JSON Schema contracts", () => {
  const built = buildJsonSchemas();

  for (const key of ["core", "ccp4i2"] as const) {
    it(`${key} schema matches the committed contract`, () => {
      const serialised = serialiseJsonSchema(built[key]);
      if (process.env.UPDATE_SCHEMA) {
        writeFileSync(CONTRACTS[key], serialised);
      }
      expect(
        existsSync(CONTRACTS[key]),
        "committed contract missing; run with UPDATE_SCHEMA=1",
      ).toBe(true);
      expect(serialised).toBe(readFileSync(CONTRACTS[key], "utf8"));
    });
  }

  it("types/moorhen-scene.md matches the generated human docs", () => {
    const md = buildSceneMarkdown();
    const mdPath = path.resolve(
      path.dirname(fileURLToPath(import.meta.url)),
      "../types/moorhen-scene.md",
    );
    if (process.env.UPDATE_SCHEMA) writeFileSync(mdPath, md);
    expect(readFileSync(mdPath, "utf8")).toBe(md);
  });

  it("moorhen-scene.system-prompt.v1.md matches the generated system prompt", () => {
    // The static system prompt is committed verbatim so a pinned submodule (e.g.
    // Materia's nlp_scene endpoint) can read it as the server-held system message
    // without importing the TS. This drift guard keeps that file byte-identical
    // to buildSceneSystemPrompt(); regenerate with UPDATE_SCHEMA=1.
    const sys = buildSceneSystemPrompt();
    const sysPath = path.join(SCENE_DIR, "moorhen-scene.system-prompt.v1.md");
    if (process.env.UPDATE_SCHEMA) writeFileSync(sysPath, sys);
    expect(existsSync(sysPath), "committed system prompt missing; run with UPDATE_SCHEMA=1").toBe(true);
    expect(readFileSync(sysPath, "utf8")).toBe(sys);
    assertServerMirror("moorhen-scene.system-prompt.v1.md", sys);
  });

  const STRUCTURED = path.join(SCENE_DIR, "moorhen-scene.structured.v1.json");

  it("moorhen-scene.structured.v1.json matches the strict profile", () => {
    const serialised = serialiseJsonSchema(buildStructuredJsonSchema());
    if (process.env.UPDATE_SCHEMA) writeFileSync(STRUCTURED, serialised);
    expect(existsSync(STRUCTURED), "committed structured profile missing; run with UPDATE_SCHEMA=1").toBe(true);
    expect(serialised).toBe(readFileSync(STRUCTURED, "utf8"));
    assertServerMirror("moorhen-scene.structured.v1.json", serialised);
  });

  it("the strict profile obeys OpenAI Structured Outputs shape rules", () => {
    // Walk every node and assert the invariants strict mode enforces: no dropped
    // vocabulary survives, oneOf is gone, and every object is closed with all
    // keys required. (The lossy constraints are re-checked by parseScene.)
    const FORBIDDEN = [
      "pattern", "format", "minimum", "maximum", "exclusiveMinimum",
      "exclusiveMaximum", "multipleOf", "minLength", "maxLength", "minItems",
      "maxItems", "uniqueItems", "default", "oneOf", "allOf", "not", "$schema",
      "prefixItems", // Azure strict rejects tuple arrays; must be lowered to items
    ];
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const walk = (node: any): void => {
      if (Array.isArray(node)) return node.forEach(walk);
      if (!node || typeof node !== "object") return;
      for (const key of FORBIDDEN) {
        expect(key in node, `forbidden keyword "${key}" leaked into the strict profile`).toBe(false);
      }
      const isObjectNode =
        node.properties !== undefined ||
        node.type === "object" ||
        (Array.isArray(node.type) && node.type.includes("object"));
      if (isObjectNode && node.properties) {
        expect(node.additionalProperties).toBe(false);
        expect([...(node.required ?? [])].sort()).toEqual(
          Object.keys(node.properties).sort(),
        );
      }
      const isArrayNode =
        node.type === "array" ||
        (Array.isArray(node.type) && node.type.includes("array"));
      if (isArrayNode) {
        // Every array must carry `items` (no bare/tuple arrays) — Azure strict.
        expect(node.items, "array node missing items").toBeDefined();
      }
      for (const v of Object.values(node)) walk(v);
    };
    walk(buildStructuredJsonSchema());
  });

  it("the strict profile stays under Azure's 100-property cap (with headroom)", () => {
    // Azure OpenAI hard-caps a strict json_schema at 100 object properties total
    // and 5 nesting levels. The authoring-core prune (STRUCTURED_PRUNE) keeps us
    // under; this guard fails loudly if a schema addition pushes it back over,
    // rather than the model 400ing in production.
    let properties = 0;
    let maxDepth = 0;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const walk = (node: any, depth: number): void => {
      if (Array.isArray(node)) return node.forEach((n) => walk(n, depth));
      if (!node || typeof node !== "object") return;
      const here = node.properties ? depth + 1 : depth;
      if (node.properties) {
        properties += Object.keys(node.properties).length;
        maxDepth = Math.max(maxDepth, here);
      }
      for (const v of Object.values(node)) walk(v, here);
    };
    walk(buildStructuredJsonSchema(), 0);
    expect(properties).toBeLessThanOrEqual(90); // ≤100 cap, ≥10 headroom
    expect(maxDepth).toBeLessThanOrEqual(5);
  });
});

// --- a minimal portable scene used across validity tests ------------------
const VALID_CORE = {
  scene: "demo",
  version: 1,
  files: [{ name: "prot", pdb: "1abc" }],
  elements: [
    {
      file: "prot",
      representations: [{ style: "CRs", colour: "by-domain", alpha: 0.8 }],
    },
  ],
  view: { centre: { file: "prot", selection: "//A" } },
};

describe("validity", () => {
  it("accepts a minimal portable scene", () => {
    expect(validateScene(VALID_CORE).errors).toEqual([]);
  });

  it("rejects unknown keys (strict)", () => {
    const { errors } = validateScene({ ...VALID_CORE, wibble: 1 });
    expect(errors.length).toBeGreaterThan(0);
  });

  it("flags a dangling activeMap cross-reference", () => {
    const { errors } = validateScene({ ...VALID_CORE, activeMap: "nope" });
    expect(errors.some((e) => e.path === "activeMap")).toBe(true);
  });

  it("flags an element referencing an unknown file", () => {
    const { errors } = validateScene({
      ...VALID_CORE,
      elements: [{ file: "ghost" }],
    });
    expect(errors.some((e) => e.path.includes("elements"))).toBe(true);
  });

  it("requires exactly one file source", () => {
    expect(
      validateScene({ ...VALID_CORE, files: [{ name: "x" }] }).errors.length,
    ).toBeGreaterThan(0);
    expect(
      validateScene({
        ...VALID_CORE,
        files: [{ name: "x", pdb: "1abc", url: "https://e/x.cif" }],
      }).errors.length,
    ).toBeGreaterThan(0);
  });

  it("validates the superpose discriminated union", () => {
    const ssm = {
      scene: "sup",
      version: 1,
      files: [
        { name: "a", pdb: "1abc" },
        { name: "b", pdb: "2def" },
      ],
      superpose: [
        { method: "ssm", move: "b", onto: "a", movChain: "A", refChain: "A" },
      ],
    };
    expect(validateScene(ssm).errors).toEqual([]);
    // unknown move file → cross-ref error
    const bad = { ...ssm, superpose: [{ ...ssm.superpose[0], move: "zzz" }] };
    expect(validateScene(bad).errors.some((e) => e.path.includes("move"))).toBe(
      true,
    );
  });
});

describe("hints + honoured geometry", () => {
  it("accepts scene-level lighting/effects hints", () => {
    const s = {
      ...VALID_CORE,
      hints: {
        lighting: { direction: [1, 1, 1], ambient: "#202020", shininess: 32 },
        effects: { ssao: { enabled: true, radius: 0.5 }, edgeDetect: { enabled: true } },
      },
    };
    expect(validateScene(s).errors).toEqual([]);
  });

  it("accepts per-representation honoured geometry (Å)", () => {
    const s = {
      ...VALID_CORE,
      elements: [
        {
          file: "prot",
          representations: [
            { style: "CBs", geometry: { bondRadius: 0.2, ballRadius: 0.4 } },
            { style: "CRs", geometry: { ribbonHelixWidth: 1.4 } },
          ],
        },
      ],
    };
    expect(validateScene(s).errors).toEqual([]);
  });

  it("rejects unknown keys in hints and geometry (strict)", () => {
    expect(
      validateScene({ ...VALID_CORE, hints: { lighting: { glow: 1 } } }).errors
        .length,
    ).toBeGreaterThan(0);
    expect(
      validateScene({
        ...VALID_CORE,
        elements: [
          { file: "prot", representations: [{ style: "CBs", geometry: { width: 2 } }] },
        ],
      }).errors.length,
    ).toBeGreaterThan(0);
  });

  it("rejects a non-positive bond radius", () => {
    const { errors } = validateScene({
      ...VALID_CORE,
      elements: [
        { file: "prot", representations: [{ style: "CBs", geometry: { bondRadius: 0 } }] },
      ],
    });
    expect(errors.length).toBeGreaterThan(0);
  });
});

describe("profiles (portable vs permissive)", () => {
  const sceneWithFileId = {
    scene: "p",
    version: 1,
    files: [{ name: "x", fileId: 7, projectId: "uuid-1" }],
    elements: [{ file: "x" }],
  };

  it("permissive accepts ccp4i2 deployment refs", () => {
    expect(validateScene(sceneWithFileId).errors).toEqual([]);
  });

  it("strict-portable rejects ccp4i2 deployment refs", () => {
    const { errors } = validateScene(sceneWithFileId, { portable: true });
    expect(errors.length).toBeGreaterThan(0);
  });
});

describe("real authored scene (fixture-corpus seed)", () => {
  it("the demo fixture validates end-to-end through the live contract", () => {
    const TEST_DIR = path.dirname(fileURLToPath(import.meta.url));
    const yaml = readFileSync(
      path.join(TEST_DIR, "fixtures/demo.scene.yaml"),
      "utf8",
    );
    const scene = parseScene(yaml);
    expect(scene.scene).toBe("gamma-demo");
    // the step-2 additions survive the parse (not stripped)
    expect(scene.hints?.lighting?.shininess).toBe(24);
    expect(scene.elements?.[0].representations?.[1].geometry?.bondRadius).toBe(
      0.18,
    );
    // pdb: ref ⇒ portable: also validates under the upstreamable core profile
    expect(() => parseScene(yaml, { portable: true })).not.toThrow();
  });
});

describe("editor markers (in-app conformance squiggles)", () => {
  const ERR = 8; // monaco.MarkerSeverity.Error

  it("parses dotted+bracket paths to segments", () => {
    expect(pathToSegments("elements[0].representations[1].style")).toEqual([
      "elements",
      0,
      "representations",
      1,
      "style",
    ]);
  });

  it("no markers for a valid scene", () => {
    const yaml = "scene: x\nversion: 1\nfiles:\n  - { name: a, pdb: 1abc }\n";
    expect(sceneMarkers(yaml, ERR)).toEqual([]);
  });

  it("locates a bad representation style at its source line", () => {
    const yaml = [
      "scene: x",
      "version: 1",
      "files:",
      "  - { name: a, pdb: 1abc }",
      "elements:",
      "  - file: a",
      "    representations:",
      "      - { style: NOPE }",
    ].join("\n");
    const markers = sceneMarkers(yaml, ERR);
    expect(markers.length).toBeGreaterThan(0);
    // the offending style is on line 8
    expect(markers.some((m) => m.startLineNumber === 8 && m.severity === ERR)).toBe(
      true,
    );
  });
});

describe("LLM authoring brief (generated from the contract)", () => {
  const brief = buildSceneBrief();

  it("is compact (≪ the 38KB prose grammar) and prints its size", () => {
    const chars = brief.length;
    const approxTokens = Math.round(chars / 4);
    // eslint-disable-next-line no-console
    console.log(
      `\n[scene brief] ${chars} chars ≈ ${approxTokens} tokens\n` +
        "──────────────────────────────────────────\n" +
        brief +
        "\n──────────────────────────────────────────",
    );
    expect(chars).toBeLessThan(8000); // vs ~38000 for the .md
  });

  it("teaches the CURRENT format (not the stale grammar)", () => {
    expect(brief).toContain("relativeUrl");
    expect(brief).toContain("geometry");
    expect(brief).toContain("hints");
    expect(brief).toContain("edgeDetect");
    expect(brief).not.toContain('"outline"'); // removed effect
  });

  it("enumerates valid representation styles + colour schemes authoritatively", () => {
    expect(brief).toContain('"CRs"');
    expect(brief).toContain('"MolecularSurface"');
    expect(brief).toContain('"by-domain"');
  });
});

describe("parseScene", () => {
  it("parses YAML and strips a code fence", () => {
    const yaml = "```yaml\nscene: x\nversion: 1\n```";
    expect(parseScene(yaml)).toMatchObject({ scene: "x", version: 1 });
  });

  it("throws SceneParseError on invalid scene", () => {
    expect(() => parseScene("scene: x\nversion: 99")).toThrow(SceneParseError);
  });
});

describe("normaliseGeneratedScene (LLM ingest tidy-up)", () => {
  it("strips explicit nulls that strict Structured Outputs emits", () => {
    // Strict mode makes every optional present-as-null; Zod optional ≠ null,
    // so these must be dropped before parseScene sees them.
    const json = JSON.stringify({
      scene: "s",
      version: 1,
      files: [{ name: "prot", pdb: "1abc" }],
      view: { centre: { file: "prot", selection: "//A" }, zoom: null, origin: null },
      hints: null,
    });
    const out = normaliseGeneratedScene(json);
    expect(out).not.toContain("null");
    const parsed = parseScene(out);
    expect(parsed).toMatchObject({ scene: "s", version: 1 });
    expect(parsed.view?.zoom).toBeUndefined();
  });

  it("strips a ```json fence and re-serialises to YAML", () => {
    const fenced = "```json\n{\"scene\":\"s\",\"version\":1}\n```";
    const out = normaliseGeneratedScene(fenced);
    expect(out).not.toContain("```");
    expect(parseScene(out)).toMatchObject({ scene: "s", version: 1 });
  });

  it("leaves plain YAML without nulls essentially unchanged in meaning", () => {
    const yaml = "scene: s\nversion: 1\n";
    expect(parseScene(normaliseGeneratedScene(yaml))).toMatchObject({
      scene: "s",
      version: 1,
    });
  });

  it("falls back to the raw text when it doesn't parse", () => {
    const junk = "this is: not: valid: yaml: [";
    expect(normaliseGeneratedScene(junk)).toBe(junk);
  });

  it("keeps array elements (only object nulls are dropped)", () => {
    const json = JSON.stringify({
      scene: "s",
      version: 1,
      domains: [{ name: "d", selection: "//A", color: "#112233" }],
    });
    const parsed = parseScene(normaliseGeneratedScene(json));
    expect(parsed.domains).toHaveLength(1);
  });
});
