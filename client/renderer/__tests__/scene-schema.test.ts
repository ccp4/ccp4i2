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
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";

import { describe, it, expect } from "vitest";

import {
  parseScene,
  validateScene,
  buildJsonSchemas,
  serialiseJsonSchema,
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
