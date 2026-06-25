/**
 * Lifter tests: Moorhen state → MoorhenScene → YAML → MoorhenScene.
 *
 * We construct minimal stand-in objects for moorhen.Molecule and
 * moorhen.MoleculeRepresentation rather than mocking the full runtime —
 * the lifter only reads a handful of fields, and using duck-typed plain
 * objects keeps the tests focused on the lift logic itself.
 */
import { describe, it, expect } from "vitest";
import type { moorhen } from "moorhen/types/moorhen";

import { liftScene, liftSceneWithHints } from "../lib/moorhen-scene-lifter";
import {
  parseScene,
  serialiseScene,
  serialiseSceneWithComments,
} from "../lib/moorhen-scene";

// Build a fake-but-valid molecule by hand. Cast to moorhen.Molecule via
// unknown — the lifter only reads `name`, `molNo`, `uniqueId`, and
// `representations`, so we don't need a real instance.
function fakeMol(opts: {
  name: string;
  molNo: number;
  uniqueId: string;
  representations?: Partial<moorhen.MoleculeRepresentation>[];
  /** Ligand info as Moorhen would surface it via mol.ligands. */
  ligands?: Partial<moorhen.LigandInfo>[];
  /** Map from comp_id → dict text, mimicking Moorhen's mol.getDict. */
  dicts?: Record<string, string>;
}): moorhen.Molecule {
  const dicts = opts.dicts ?? {};
  return {
    name: opts.name,
    molNo: opts.molNo,
    uniqueId: opts.uniqueId,
    representations: (opts.representations ?? []) as moorhen.MoleculeRepresentation[],
    ligands: (opts.ligands ?? []) as moorhen.LigandInfo[],
    getDict: (comp_id: string) => dicts[comp_id] ?? "",
  } as unknown as moorhen.Molecule;
}

const fakeGlRef = {
  origin: [10.0, 20.0, 30.0],
  quat: [0, 0, 0, -1],
  zoom: 1.5,
};

describe("liftScene", () => {
  it("captures camera and a single file with a single visible representation", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "thing",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/287/download/",
          representations: [
            { style: "CRs", cid: "/*/*/*/*", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
      projectId: "proj-uuid",
      projectName: "kinase-2026",
      sceneName: "test",
      createdBy: "user@example.com",
    });

    expect(scene.scene).toBe("test");
    expect(scene.authoredIn?.projectId).toBe("proj-uuid");
    expect(scene.authoredIn?.createdBy).toBe("user@example.com");
    expect(scene.files).toHaveLength(1);
    expect(scene.files![0]).toMatchObject({
      name: "thing",
      projectId: "proj-uuid",
      projectName: "kinase-2026",
      fileId: 287,
    });
    expect(scene.view?.origin).toEqual([10, 20, 30]);
    expect(scene.view?.zoom).toBe(1.5);
    expect(scene.elements).toHaveLength(1);
    expect(scene.elements![0].representations![0].style).toBe("CRs");
    // CID == default wildcard is omitted from the serialised form.
    expect(scene.elements![0].representations![0].selection).toBeUndefined();
  });

  it("falls back to url: when the molecule wasn't loaded via ccp4i2", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "pdb",
          molNo: 0,
          uniqueId: "https://www.ebi.ac.uk/pdbe/entry-files/download/3jbt.cif",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.files![0]).toEqual({
      name: "pdb",
      url: "https://www.ebi.ac.uk/pdbe/entry-files/download/3jbt.cif",
    });
  });

  it("hides representations that are explicitly hidden", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
            { style: "MolecularSurface", visible: false, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations).toHaveLength(1);
    expect(scene.elements![0].representations![0].style).toBe("CRs");
  });

  it("lifts a named multi-rule colour scheme", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            {
              style: "CRs",
              visible: true,
              colourRules: [{ ruleType: "b-factor", isMultiColourRule: true, args: [] } as unknown as moorhen.ColourRule],
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations![0].colour).toBe("b-factor");
  });

  it("lifts a single-colour hex rule", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            {
              style: "ligands",
              visible: true,
              colourRules: [
                {
                  ruleType: "molecule",
                  isMultiColourRule: false,
                  color: "#2ecc71",
                  args: [],
                } as unknown as moorhen.ColourRule,
              ],
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations![0].colour).toBe("#2ecc71");
  });

  it("captures a non-default nonCustomOpacity as alpha", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            {
              style: "MolecularSurface",
              visible: true,
              colourRules: [],
              nonCustomOpacity: 0.4,
            } as Partial<moorhen.MoleculeRepresentation>,
            // Fully opaque (default) → alpha omitted.
            {
              style: "CRs",
              visible: true,
              colourRules: [],
              nonCustomOpacity: 1,
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations![0].alpha).toBe(0.4);
    expect(scene.elements![0].representations![1].alpha).toBeUndefined();
  });

  it("recognises by-domain pipe-delimited args", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            {
              style: "CRs",
              visible: true,
              colourRules: [
                {
                  ruleType: "by-domain",
                  isMultiColourRule: true,
                  args: ["//A/1-120^#4b8bbe|//A/121-130^#f1c40f|//A/131-300^#e74c3c"],
                } as unknown as moorhen.ColourRule,
              ],
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations![0].colour).toBe("by-domain");
  });

  it("falls back to the raw escape hatch for unrecognised rules", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [
            {
              style: "CRs",
              visible: true,
              colourRules: [
                {
                  ruleType: "some-bespoke-thing",
                  isMultiColourRule: true,
                  args: ["weird", 42],
                  applyColourToNonCarbonAtoms: true,
                } as unknown as moorhen.ColourRule,
              ],
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
    });
    expect(scene.elements![0].representations![0].colour).toEqual({
      raw: {
        ruleType: "some-bespoke-thing",
        args: ["weird", 42],
        isMultiColourRule: true,
        applyColourToNonCarbonAtoms: true,
      },
    });
  });

  it("round-trips: lift → serialise → parse produces an equivalent scene", () => {
    const lifted = liftScene({
      molecules: [
        fakeMol({
          name: "thing",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/42/download/",
          representations: [
            {
              style: "CRs",
              cid: "//A",
              visible: true,
              colourRules: [{ ruleType: "b-factor", isMultiColourRule: true, args: [] } as unknown as moorhen.ColourRule],
            } as Partial<moorhen.MoleculeRepresentation>,
          ],
        }),
      ],
      glRef: fakeGlRef,
      projectId: "p",
      sceneName: "rt",
    });
    const yaml = serialiseScene(lifted);
    const parsed = parseScene(yaml);
    expect(parsed.scene).toBe("rt");
    expect(parsed.files![0]).toMatchObject({ name: "thing", fileId: 42, projectId: "p" });
    expect(parsed.elements![0].representations![0]).toEqual({
      style: "CRs",
      selection: "//A",
      colour: "b-factor",
    });
  });
});

describe("serialiseSceneWithComments", () => {
  it("attaches a per-file comment above the corresponding entry", () => {
    const { scene, hints } = liftSceneWithHints({
      molecules: [
        fakeMol({
          name: "thing",
          molNo: 0,
          uniqueId: "https://example.org/foo.pdb",
          representations: [{ style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>],
        }),
      ],
      glRef: fakeGlRef,
    });
    const yaml = serialiseSceneWithComments(scene, hints.fileComments);
    expect(yaml).toContain("# loaded from https://example.org/foo.pdb");
    // And it still parses cleanly.
    expect(() => parseScene(yaml)).not.toThrow();
  });

  it("falls back to plain serialisation when there are no comments", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "m",
          molNo: 0,
          uniqueId: "x",
          representations: [{ style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>],
        }),
      ],
      glRef: fakeGlRef,
    });
    const withComments = serialiseSceneWithComments(scene, {});
    const plain = serialiseScene(scene);
    expect(withComments).toBe(plain);
  });
});

describe("liftScene — dictionary lifting", () => {
  it("emits a kind: dictionary file ref for each non-standard ligand", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "complex",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/482/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [
            { resName: "LIG", chainName: "A", resNum: "401" } as Partial<moorhen.LigandInfo>,
          ],
          dicts: {
            LIG: "data_comp_LIG\n_chem_comp.id LIG\n",
          },
        }),
      ],
      glRef: fakeGlRef,
      projectId: "p-uuid",
    });

    expect(scene.files).toHaveLength(2); // complex + dict
    const dictRef = scene.files!.find((f) => f.kind === "dictionary");
    expect(dictRef).toBeDefined();
    expect(dictRef!.cifText).toContain("data_comp_LIG");
    // Element should reference the dict by name.
    expect(scene.elements![0].dictionaries).toEqual([dictRef!.name]);
  });

  it("skips standard amino acids and water", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "protein",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/100/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [
            { resName: "HOH" } as Partial<moorhen.LigandInfo>,
            { resName: "ALA" } as Partial<moorhen.LigandInfo>,
            { resName: "GLU" } as Partial<moorhen.LigandInfo>,
          ],
          dicts: {
            HOH: "data_comp_HOH\n",
            ALA: "data_comp_ALA\n",
            GLU: "data_comp_GLU\n",
          },
        }),
      ],
      glRef: fakeGlRef,
    });
    // Only one file (the coords), no dict refs.
    expect(scene.files).toHaveLength(1);
    expect(scene.files![0].kind).toBeUndefined();
  });

  it("dedupes by comp_id across multiple ligand instances", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "complex",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/100/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [
            { resName: "LIG", resNum: "401" } as Partial<moorhen.LigandInfo>,
            { resName: "LIG", resNum: "402" } as Partial<moorhen.LigandInfo>,
            { resName: "LIG", resNum: "403" } as Partial<moorhen.LigandInfo>,
          ],
          dicts: { LIG: "data_comp_LIG\n" },
        }),
      ],
      glRef: fakeGlRef,
    });
    const dictRefs = scene.files!.filter((f) => f.kind === "dictionary");
    expect(dictRefs).toHaveLength(1);
  });

  it("emits separate dict refs for two molecules with same-named ligands", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "complex-A",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/100/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [{ resName: "LIG" } as Partial<moorhen.LigandInfo>],
          dicts: { LIG: "data_comp_LIG\nA chemistry\n" },
        }),
        fakeMol({
          name: "complex-B",
          molNo: 1,
          uniqueId: "/api/proxy/ccp4i2/files/200/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [{ resName: "LIG" } as Partial<moorhen.LigandInfo>],
          dicts: { LIG: "data_comp_LIG\nB chemistry\n" },
        }),
      ],
      glRef: fakeGlRef,
    });
    const dictRefs = scene.files!.filter((f) => f.kind === "dictionary");
    expect(dictRefs).toHaveLength(2);
    // Different dict bodies — that's the whole point of per-molecule scoping.
    expect(dictRefs[0].cifText).toContain("A chemistry");
    expect(dictRefs[1].cifText).toContain("B chemistry");
    // Each element references its own dict.
    expect(scene.elements![0].dictionaries).toEqual([dictRefs[0].name]);
    expect(scene.elements![1].dictionaries).toEqual([dictRefs[1].name]);
  });

  it("skips non-standard ligands when no dict is stored on the molecule", () => {
    const scene = liftScene({
      molecules: [
        fakeMol({
          name: "complex",
          molNo: 0,
          uniqueId: "/api/proxy/ccp4i2/files/100/download/",
          representations: [
            { style: "CRs", visible: true, colourRules: [] } as Partial<moorhen.MoleculeRepresentation>,
          ],
          ligands: [{ resName: "MYS" } as Partial<moorhen.LigandInfo>],
          // No dicts entry for MYS — getDict returns "" so we skip emission.
          dicts: {},
        }),
      ],
      glRef: fakeGlRef,
    });
    const dictRefs = scene.files!.filter((f) => f.kind === "dictionary");
    expect(dictRefs).toHaveLength(0);
  });
});

function fakeMap(opts: {
  name: string;
  molNo: number;
  uniqueId: string;
  isDifference?: boolean;
  selectedColumns?: Partial<moorhen.selectedMtzColumns>;
}): moorhen.Map {
  return {
    name: opts.name,
    molNo: opts.molNo,
    uniqueId: opts.uniqueId,
    isDifference: !!opts.isDifference,
    selectedColumns: (opts.selectedColumns ?? { F: "FWT", PHI: "PHWT" }) as moorhen.selectedMtzColumns,
  } as unknown as moorhen.Map;
}

describe("liftScene — maps", () => {
  it("emits an MTZ file ref + SceneMap with columns + isDifference", () => {
    const scene = liftScene({
      molecules: [],
      glRef: fakeGlRef,
      maps: [
        fakeMap({
          name: "best",
          molNo: 5,
          uniqueId: "https://example/best.mtz",
          selectedColumns: { F: "FWT", PHI: "PHWT" },
        }),
        fakeMap({
          name: "diff",
          molNo: 6,
          uniqueId: "https://example/diff.mtz",
          isDifference: true,
          selectedColumns: { F: "DELFWT", PHI: "PHDELWT" },
        }),
      ],
    });
    const mtzFiles = scene.files!.filter((f) => f.kind === "mtz");
    expect(mtzFiles).toHaveLength(2);
    expect(scene.maps).toHaveLength(2);
    expect(scene.maps![0]).toMatchObject({
      file: mtzFiles[0].name,
      columns: { F: "FWT", PHI: "PHWT" },
    });
    expect(scene.maps![0].isDifference).toBeUndefined();
    expect(scene.maps![1].isDifference).toBe(true);
  });

  it("uses fileId+projectId for ccp4i2-loaded MTZs", () => {
    const scene = liftScene({
      molecules: [],
      glRef: fakeGlRef,
      projectId: "proj-uuid",
      projectName: "Demo",
      maps: [
        fakeMap({
          name: "best",
          molNo: 5,
          uniqueId: "/api/proxy/ccp4i2/files/777/download/",
        }),
      ],
    });
    const mtzFile = scene.files!.find((f) => f.kind === "mtz");
    expect(mtzFile).toMatchObject({
      fileId: 777,
      projectId: "proj-uuid",
      projectName: "Demo",
    });
  });

  it("promotes activeMap by molNo to scene.activeMap", () => {
    const scene = liftScene({
      molecules: [],
      glRef: fakeGlRef,
      maps: [
        fakeMap({ name: "best", molNo: 5, uniqueId: "https://example/best.mtz" }),
        fakeMap({ name: "diff", molNo: 6, uniqueId: "https://example/diff.mtz", isDifference: true }),
      ],
      activeMapMolNo: 5,
    });
    expect(scene.activeMap).toBe(scene.maps![0].name);
  });

  it("captures render state from mapState keyed by molNo", () => {
    const scene = liftScene({
      molecules: [],
      glRef: fakeGlRef,
      maps: [
        fakeMap({ name: "best", molNo: 5, uniqueId: "https://example/best.mtz" }),
      ],
      mapState: {
        5: {
          contourLevel: 1.234,
          radius: 17.5,
          alpha: 0.9,
          style: "solid",
          colour: "#336699",
          visible: false,
        },
      },
    });
    expect(scene.maps![0]).toMatchObject({
      contourLevel: 1.234,
      radius: 17.5,
      alpha: 0.9,
      style: "solid",
      colour: "#336699",
      visible: false,
    });
  });

  it("emits positive/negative colours only for difference maps", () => {
    const scene = liftScene({
      molecules: [],
      glRef: fakeGlRef,
      maps: [
        fakeMap({ name: "diff", molNo: 5, uniqueId: "https://example/diff.mtz", isDifference: true }),
      ],
      mapState: {
        5: { colour: "#000000", positiveColour: "#00ff00", negativeColour: "#ff0000" },
      },
    });
    const m = scene.maps![0];
    expect(m.positiveColour).toBe("#00ff00");
    expect(m.negativeColour).toBe("#ff0000");
    expect(m.colour).toBeUndefined();
  });
});
