/**
 * Annotations: scene `annotations[]` → Moorhen vectors (resolver) and back
 * (lifter). The geometry rules are pure, so they are tested without a store.
 */
import { describe, it, expect } from "vitest";

import {
  ANNOTATION_TAG,
  ANNOTATION_VECTOR_PREFIX,
  annotationFromVector,
  annotationsToVectors,
} from "../lib/moorhen-scene-resolver";
import { liftScene } from "../lib/moorhen-scene-lifter";
import { parseScene } from "../lib/scene";
import type { moorhen } from "moorhen/types/moorhen";
import type { SceneAnnotation, ScenePoint } from "../types/moorhen-scene";

type Xyz = [number, number, number];

/** Anchors resolve to a fixed point; tuples pass through. */
const resolveAll = (pt: ScenePoint): Xyz | null => (Array.isArray(pt) ? pt : [1, 2, 3]);

describe("annotationsToVectors", () => {
  it("a label is a zero-length vector with its banner at the end", () => {
    const [v] = annotationsToVectors(
      [{ kind: "label", at: [5, 6, 7], text: "here", colour: "#ff0000" }],
      resolveAll,
    );
    expect([v.xFrom, v.yFrom, v.zFrom]).toEqual([5, 6, 7]);
    expect([v.xTo, v.yTo, v.zTo]).toEqual([5, 6, 7]);
    expect(v.labelMode).toBe("end");
    expect(v.labelText).toBe("here");
    expect(v.arrowMode).toBe("none");          // an arrow on a zero-length vector divides by zero
    expect(v.radius).toBeGreaterThan(0);        // a falsy radius takes Moorhen's 0.07 default
    expect(v.vectorColour).toEqual({ r: 255, g: 0, b: 0 });
    expect(v.textColour).toEqual({ r: 255, g: 0, b: 0 });
    expect(v.uniqueId.startsWith(ANNOTATION_VECTOR_PREFIX)).toBe(true);
    expect(v.uniqueId).toContain("__TAG");      // hidden from Moorhen's Vectors dialog
  });

  it("a line keeps its ends, dashes, arrows and banner placement", () => {
    const [v] = annotationsToVectors(
      [{ kind: "line", from: [0, 0, 0], to: [0, 0, 10], dashed: true, arrows: "end",
         text: "10 Å", textAt: "end", radius: 0.25, textColour: "#00ff00" }],
      resolveAll,
    );
    expect(v.drawMode).toBe("dashedcylinder");
    expect(v.arrowMode).toBe("end");
    expect(v.labelMode).toBe("end");
    expect(v.radius).toBe(0.25);
    expect(v.textColour).toEqual({ r: 0, g: 255, b: 0 });
  });

  it("a line with no text has no label; defaults are applied", () => {
    const [v] = annotationsToVectors([{ kind: "line", from: [0, 0, 0], to: [1, 0, 0] }], resolveAll);
    expect(v.labelMode).toBe("none");
    expect(v.drawMode).toBe("cylinder");
    expect(v.arrowMode).toBe("none");
    expect(v.radius).toBe(0.1);
  });

  it("a disc is a cylinder of its radius, `thickness` long, centred along a unit normal", () => {
    const [v] = annotationsToVectors(
      [{ kind: "disc", at: [10, 10, 10], radius: 2, thickness: 1, normal: [0, 3, 0], text: "7 hits" }],
      resolveAll,
    );
    expect([v.xFrom, v.yFrom, v.zFrom]).toEqual([10, 9.5, 10]);
    expect([v.xTo, v.yTo, v.zTo]).toEqual([10, 10.5, 10]);
    expect(v.radius).toBe(2);
    expect(v.labelMode).toBe("end");
    expect(v.labelText).toBe("7 hits");
  });

  it("a disc defaults to 0.3 Å thick along z, and a zero normal falls back to z", () => {
    const [v] = annotationsToVectors(
      [{ kind: "disc", at: [0, 0, 0], radius: 1, normal: [0, 0, 0] }],
      resolveAll,
    );
    expect(v.zFrom).toBeCloseTo(-0.15);
    expect(v.zTo).toBeCloseTo(0.15);
    expect(v.labelMode).toBe("none");
  });

  it("a sphere is a bristle of arms through its centre, stamped and labelled once", () => {
    const vs = annotationsToVectors(
      [{ kind: "sphere", at: [10, 20, 30], radius: 2, text: "7 hits", colour: "#0000ff" }],
      resolveAll,
    );
    expect(vs).toHaveLength(7);
    // Each arm spans 2r through the point: the axes first, then the diagonals.
    expect([vs[0].xFrom, vs[0].xTo, vs[0].yFrom, vs[0].zFrom]).toEqual([8, 12, 20, 30]);
    expect([vs[1].yFrom, vs[1].yTo, vs[1].xFrom]).toEqual([18, 22, 10]);
    expect([vs[2].zFrom, vs[2].zTo]).toEqual([28, 32]);
    const d = 2 / Math.sqrt(3);
    expect([vs[3].xFrom, vs[3].yFrom, vs[3].zFrom]).toEqual([10 - d, 20 - d, 30 - d]);
    expect([vs[6].xTo, vs[6].yTo, vs[6].zTo]).toEqual([10 + d, 20 - d, 30 - d]);
    expect(vs.every((v) => v.radius === 2 * 0.55)).toBe(true);
    expect(vs.every((v) => v.vectorColour.b === 255)).toBe(true);
    // The banner sits at the centre of the first arm; the others are mute.
    expect(vs.map((v) => v.labelMode)).toEqual(["middle", ...Array(6).fill("none")]);
    expect(vs[0].labelText).toBe("7 hits");
    // Same replace-on-apply prefix on every arm; one stamp, so it lifts once.
    expect(vs.every((v) => v.uniqueId.startsWith(`${ANNOTATION_VECTOR_PREFIX}:0`))).toBe(true);
    expect(vs.map((v) => annotationFromVector(v) !== null)).toEqual([true, ...Array(6).fill(false)]);
  });

  it("anchored points are resolved; an unresolvable one drops just its annotation", () => {
    const resolve = (pt: ScenePoint): Xyz | null => {
      if (Array.isArray(pt)) return pt;
      return pt.selection === "//A/1" ? [1, 1, 1] : null;
    };
    const vs = annotationsToVectors(
      [
        { kind: "label", at: { file: "p", selection: "//A/1" }, text: "a" },
        { kind: "label", at: { file: "p", selection: "//A/999" }, text: "b" },
        { kind: "line", from: { file: "p", selection: "//A/1" }, to: [4, 4, 4] },
      ],
      resolve,
    );
    expect(vs.map((v) => v.labelText)).toEqual(["a", ""]);
    expect([vs[0].xFrom, vs[1].xTo]).toEqual([1, 4]);
    // Each vector keeps its own index, so a dropped one leaves a gap, not a shift.
    expect(vs.map((v) => v.uniqueId)).toEqual([
      `${ANNOTATION_VECTOR_PREFIX}:0`, `${ANNOTATION_VECTOR_PREFIX}:2`,
    ]);
  });

  it("stamps each vector with the annotation as authored and as lowered", () => {
    const authored: SceneAnnotation = {
      kind: "line", from: { file: "p", selection: "//A/1" }, to: [4, 4, 4], text: "x",
    };
    const [v] = annotationsToVectors([authored], resolveAll);
    expect(v.customTags?.[0]).toBe(ANNOTATION_TAG);
    const read = annotationFromVector(v);
    expect(read?.authored).toEqual(authored);
    expect(read?.lowered).toEqual({ ...authored, from: [1, 2, 3] });
  });

  it("annotationFromVector ignores vectors that are not the scene's", () => {
    expect(annotationFromVector({ uniqueId: "user-made", customTags: [] })).toBeNull();
    expect(annotationFromVector({ uniqueId: `${ANNOTATION_VECTOR_PREFIX}:0`, customTags: ["other"] })).toBeNull();
    expect(
      annotationFromVector({ uniqueId: `${ANNOTATION_VECTOR_PREFIX}:0`, customTags: [ANNOTATION_TAG, "{not json"] }),
    ).toBeNull();
  });
});

describe("liftScene: annotations", () => {
  const mol = {
    name: "prot",
    molNo: 0,
    uniqueId: "https://example.org/structures/prot.pdb",
    representations: [],
    ligands: [],
    getDict: () => "",
  } as unknown as moorhen.Molecule;
  const glRef = { origin: [0, 0, 0], quat: [0, 0, 0, -1], zoom: 1 };

  const authoredOn = (file: string): SceneAnnotation => ({
    kind: "label", at: { file, selection: "//A/1" }, text: "a",
  });

  it("lifts a scene-made vector as authored when its anchor file is in the lift", () => {
    const [v] = annotationsToVectors([authoredOn("prot")], resolveAll);
    const scene = liftScene({ molecules: [mol], glRef, vectors: [v] });
    expect(scene.annotations).toEqual([authoredOn("prot")]);
    expect(() => parseScene(JSON.stringify(scene))).not.toThrow();
  });

  it("falls back to the lowered form when the anchor names a file the lift does not have", () => {
    const [v] = annotationsToVectors([authoredOn("other-name")], resolveAll);
    const scene = liftScene({ molecules: [mol], glRef, vectors: [v] });
    expect(scene.annotations).toEqual([{ kind: "label", at: [1, 2, 3], text: "a" }]);
    expect(() => parseScene(JSON.stringify(scene))).not.toThrow();
  });

  it("lifts a sphere once, not once per arm", () => {
    const sphere: SceneAnnotation = { kind: "sphere", at: [1, 2, 3], radius: 1.2, text: "site" };
    const arms = annotationsToVectors([sphere], resolveAll);
    const scene = liftScene({ molecules: [mol], glRef, vectors: arms });
    expect(scene.annotations).toEqual([sphere]);
  });

  it("leaves the user's own vectors out, and omits the block when there is nothing", () => {
    const scene = liftScene({
      molecules: [mol], glRef,
      vectors: [{ uniqueId: "mine", customTags: [] }],
    });
    expect(scene.annotations).toBeUndefined();
  });
});
