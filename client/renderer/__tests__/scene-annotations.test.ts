/**
 * Annotations: scene `annotations[]` → Moorhen vectors and 3D objects
 * (resolver) and back (lifter). The geometry rules are pure, so they are
 * tested without a store.
 */
import { describe, it, expect } from "vitest";

import {
  ANNOTATION_OBJECT_TAGS,
  ANNOTATION_TAG,
  ANNOTATION_VECTOR_PREFIX,
  annotationFromObject,
  annotationFromVector,
  annotationIndexOf,
  annotationsToDrawables,
} from "../lib/moorhen-scene-resolver";
import { liftScene } from "../lib/moorhen-scene-lifter";
import { parseScene } from "../lib/scene";
import type { moorhen } from "moorhen/types/moorhen";
import type { SceneAnnotation, ScenePoint } from "../types/moorhen-scene";

type Xyz = [number, number, number];

/** Anchors resolve to a fixed point; tuples pass through. */
const resolveAll = (pt: ScenePoint): Xyz | null => (Array.isArray(pt) ? pt : [1, 2, 3]);

/** Lower, and give back just the vectors: the label and line tests' shape. */
const annotationsToVectors: typeof annotationsToDrawables extends (...a: infer A) => unknown
  ? (...a: A) => ReturnType<typeof annotationsToDrawables>["vectors"]
  : never = (...a) => annotationsToDrawables(...a).vectors;

describe("annotationsToDrawables", () => {
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

  it("a disc is a cylinder object of its radius, `thickness` long, centred along a unit normal", () => {
    const { vectors, objects } = annotationsToDrawables(
      [{ kind: "disc", at: [10, 10, 10], radius: 2, thickness: 1, normal: [0, 3, 0], text: "7 hits",
         colour: "#E07A5F" }],
      resolveAll,
    );
    expect(objects).toHaveLength(1);
    const [o] = objects;
    expect(o.type).toBe("cylinder");
    if (o.type !== "cylinder") return;
    expect(o.origin).toEqual([10, 9.5, 10]);
    expect(o.end).toEqual([10, 10.5, 10]);
    expect(o.radius).toBe(2);
    expect(o.colour).toBe("#e07a5fff");          // objects take #rrggbbaa
    expect(o.locked).toBe(true);
    expect(o.tags?.source).toBe(ANNOTATION_OBJECT_TAGS.source);
    expect(o.uniqueId).toBe(`${ANNOTATION_VECTOR_PREFIX}:0`);
    // The banner is a zero-length vector on the +normal face, not stamped:
    // the object carries the annotation.
    expect(vectors).toHaveLength(1);
    const [v] = vectors;
    expect([v.xFrom, v.yFrom, v.zFrom]).toEqual([10, 10.5, 10]);
    expect([v.xTo, v.yTo, v.zTo]).toEqual([10, 10.5, 10]);
    expect(v.labelMode).toBe("end");
    expect(v.labelText).toBe("7 hits");
    expect(v.uniqueId).toBe(`${ANNOTATION_VECTOR_PREFIX}:0#label`);
    expect(annotationFromVector(v)).toBeNull();
    expect(annotationFromObject(o)).not.toBeNull();
  });

  it("a disc defaults to 0.3 Å thick along z, a zero normal falls back to z, and no text means no banner", () => {
    const { vectors, objects } = annotationsToDrawables(
      [{ kind: "disc", at: [0, 0, 0], radius: 1, normal: [0, 0, 0] }],
      resolveAll,
    );
    const [o] = objects;
    if (o.type !== "cylinder") throw new Error("expected a cylinder");
    expect(o.origin[2]).toBeCloseTo(-0.15);
    expect(o.end[2]).toBeCloseTo(0.15);
    expect(o.colour).toBe("#f2c14eff");           // the default colour
    expect(vectors).toHaveLength(0);
  });

  it("a sphere is a sphere object at its centre, with its banner there too", () => {
    const { vectors, objects } = annotationsToDrawables(
      [{ kind: "sphere", at: [10, 20, 30], radius: 2, text: "7 hits", colour: "#0000ff" }],
      resolveAll,
    );
    expect(objects).toHaveLength(1);
    const [o] = objects;
    expect(o.type).toBe("sphere");
    if (o.type !== "sphere") return;
    expect(o.origin).toEqual([10, 20, 30]);
    expect(o.radius).toBe(2);
    expect(o.colour).toBe("#0000ffff");
    expect(o.uniqueId).toBe(`${ANNOTATION_VECTOR_PREFIX}:0`);
    expect(vectors).toHaveLength(1);
    expect([vectors[0].xFrom, vectors[0].yFrom, vectors[0].zFrom]).toEqual([10, 20, 30]);
    expect(vectors[0].labelText).toBe("7 hits");
    expect(vectors[0].vectorColour.b).toBe(255);
    expect(vectors[0].uniqueId).toBe(`${ANNOTATION_VECTOR_PREFIX}:0#label`);
    expect(annotationFromVector(vectors[0])).toBeNull();
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
    expect(vs.map((v) => annotationIndexOf(v.uniqueId))).toEqual([0, 2]);
    expect(annotationIndexOf(`${ANNOTATION_VECTOR_PREFIX}:7#label`)).toBe(7);
    expect(annotationIndexOf("user-made")).toBeNull();
  });

  it("stamps each drawable with the annotation as authored and as lowered", () => {
    const authored: SceneAnnotation = {
      kind: "line", from: { file: "p", selection: "//A/1" }, to: [4, 4, 4], text: "x",
    };
    const [v] = annotationsToVectors([authored], resolveAll);
    expect(v.customTags?.[0]).toBe(ANNOTATION_TAG);
    const read = annotationFromVector(v);
    expect(read?.authored).toEqual(authored);
    expect(read?.lowered).toEqual({ ...authored, from: [1, 2, 3] });

    const sphere: SceneAnnotation = { kind: "sphere", at: { file: "p", selection: "//A/1" }, radius: 1 };
    const [o] = annotationsToDrawables([sphere], resolveAll).objects;
    const readObj = annotationFromObject(o);
    expect(readObj?.authored).toEqual(sphere);
    expect(readObj?.lowered).toEqual({ ...sphere, at: [1, 2, 3] });
  });

  it("annotationFromVector and annotationFromObject ignore drawables that are not the scene's", () => {
    expect(annotationFromVector({ uniqueId: "user-made", customTags: [] })).toBeNull();
    expect(annotationFromVector({ uniqueId: `${ANNOTATION_VECTOR_PREFIX}:0`, customTags: ["other"] })).toBeNull();
    expect(
      annotationFromVector({ uniqueId: `${ANNOTATION_VECTOR_PREFIX}:0`, customTags: [ANNOTATION_TAG, "{not json"] }),
    ).toBeNull();
    expect(annotationFromObject({ uniqueId: "user-made" })).toBeNull();
    expect(annotationFromObject({ uniqueId: "x", tags: { source: "xpid" } })).toBeNull();
    expect(annotationFromObject({ uniqueId: "x", tags: { ...ANNOTATION_OBJECT_TAGS } })).toBeNull();
    expect(annotationFromObject({ uniqueId: "x", tags: { ...ANNOTATION_OBJECT_TAGS, annotation: "{no" } })).toBeNull();
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

  it("lifts a sphere from its object once, not again from its banner", () => {
    const sphere: SceneAnnotation = { kind: "sphere", at: [1, 2, 3], radius: 1.2, text: "site" };
    const { vectors, objects } = annotationsToDrawables([sphere], resolveAll);
    expect(vectors).toHaveLength(1);
    const scene = liftScene({ molecules: [mol], glRef, vectors, objects });
    expect(scene.annotations).toEqual([sphere]);
  });

  it("restores the authored order across the vectors and objects stores", () => {
    const authored: SceneAnnotation[] = [
      { kind: "sphere", at: [1, 2, 3], radius: 1 },
      { kind: "label", at: [4, 5, 6], text: "a" },
      { kind: "disc", at: [7, 8, 9], radius: 2 },
      { kind: "line", from: [0, 0, 0], to: [1, 1, 1] },
    ];
    const { vectors, objects } = annotationsToDrawables(authored, resolveAll);
    const scene = liftScene({ molecules: [mol], glRef, vectors, objects });
    expect(scene.annotations).toEqual(authored);
  });

  it("leaves the user's own objects out", () => {
    const scene = liftScene({
      molecules: [mol], glRef,
      objects: [{ uniqueId: "mine", tags: { source: "xpid" } }, { uniqueId: "untagged" }],
    });
    expect(scene.annotations).toBeUndefined();
  });

  it("leaves the user's own vectors out, and omits the block when there is nothing", () => {
    const scene = liftScene({
      molecules: [mol], glRef,
      vectors: [{ uniqueId: "mine", customTags: [] }],
    });
    expect(scene.annotations).toBeUndefined();
  });
});
