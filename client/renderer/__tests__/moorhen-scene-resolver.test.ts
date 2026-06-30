/**
 * Tests for the residue-range clamping logic in the scene resolver.
 *
 * The rest of the resolver touches Moorhen runtime objects (molecules,
 * Redux dispatch, ColourRule construction) and is exercised end-to-end
 * via manual verification in the browser. The clamp logic is pure and
 * the most failure-prone bit, so it's worth a focused unit test.
 */
import { describe, it, expect } from "vitest";
import {
  buildPendingRules,
  clampRangeToPresent,
  directionToLightPosition,
  expandLsqMatches,
  geometryToM2tParams,
  isFetchable,
  resolveChainSelector,
  resolveClipFogPlanes,
  splitMultiCid,
} from "../lib/moorhen-scene-resolver";

const present = (...nums: number[]) => new Set(nums);

describe("resolveClipFogPlanes (broad → fine precedence)", () => {
  const fco = 250;
  it("slab alone sets all four planes (absolute Å about fco), freezes", () => {
    const p = resolveClipFogPlanes({}, 20, 1, fco);
    expect(p).toEqual({ clipStart: 20, clipEnd: 20, fogStart: 230, fogEnd: 270, reset: false });
  });
  it("explicit fogEnd OVERRIDES the slab-derived back fog plane (the fix)", () => {
    const p = resolveClipFogPlanes({ fogEnd: 999 }, 20, 1, fco);
    expect(p.clipStart).toBe(20); // slab kept
    expect(p.fogStart).toBe(230); // slab kept
    expect(p.fogEnd).toBe(999); // explicit wins
    expect(p.reset).toBe(false);
  });
  it("clip:{front,back} sets zoom-scaled field-depth planes", () => {
    const p = resolveClipFogPlanes({ clip: { front: 8, back: 21 } }, undefined, 2, fco);
    expect(p).toEqual({ clipStart: 16, clipEnd: 42, fogStart: 234, fogEnd: 292, reset: false });
  });
  it("clip:\"auto\" alone = recompute (reset true), no planes", () => {
    expect(resolveClipFogPlanes({ clip: "auto" }, undefined, 1, fco)).toEqual({ reset: true });
  });
  it("explicit overrides clip:\"auto\" (finer wins) and freezes", () => {
    const p = resolveClipFogPlanes({ clip: "auto", fogEnd: 300 }, undefined, 1, fco);
    expect(p.fogEnd).toBe(300);
    expect(p.reset).toBe(false);
  });
  it("nothing specified → leave untouched (reset undefined)", () => {
    expect(resolveClipFogPlanes({ origin: [0, 0, 0] }, undefined, 1, fco)).toEqual({ reset: undefined });
  });
  it("clip:\"lock\" alone freezes current with no planes", () => {
    expect(resolveClipFogPlanes({ clip: "lock" }, undefined, 1, fco)).toEqual({ reset: false });
  });
});

describe("buildPendingRules cascade (element.colour ↔ representation.colour)", () => {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const base: any = {
    molecule: {},
    domains: [],
    fileName: "f",
    log: [],
    policy: "clamp-and-log",
    dispatch: () => {},
  };

  it("falls back to element.colour when the representation has none", () => {
    const rules = buildPendingRules(
      { ...base, rep: { style: "CRs" }, elementColour: "#abcdef" },
      "//A",
    );
    expect(rules[0].color).toBe("#abcdef");
  });

  it("representation.colour overrides element.colour", () => {
    const rules = buildPendingRules(
      { ...base, rep: { style: "CRs", colour: "#123456" }, elementColour: "#abcdef" },
      "//A",
    );
    expect(rules[0].color).toBe("#123456");
  });

  it("no colour at either level → no rules", () => {
    expect(buildPendingRules({ ...base, rep: { style: "CRs" } }, "//A")).toEqual([]);
  });
});

describe("geometryToM2tParams (honoured geometry → Moorhen m2tParameters)", () => {
  it("maps each field to its m2tParameters key, merged onto the base", () => {
    const base = { ribbonStyleHelixWidth: 1.0, cylindersStyleCylinderRadius: 0.1 };
    const out = geometryToM2tParams(
      { bondRadius: 0.18, probeRadius: 1.4, ribbonHelixWidth: 1.4 },
      base,
    );
    expect(out.cylindersStyleCylinderRadius).toBe(0.18); // overrides base
    expect(out.surfaceStyleProbeRadius).toBe(1.4);
    expect(out.ribbonStyleHelixWidth).toBe(1.4);
  });

  it("leaves unspecified params untouched", () => {
    const out = geometryToM2tParams({ bondRadius: 0.2 }, { foo: 1 });
    expect(out.foo).toBe(1);
    expect(out.cylindersStyleCylinderRadius).toBe(0.2);
    expect("surfaceStyleProbeRadius" in out).toBe(false);
  });
});

describe("directionToLightPosition (scene direction → Moorhen lightPosition)", () => {
  it("normalises to Moorhen's default magnitude with w=1", () => {
    expect(directionToLightPosition([0, 0, 1])).toEqual([0, 0, 60, 1]);
  });
  it("preserves direction but fixes the magnitude (not passed raw)", () => {
    const p = directionToLightPosition([0.2, 0.2, 1.0]);
    expect(Math.hypot(p[0], p[1], p[2])).toBeCloseTo(60, 5);
    expect(p[3]).toBe(1);
    // direction preserved: ratios match the input
    expect(p[0] / p[2]).toBeCloseTo(0.2, 5);
  });
});

describe("clampRangeToPresent", () => {
  it("returns the request unchanged when all residues are present", () => {
    const p = present(1, 2, 3, 4, 5);
    expect(clampRangeToPresent(1, 5, p)).toEqual([[1, 5]]);
  });

  it("clamps the start inward when leading residues are missing", () => {
    const p = present(5, 6, 7, 8, 9, 10);
    expect(clampRangeToPresent(1, 10, p)).toEqual([[5, 10]]);
  });

  it("clamps the end inward when trailing residues are missing", () => {
    const p = present(1, 2, 3, 4, 5);
    expect(clampRangeToPresent(1, 10, p)).toEqual([[1, 5]]);
  });

  it("splits across internal gaps", () => {
    const p = present(1, 2, 3, 7, 8, 9);
    expect(clampRangeToPresent(1, 10, p)).toEqual([
      [1, 3],
      [7, 9],
    ]);
  });

  it("returns empty when no requested residue is present", () => {
    const p = present(50, 51, 52);
    expect(clampRangeToPresent(1, 10, p)).toEqual([]);
  });

  it("handles a single-residue request", () => {
    expect(clampRangeToPresent(42, 42, present(42))).toEqual([[42, 42]]);
    expect(clampRangeToPresent(42, 42, present(41, 43))).toEqual([]);
  });

  it("handles multiple gaps", () => {
    const p = present(1, 2, 4, 5, 7, 9, 10);
    expect(clampRangeToPresent(1, 10, p)).toEqual([
      [1, 2],
      [4, 5],
      [7, 7],
      [9, 10],
    ]);
  });

  it("handles a request narrower than the present set", () => {
    const p = present(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    expect(clampRangeToPresent(3, 6, p)).toEqual([[3, 6]]);
  });
});

describe("resolveChainSelector", () => {
  const presentByChain = new Map<string, Set<number>>([
    ["A", present(1, 2)],
    ["B", present(1, 2)],
    ["C", present(1, 2)],
  ]);

  it("returns a single-element array for a non-wildcard string", () => {
    expect(resolveChainSelector("A", presentByChain)).toEqual(["A"]);
  });

  it('expands "*" to every present chain, sorted', () => {
    expect(resolveChainSelector("*", presentByChain)).toEqual(["A", "B", "C"]);
  });

  it('returns empty for "*" against an unloaded molecule', () => {
    expect(resolveChainSelector("*", new Map())).toEqual([]);
  });

  it("returns explicit lists verbatim", () => {
    expect(resolveChainSelector(["B", "A"], presentByChain)).toEqual(["B", "A"]);
  });

  it("preserves single-chain references to chains not present", () => {
    // The resolver is responsible for surfacing absent-chain errors later;
    // selector resolution only decides which chains we're targeting.
    expect(resolveChainSelector("Z", presentByChain)).toEqual(["Z"]);
  });
});

describe("isFetchable", () => {
  it("returns true for a pdb: ref", () => {
    expect(isFetchable({ name: "x", pdb: "1m17" })).toBe(true);
  });

  it("returns true for a url: ref", () => {
    expect(isFetchable({ name: "x", url: "https://example.com/foo.cif" })).toBe(true);
  });

  it("returns true for a fileId + projectId ref", () => {
    expect(
      isFetchable({ name: "x", fileId: 42, projectId: "uuid" }),
    ).toBe(true);
  });

  it("returns false for a fileId without projectId", () => {
    expect(isFetchable({ name: "x", fileId: 42 })).toBe(false);
  });

  it("returns false for a relativeUrl-only ref (resolver can't fetch it standalone)", () => {
    expect(isFetchable({ name: "x", relativeUrl: "/api/proxy/pdbe/x.cif" })).toBe(false);
  });

  it("returns false for an empty ref", () => {
    expect(isFetchable({ name: "x" })).toBe(false);
  });
});

describe("expandLsqMatches", () => {
  it("expands chain+range shorthand into a single same-on-both-sides match", () => {
    expect(expandLsqMatches({ chain: "A", range: "104-260" })).toEqual([
      { refChain: "A", movChain: "A", refRange: "104-260", movRange: "104-260" },
    ]);
  });

  it("returns explicit matches verbatim when no shorthand is set", () => {
    const matches = [
      { refChain: "A", refRange: "1-100", movChain: "A", movRange: "1-100" },
      { refChain: "B", refRange: "50-150", movChain: "B", movRange: "60-160" },
    ];
    expect(expandLsqMatches({ matches })).toEqual(matches);
  });

  it("prefers the shorthand when both are present (shouldn't happen — validator rejects this)", () => {
    // Defensive: if the validator ever lets both through, shorthand wins
    // because that's the more specific instruction.
    const result = expandLsqMatches({
      chain: "A",
      range: "104-260",
      matches: [
        { refChain: "B", refRange: "1-50", movChain: "B", movRange: "1-50" },
      ],
    });
    expect(result).toEqual([
      { refChain: "A", movChain: "A", refRange: "104-260", movRange: "104-260" },
    ]);
  });

  it("returns an empty list when neither shorthand nor matches is set", () => {
    expect(expandLsqMatches({})).toEqual([]);
  });
});

describe("splitMultiCid", () => {
  it("returns a single-element array for a plain CID", () => {
    expect(splitMultiCid("//A/100-200")).toEqual(["//A/100-200"]);
  });

  it("splits a ||-joined multi-CID into chunks", () => {
    expect(splitMultiCid("//A/115||//A/116||//A/121-122")).toEqual([
      "//A/115",
      "//A/116",
      "//A/121-122",
    ]);
  });

  it("trims whitespace and drops empty chunks", () => {
    expect(splitMultiCid("//A/1 ||  //A/2 ||  || //A/3")).toEqual([
      "//A/1",
      "//A/2",
      "//A/3",
    ]);
  });

  it("returns the wildcard verbatim", () => {
    expect(splitMultiCid("/*/*/*/*")).toEqual(["/*/*/*/*"]);
  });

  it("returns [] for an all-empty multi-CID", () => {
    expect(splitMultiCid("||||")).toEqual([]);
  });
});
