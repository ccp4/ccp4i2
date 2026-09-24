import { describe, expect, it } from "vitest";
import {
  describeEmdbFile,
  halfMapSibling,
  otherHalfMap,
  pickDefaultEmdbFile,
  requiredSubTypes,
  type EmdbFile,
} from "../lib/emdb-fetch";

const f = (kind: EmdbFile["kind"], file: string, sub_type: number, index = 1): EmdbFile => ({
  kind, file, sub_type, index, label: kind, pixel_spacing: 0.5332, dimensions: [256, 256, 256],
  size_kbytes: 67109, contour_level: null,
});
const FILES = [f("map", "emd_1.map.gz", 1), f("half_map", "emd_1_half_map_1.map.gz", 5, 1),
  f("half_map", "emd_1_half_map_2.map.gz", 5, 2), f("mask", "emd_1_msk_1.map", 4)];

describe("requiredSubTypes", () => {
  it("reads every shape the qualifier arrives in", () => {
    expect(requiredSubTypes(5)).toEqual([5]);
    expect(requiredSubTypes("5")).toEqual([5]);
    expect(requiredSubTypes("4,5")).toEqual([4, 5]);
    expect(requiredSubTypes([1, 2])).toEqual([1, 2]);
    expect(requiredSubTypes(undefined)).toEqual([]);
    expect(requiredSubTypes(0)).toEqual([]);
  });
});

describe("pickDefaultEmdbFile", () => {
  it("preselects by the parameter's required subtype", () => {
    expect(pickDefaultEmdbFile(FILES, 5)?.file).toBe("emd_1_half_map_1.map.gz");
    expect(pickDefaultEmdbFile(FILES, 4)?.file).toBe("emd_1_msk_1.map");
  });
  it("falls back to the main map, then whatever is listed", () => {
    expect(pickDefaultEmdbFile(FILES, undefined)?.kind).toBe("map");
    expect(pickDefaultEmdbFile(FILES, 4)?.kind).toBe("mask");
    expect(pickDefaultEmdbFile(FILES.slice(1, 2), 4)?.kind).toBe("half_map");
    expect(pickDefaultEmdbFile([], 5)).toBeNull();
  });
});

describe("otherHalfMap", () => {
  it("pairs a half map with the other one and nothing else", () => {
    expect(otherHalfMap(FILES, FILES[1])?.file).toBe("emd_1_half_map_2.map.gz");
    expect(otherHalfMap(FILES, FILES[2])?.file).toBe("emd_1_half_map_1.map.gz");
    expect(otherHalfMap(FILES, FILES[0])).toBeNull();
    expect(otherHalfMap(FILES.slice(0, 2), FILES[1])).toBeNull();
  });
});

describe("halfMapSibling", () => {
  const lookup = {
    "servalcat.inputData.MAPIN1": { _class: "CMapDataFile", _qualifiers: { requiredSubType: 5 } },
    "servalcat.inputData.MAPIN2": { _class: "CMapDataFile", _qualifiers: { requiredSubType: "5" } },
    "servalcat.inputData.MAPMASK": { _class: "CMapDataFile", _qualifiers: { requiredSubType: 4 } },
    "servalcat.inputData.XYZIN": { _class: "CPdbDataFile", _qualifiers: {} },
    "servalcat.controlParameters.X": { _class: "CMapDataFile", _qualifiers: { requiredSubType: 5 } },
  };
  it("finds the other half-map input beside the one being fetched", () => {
    expect(halfMapSibling(lookup, "servalcat.inputData.MAPIN1")).toEqual({ objectPath: "servalcat.inputData.MAPIN2", name: "MAPIN2" });
    expect(halfMapSibling(lookup, "servalcat.inputData.MAPIN2")?.name).toBe("MAPIN1");
  });
  it("ignores masks, other classes, other containers, and a missing lookup", () => {
    expect(halfMapSibling(lookup, "servalcat.inputData.MAPMASK")?.name).toBe("MAPIN1");
    expect(halfMapSibling({ "t.inputData.MAPMASK": lookup["servalcat.inputData.MAPMASK"] }, "t.inputData.MAPIN")).toBeNull();
    expect(halfMapSibling(undefined, "t.inputData.MAPIN")).toBeNull();
  });
});

describe("describeEmdbFile", () => {
  it("says size, spacing and box", () => {
    expect(describeEmdbFile(FILES[1])).toBe("66 MB, 0.53 Å/px, 256³");
  });
});
