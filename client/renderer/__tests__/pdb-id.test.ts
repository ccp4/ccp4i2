import { describe, it, expect } from "vitest";
import {
  toExtendedPdbId,
  toShortPdbId,
  toDisplayPdbId,
  pdbEntryPayload,
} from "../components/task/task-elements/pdb-id";

describe("toExtendedPdbId", () => {
  it("expands a legacy 4-char code", () => {
    expect(toExtendedPdbId("1jst")).toBe("pdb_00001jst");
    expect(toExtendedPdbId("1JST")).toBe("pdb_00001jst");
    expect(toExtendedPdbId("  1jst  ")).toBe("pdb_00001jst");
  });
  it("passes an already-extended id through (lowercased)", () => {
    expect(toExtendedPdbId("pdb_00001jst")).toBe("pdb_00001jst");
    expect(toExtendedPdbId("PDB_00001JST")).toBe("pdb_00001jst");
    expect(toExtendedPdbId("pdb_00abcdef")).toBe("pdb_00abcdef");
  });
  it("leaves unrecognised input trimmed/lowercased", () => {
    expect(toExtendedPdbId("notapdb")).toBe("notapdb");
  });
});

describe("toShortPdbId", () => {
  it("returns the 4-char form when derivable", () => {
    expect(toShortPdbId("1jst")).toBe("1jst");
    expect(toShortPdbId("PDB_00001JST")).toBe("1jst");
  });
  it("returns null for a genuinely-extended id beyond the 4-char namespace", () => {
    expect(toShortPdbId("pdb_00abcdef")).toBeNull();
    expect(toShortPdbId("notapdb")).toBeNull();
  });
});

describe("toDisplayPdbId", () => {
  it("prefers the compact form, falling back to extended", () => {
    expect(toDisplayPdbId("pdb_00001jst")).toBe("1jst");
    expect(toDisplayPdbId("pdb_00abcdef")).toBe("pdb_00abcdef");
  });
});

describe("pdbEntryPayload", () => {
  it("resolves a response keyed by the legacy id when queried with the extended id", () => {
    const data = { "1jst": { PDB: { downloads: [] } } };
    expect(pdbEntryPayload(data, "pdb_00001jst")).toBe(data["1jst"]);
  });
  it("resolves a response keyed by the extended id", () => {
    const data = { pdb_00abcdef: { x: 1 } };
    expect(pdbEntryPayload(data, "pdb_00abcdef")).toBe(data["pdb_00abcdef"]);
  });
  it("falls back to the sole value when the key does not match any variant", () => {
    const data = { weirdkey: { y: 2 } };
    expect(pdbEntryPayload(data, "1jst")).toBe(data["weirdkey"]);
  });
  it("returns undefined for empty or ambiguous responses", () => {
    expect(pdbEntryPayload({}, "1jst")).toBeUndefined();
    expect(pdbEntryPayload(null, "1jst")).toBeUndefined();
    expect(pdbEntryPayload({ a: 1, b: 2 }, "1jst")).toBeUndefined();
  });
});
