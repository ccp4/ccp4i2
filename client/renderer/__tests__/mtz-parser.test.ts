import { describe, it, expect } from "vitest";
import { existsSync, readFileSync } from "fs";
import path from "path";
import { parseMtzHeader, getColumnNamesByType } from "../lib/mtz-parser";

/**
 * Build a minimal little-endian MTZ: magic, header pointer, machine stamp,
 * no reflections, then the header as 80-character records.
 */
function mtzWithHeader(records: string[]): ArrayBuffer {
  const headerText = records
    .concat(["END", "MTZENDOFHEADERS"])
    .map((r) => r.padEnd(80))
    .join("");
  const headerByteOffset = 80; // after the 20-byte stamp block, word-aligned
  const buf = new ArrayBuffer(headerByteOffset + headerText.length);
  const view = new DataView(buf);
  const bytes = new Uint8Array(buf);
  bytes.set([0x4d, 0x54, 0x5a, 0x20], 0); // "MTZ "
  view.setInt32(4, headerByteOffset / 4 + 1, true); // 1-indexed word offset
  bytes.set([0x44, 0x41, 0x00, 0x00], 8); // little-endian machine stamp
  for (let i = 0; i < headerText.length; i++) {
    bytes[headerByteOffset + i] = headerText.charCodeAt(i);
  }
  return buf;
}

const PREAMBLE = [
  "VERS MTZ:V1.1",
  "TITLE test",
  "NCOL        5            3        0",
  "CELL    75.1100   75.1100  133.3100   90.0000   90.0000  120.0000",
  "SORT    1   2   3   0   0",
  "SYMINF   6  6 P   154 P3221      PG321",
  "SYMM X,   Y,   Z",
  "RESO      0.00448     0.11084",
  "VALM  NAN",
];

describe("parseMtzHeader COLUMN records", () => {
  it("reads five-field records, with their dataset ids", () => {
    const header = parseMtzHeader(
      mtzWithHeader([
        ...PREAMBLE,
        "COLUMN H                              H            0.0000          24.0000    0",
        "COLUMN K                              H          -12.0000           0.0000    0",
        "COLUMN L                              H          -44.0000          44.0000    0",
        "COLUMN F                              F            0.3200         137.1900    1",
        "COLUMN SIGF                           Q            0.1100          38.9500    1",
      ])
    );
    expect(header.columns.map((c) => [c.label, c.type, c.datasetId])).toEqual([
      ["H", "H", 0], ["K", "H", 0], ["L", "H", 0], ["F", "F", 1], ["SIGF", "Q", 1],
    ]);
  });

  it("reads pre-dataset four-field records as dataset 0", () => {
    // beta_blip_P3221.mtz in demo_data is written this way. With the dataset
    // id required, it parsed to no columns at all, and the column dialog
    // never opened for it.
    const header = parseMtzHeader(
      mtzWithHeader([
        ...PREAMBLE,
        "COLUMN H                              H            0.0000          24.0000",
        "COLUMN K                              H          -12.0000           0.0000",
        "COLUMN L                              H          -44.0000          44.0000",
        "COLUMN Fobs                           F            0.3200         137.1900",
        "COLUMN Sigma                          Q            0.1100          38.9500",
      ])
    );
    expect(header.columns.map((c) => [c.label, c.type, c.datasetId])).toEqual([
      ["H", "H", 0], ["K", "H", 0], ["L", "H", 0], ["Fobs", "F", 0], ["Sigma", "Q", 0],
    ]);
    expect(getColumnNamesByType(header)).toEqual({
      H: "H", K: "H", L: "H", Fobs: "F", Sigma: "Q",
    });
  });

  it("reads the real beta_blip_P3221.mtz when the server tree is present", () => {
    const file = path.resolve(
      __dirname,
      "../../../server/ccp4i2/demo_data/beta_blip/beta_blip_P3221.mtz"
    );
    if (!existsSync(file)) return;
    const buf = readFileSync(file);
    const header = parseMtzHeader(
      buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)
    );
    expect(header.columns.map((c) => c.label)).toEqual(["H", "K", "L", "Fobs", "Sigma"]);
    expect(header.columns.every((c) => c.datasetId === 0)).toBe(true);
  });
});
