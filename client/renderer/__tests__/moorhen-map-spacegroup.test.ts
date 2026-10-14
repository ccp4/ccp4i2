import { describe, expect, it } from "vitest";
import { ccp4RestoreSpaceGroup } from "../lib/moorhen-map-file";

/** A bare little-endian CCP4 header with a cell and a space group. */
function header(ispg: number, cell = [82.1, 96.6, 57.9, 90, 90, 90]): ArrayBuffer {
  const buf = new ArrayBuffer(1024);
  const dv = new DataView(buf);
  cell.forEach((v, i) => dv.setFloat32(40 + 4 * i, v, true));
  dv.setInt32(88, ispg, true);
  dv.setUint8(212, 0x44);
  return buf;
}
const ispgOf = (b: ArrayBuffer) => new DataView(b).getInt32(88, true);
const BAZ2B = { ccp4: 20, cell: [82.1, 96.6, 57.9, 90, 90, 90] as [number, number, number, number, number, number] };

describe("ccp4RestoreSpaceGroup", () => {
  it("gives a P1 map on a model's cell the model's space group", () => {
    expect(ispgOf(ccp4RestoreSpaceGroup(header(1), [BAZ2B]))).toBe(20);
    expect(ispgOf(ccp4RestoreSpaceGroup(header(0), [BAZ2B]))).toBe(20);
  });

  it("leaves a map on another cell alone (an EM map, another crystal)", () => {
    const em = header(1, [200, 200, 200, 90, 90, 90]);
    expect(ccp4RestoreSpaceGroup(em, [BAZ2B])).toBe(em);
  });

  it("leaves a map that already has a space group alone", () => {
    const xtal = header(19);
    expect(ccp4RestoreSpaceGroup(xtal, [BAZ2B])).toBe(xtal);
  });

  it("does nothing with no crystal to match", () => {
    const m = header(1);
    expect(ccp4RestoreSpaceGroup(m, [])).toBe(m);
  });
});
