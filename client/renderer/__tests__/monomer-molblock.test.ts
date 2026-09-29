import { describe, it, expect } from "vitest";
import {
  buildMolblock,
  readCoordinates,
  type MonomerAtomDetail,
  type MonomerBond,
} from "../lib/monomer-molblock";

// A fragment of LYS, enough to exercise order, charge and bond orders.
const ATOMS: MonomerAtomDetail[] = [
  { name: "N", element: "N", charge: 1 },
  { name: "CA", element: "C", charge: 0 },
  { name: "C", element: "C", charge: 0 },
  { name: "O", element: "O", charge: 0 },
  { name: "OXT", element: "O", charge: -1 },
];
const BONDS: MonomerBond[] = [
  { atom1: "N", atom2: "CA", type: "single" },
  { atom1: "CA", atom2: "C", type: "single" },
  { atom1: "C", atom2: "O", type: "double" },
  { atom1: "C", atom2: "OXT", type: "single" },
];

const lines = (mb: string) => mb.split("\n");

describe("buildMolblock", () => {
  it("writes a counts line agreeing with what follows", () => {
    const mb = buildMolblock(ATOMS, BONDS);
    const counts = lines(mb)[3];
    expect(counts.slice(0, 3)).toBe(padded(ATOMS.length));
    expect(counts.slice(3, 6)).toBe(padded(BONDS.length));
    expect(counts).toContain("V2000");
  });

  it("keeps atoms in the order given, which is the whole point", () => {
    // Atom i of the molblock must be ATOMS[i]: that is what lets the picker
    // label a depiction with dictionary names instead of indices.
    const atomLines = lines(buildMolblock(ATOMS, BONDS)).slice(4, 4 + ATOMS.length);
    expect(atomLines.map((l) => l.slice(31, 34).trim())).toEqual(["N", "C", "C", "O", "O"]);
  });

  it("indexes bonds 1-based against that order", () => {
    const bondLines = lines(buildMolblock(ATOMS, BONDS)).slice(
      4 + ATOMS.length,
      4 + ATOMS.length + BONDS.length
    );
    // N-CA is atoms 1-2 single; C=O is atoms 3-4 double.
    expect(bondLines[0].slice(0, 9)).toBe(padded(1) + padded(2) + padded(1));
    expect(bondLines[2].slice(0, 9)).toBe(padded(3) + padded(4) + padded(2));
  });

  it("carries formal charges in an M CHG block, in FOUR-wide fields", () => {
    // Four, not three like the counts and bond blocks. Getting this wrong
    // does not fail loudly: RDKit rejects the whole molecule, and only for
    // some inputs -- LYS and ATP came back null while PLP and TRP parsed.
    const chg = lines(buildMolblock(ATOMS, BONDS)).find((l) => l.startsWith("M  CHG"));
    // Two charged atoms: N (+1, index 1) and OXT (-1, index 5).
    expect(chg).toBe("M  CHG  2   1   1   5  -1");
  });

  it("maps delocalised bond types to the aromatic code, not to single", () => {
    // CCP4 dictionaries use both spellings; dropping them to single draws
    // rings wrongly.
    for (const type of ["aromatic", "deloc"]) {
      const mb = buildMolblock(ATOMS.slice(0, 2), [{ atom1: "N", atom2: "CA", type }]);
      const bond = lines(mb)[4 + 2];
      expect(bond.slice(6, 9)).toBe(padded(4));
    }
  });

  it("skips a bond naming an atom that is not in the list", () => {
    // Hydrogens are dropped server-side; their bonds must not survive as
    // references to atom 0.
    const mb = buildMolblock(ATOMS, [...BONDS, { atom1: "N", atom2: "H1", type: "single" }]);
    expect(lines(mb)[3].slice(3, 6)).toBe(padded(BONDS.length));
  });

  it("ends with M END", () => {
    expect(lines(buildMolblock(ATOMS, BONDS)).filter(Boolean).pop()).toBe("M  END");
  });
});

describe("readCoordinates", () => {
  it("reads back one coordinate per atom", () => {
    const coords = readCoordinates(buildMolblock(ATOMS, BONDS));
    expect(coords).toHaveLength(ATOMS.length);
    expect(coords.every((c) => Number.isFinite(c.x) && Number.isFinite(c.y))).toBe(true);
  });

  it("returns nothing for a molblock it cannot parse", () => {
    expect(readCoordinates("")).toEqual([]);
    expect(readCoordinates("not\na\nmolblock\n")).toEqual([]);
  });
});

function padded(value: number): string {
  const s = String(value);
  return " ".repeat(Math.max(0, 3 - s.length)) + s;
}
