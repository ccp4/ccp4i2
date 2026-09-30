/**
 * The arithmetic between RDKit's coordinates and the SVG: scaling, the
 * y-axis flip, and the parallel lines of a multiple bond.
 */
import { describe, it, expect } from "vitest";
import {
  layoutMonomer,
  bondLines,
  multiplicity,
  bondKey,
  type PlacedBond,
} from "../components/monomer/monomer-layout";
import type { MonomerAtomDetail, MonomerBond } from "../lib/monomer-molblock";

const BOX = { width: 200, height: 100, padding: 10 };

const atoms = (...names: string[]): MonomerAtomDetail[] =>
  names.map((name) => ({ name, element: name[0], charge: 0 }));

describe("layoutMonomer", () => {
  it("fits the molecule inside the box", () => {
    const coords = [
      { x: 0, y: 0 },
      { x: 10, y: 5 },
      { x: 5, y: -5 },
    ];
    const out = layoutMonomer(atoms("N", "CA", "CB"), [], coords, BOX);
    for (const atom of out.atoms) {
      expect(atom.x).toBeGreaterThanOrEqual(BOX.padding - 0.001);
      expect(atom.x).toBeLessThanOrEqual(BOX.width - BOX.padding + 0.001);
      expect(atom.y).toBeGreaterThanOrEqual(BOX.padding - 0.001);
      expect(atom.y).toBeLessThanOrEqual(BOX.height - BOX.padding + 0.001);
    }
  });

  it("scales both axes equally, so bond lengths stay comparable", () => {
    // A square: if width and height were fitted independently it would come
    // out as a rectangle, and every bond length would be a lie.
    const coords = [
      { x: 0, y: 0 },
      { x: 1, y: 0 },
      { x: 1, y: 1 },
      { x: 0, y: 1 },
    ];
    const out = layoutMonomer(atoms("A", "B", "C", "D"), [], coords, BOX);
    const side1 = Math.hypot(out.atoms[1].x - out.atoms[0].x, out.atoms[1].y - out.atoms[0].y);
    const side2 = Math.hypot(out.atoms[2].x - out.atoms[1].x, out.atoms[2].y - out.atoms[1].y);
    expect(side1).toBeCloseTo(side2, 6);
  });

  it("flips y, because SVG grows downwards and chemistry drawings do not", () => {
    const coords = [
      { x: 0, y: 0 },
      { x: 0, y: 10 },
    ];
    const out = layoutMonomer(atoms("LOW", "HIGH"), [], coords, BOX);
    // The chemically higher atom must be drawn nearer the top, i.e. smaller y.
    expect(out.atoms[1].y).toBeLessThan(out.atoms[0].y);
  });

  it("survives a molecule with no extent on an axis", () => {
    // A linear molecule has zero y-span; an unguarded divide makes NaN.
    const coords = [
      { x: 0, y: 0 },
      { x: 1, y: 0 },
      { x: 2, y: 0 },
    ];
    const out = layoutMonomer(atoms("A", "B", "C"), [], coords, BOX);
    expect(out.atoms.every((a) => Number.isFinite(a.x) && Number.isFinite(a.y))).toBe(true);
  });

  it("survives a single atom", () => {
    const out = layoutMonomer(atoms("NZ"), [], [{ x: 3, y: 4 }], BOX);
    expect(out.atoms).toHaveLength(1);
    expect(Number.isFinite(out.atoms[0].x)).toBe(true);
    expect(Number.isFinite(out.atoms[0].y)).toBe(true);
  });

  it("drops a bond whose atom was not laid out", () => {
    const bonds: MonomerBond[] = [
      { atom1: "A", atom2: "B", type: "single" },
      { atom1: "A", atom2: "GHOST", type: "single" },
    ];
    const out = layoutMonomer(atoms("A", "B"), bonds, [{ x: 0, y: 0 }, { x: 1, y: 1 }], BOX);
    expect(out.bonds).toHaveLength(1);
  });

  it("returns an empty layout rather than throwing when given nothing", () => {
    expect(layoutMonomer([], [], [], BOX).atoms).toEqual([]);
  });

  it("ignores coordinates beyond the atoms it was given", () => {
    // RDKit adding an atom (it should not, but the contract is order-based)
    // must not produce a phantom.
    const out = layoutMonomer(atoms("A"), [], [{ x: 0, y: 0 }, { x: 5, y: 5 }], BOX);
    expect(out.atoms).toHaveLength(1);
  });
});

describe("bondLines", () => {
  const bond = (type: string): PlacedBond => ({
    atom1: "A",
    atom2: "B",
    type,
    x1: 0,
    y1: 0,
    x2: 10,
    y2: 0,
  });

  it("draws one line for a single bond", () => {
    expect(bondLines(bond("single"), 4)).toHaveLength(1);
  });

  it("draws two parallel lines for a double bond, offset perpendicular", () => {
    const lines = bondLines(bond("double"), 4);
    expect(lines).toHaveLength(2);
    // The bond runs along x, so the offset must be in y, and symmetric.
    expect(lines[0].y1).toBeCloseTo(-2);
    expect(lines[1].y1).toBeCloseTo(2);
    expect(lines[0].x1).toBeCloseTo(lines[1].x1);
  });

  it("draws three for a triple, one of them down the centre", () => {
    const lines = bondLines(bond("triple"), 4);
    expect(lines).toHaveLength(3);
    expect(lines.some((l) => Math.abs(l.y1) < 1e-9)).toBe(true);
  });

  it("does not divide by zero for two atoms in the same place", () => {
    const degenerate: PlacedBond = { ...bond("double"), x2: 0, y2: 0 };
    const lines = bondLines(degenerate, 4);
    expect(lines.every((l) => Number.isFinite(l.x1) && Number.isFinite(l.y1))).toBe(true);
  });
});

describe("multiplicity", () => {
  it("draws delocalised bonds as double, not as single", () => {
    // CCP4 dictionaries use both spellings; drawing them as single makes an
    // aromatic ring look saturated.
    expect(multiplicity("aromatic")).toBe(2);
    expect(multiplicity("deloc")).toBe(2);
  });

  it("is case-insensitive and falls back to single", () => {
    expect(multiplicity("DOUBLE")).toBe(2);
    expect(multiplicity("metal")).toBe(1);
    expect(multiplicity("nonsense")).toBe(1);
  });
});

describe("bondKey", () => {
  it("identifies a bond regardless of the order its atoms are given", () => {
    expect(bondKey("C", "O")).toBe(bondKey("O", "C"));
    expect(bondKey("C", "O")).not.toBe(bondKey("C", "N"));
  });
});
