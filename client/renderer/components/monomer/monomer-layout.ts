/**
 * Laying a monomer out in 2D, and the geometry the drawing needs.
 *
 * Kept apart from the component so it can be tested without a DOM or a
 * WASM module: everything here is arithmetic over the coordinates RDKit
 * returned.
 */
import type { MonomerAtomDetail, MonomerBond } from "../../lib/monomer-molblock";

export interface PlacedAtom {
  name: string;
  element: string;
  charge: number;
  x: number;
  y: number;
}

export interface PlacedBond {
  atom1: string;
  atom2: string;
  type: string;
  x1: number;
  y1: number;
  x2: number;
  y2: number;
}

export interface MonomerLayout {
  atoms: PlacedAtom[];
  bonds: PlacedBond[];
  width: number;
  height: number;
}

/**
 * Fit raw 2D coordinates into a box.
 *
 * RDKit works in bond-length units and puts the origin wherever it likes,
 * with y increasing upwards; SVG wants pixels in a known box with y
 * increasing downwards. A single uniform scale keeps bond lengths equal,
 * which is the whole point of a chemical depiction -- fitting width and
 * height independently would stretch the molecule.
 */
export function layoutMonomer(
  atomDetails: MonomerAtomDetail[],
  bonds: MonomerBond[],
  coordinates: { x: number; y: number }[],
  box: { width: number; height: number; padding: number }
): MonomerLayout {
  const usable = Math.min(atomDetails.length, coordinates.length);
  if (usable === 0) {
    return { atoms: [], bonds: [], width: box.width, height: box.height };
  }

  const xs = coordinates.slice(0, usable).map((c) => c.x);
  const ys = coordinates.slice(0, usable).map((c) => c.y);
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);

  const spanX = maxX - minX;
  const spanY = maxY - minY;
  const usableWidth = Math.max(1, box.width - 2 * box.padding);
  const usableHeight = Math.max(1, box.height - 2 * box.padding);
  // A single atom, or a perfectly linear molecule, has no span on one axis.
  // Guard both, or the scale is Infinity and every atom lands on NaN.
  const scale = Math.min(
    spanX > 1e-6 ? usableWidth / spanX : Infinity,
    spanY > 1e-6 ? usableHeight / spanY : Infinity
  );
  const finiteScale = Number.isFinite(scale) ? scale : 1;

  // Centre whatever is left over, so a wide flat molecule sits in the middle
  // rather than hugging the top.
  const offsetX = box.padding + (usableWidth - spanX * finiteScale) / 2;
  const offsetY = box.padding + (usableHeight - spanY * finiteScale) / 2;

  const atoms: PlacedAtom[] = [];
  for (let i = 0; i < usable; i += 1) {
    const detail = atomDetails[i];
    atoms.push({
      name: detail.name,
      element: detail.element,
      charge: detail.charge ?? 0,
      x: offsetX + (coordinates[i].x - minX) * finiteScale,
      // SVG y grows downwards; chemistry drawings do not.
      y: offsetY + (maxY - coordinates[i].y) * finiteScale,
    });
  }

  const byName = new Map(atoms.map((a) => [a.name, a]));
  const placedBonds: PlacedBond[] = [];
  for (const bond of bonds) {
    const a = byName.get(bond.atom1);
    const b = byName.get(bond.atom2);
    if (!a || !b) continue; // a bond to an atom that was not laid out
    placedBonds.push({
      atom1: bond.atom1,
      atom2: bond.atom2,
      type: bond.type,
      x1: a.x,
      y1: a.y,
      x2: b.x,
      y2: b.y,
    });
  }

  return { atoms, bonds: placedBonds, width: box.width, height: box.height };
}

/**
 * The parallel lines of a multiple bond, offset perpendicular to it.
 * Returns one line for a single bond, two for a double, three for a triple.
 */
export function bondLines(
  bond: PlacedBond,
  separation: number
): { x1: number; y1: number; x2: number; y2: number }[] {
  const order = multiplicity(bond.type);
  const line = { x1: bond.x1, y1: bond.y1, x2: bond.x2, y2: bond.y2 };
  if (order <= 1) return [line];

  const dx = bond.x2 - bond.x1;
  const dy = bond.y2 - bond.y1;
  const length = Math.hypot(dx, dy);
  // Two atoms on top of each other have no direction to offset along.
  if (length < 1e-6) return [line];
  const nx = -dy / length;
  const ny = dx / length;

  // Two lines straddle the centre; three put one on it.
  const offsets =
    order === 2 ? [-separation / 2, separation / 2] : [-separation, 0, separation];
  return offsets.map((offset) => ({
    x1: bond.x1 + nx * offset,
    y1: bond.y1 + ny * offset,
    x2: bond.x2 + nx * offset,
    y2: bond.y2 + ny * offset,
  }));
}

/** How many lines a bond type draws as. Aromatic and deloc draw as two. */
export function multiplicity(type: string): number {
  switch (String(type).toLowerCase()) {
    case "double":
    case "aromatic":
    case "deloc":
      return 2;
    case "triple":
      return 3;
    default:
      return 1;
  }
}

/**
 * A stable identity for a bond, order-independent: C-O and O-C are one bond.
 * The separator is a character CIF atom names cannot contain.
 */
export function bondKey(atom1: string, atom2: string): string {
  return [atom1, atom2].sort().join("|");
}
