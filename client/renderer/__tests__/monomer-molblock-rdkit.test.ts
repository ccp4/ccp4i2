/**
 * The molblock we build must be one RDKit will actually accept.
 *
 * The hand-written unit tests passed while RDKit rejected two monomers out
 * of four: the M CHG fields were a character too narrow, which is invisible
 * to any assertion that does not put the result in front of a real parser.
 * So this puts real dictionary-shaped input through the real library.
 */
import { describe, it, expect, beforeAll } from "vitest";
import { buildMolblock, readCoordinates, type MonomerAtomDetail, type MonomerBond } from "../lib/monomer-molblock";

interface Monomer {
  atom_details: MonomerAtomDetail[];
  bonds: MonomerBond[];
}

// Shapes taken from the CCP4 monomer library: a charged amino acid (both
// signs at once), an aromatic/fused ring, and a phosphate chain.
const MONOMERS: Record<string, Monomer> = {
  // The real CCP4 LYS, all ten heavy atoms. The size matters: its charged
  // OXT is atom 10, and a two-digit index is what exposes a mis-sized M CHG
  // field. A trimmed six-atom version parses happily even when broken.
  LYS: {
    atom_details: [
      { name: "N", element: "N", charge: 1 }, { name: "CA", element: "C", charge: 0 },
      { name: "C", element: "C", charge: 0 }, { name: "O", element: "O", charge: 0 },
      { name: "CB", element: "C", charge: 0 }, { name: "CG", element: "C", charge: 0 },
      { name: "CD", element: "C", charge: 0 }, { name: "CE", element: "C", charge: 0 },
      { name: "NZ", element: "N", charge: 1 }, { name: "OXT", element: "O", charge: -1 },
    ],
    bonds: [
      { atom1: "N", atom2: "CA", type: "single" }, { atom1: "CA", atom2: "C", type: "single" },
      { atom1: "CA", atom2: "CB", type: "single" }, { atom1: "C", atom2: "O", type: "double" },
      { atom1: "C", atom2: "OXT", type: "single" }, { atom1: "CB", atom2: "CG", type: "single" },
      { atom1: "CG", atom2: "CD", type: "single" }, { atom1: "CD", atom2: "CE", type: "single" },
      { atom1: "CE", atom2: "NZ", type: "single" },
    ],
  },
  BENZENE_LIKE: {
    atom_details: Array.from({ length: 6 }, (_, i) => ({
      name: `C${i + 1}`,
      element: "C",
      charge: 0,
    })),
    bonds: Array.from({ length: 6 }, (_, i) => ({
      atom1: `C${i + 1}`,
      atom2: `C${((i + 1) % 6) + 1}`,
      type: "aromatic",
    })),
  },
  PHOSPHATE: {
    atom_details: [
      { name: "P", element: "P", charge: 0 },
      { name: "O1P", element: "O", charge: -1 },
      { name: "O2P", element: "O", charge: -1 },
      { name: "O3P", element: "O", charge: 0 },
      { name: "O5'", element: "O", charge: 0 },
    ],
    bonds: [
      { atom1: "P", atom2: "O1P", type: "single" },
      { atom1: "P", atom2: "O2P", type: "single" },
      { atom1: "P", atom2: "O3P", type: "double" },
      { atom1: "P", atom2: "O5'", type: "single" },
    ],
  },
};

let RDKit: any = null;

beforeAll(async () => {
  try {
    const mod = await import("@rdkit/rdkit/dist/RDKit_minimal.js");
    RDKit = await (mod.default as any)();
    RDKit.prefer_coordgen(true);
  } catch {
    RDKit = null; // reported per-test below
  }
}, 60000);

describe("the molblock is one RDKit accepts", () => {
  for (const [code, monomer] of Object.entries(MONOMERS)) {
    it(`${code} parses and lays out`, () => {
      if (!RDKit) return expect.unreachable("RDKit WASM failed to load");
      const mol = RDKit.get_mol(buildMolblock(monomer.atom_details, monomer.bonds, code));
      expect(mol, `${code}: RDKit rejected the molblock`).toBeTruthy();
      try {
        mol.set_new_coords(true);
        const coords = readCoordinates(mol.get_molblock());
        // One coordinate per atom, in the order we supplied: that ordering
        // is what lets the picker label atoms by dictionary name.
        expect(coords).toHaveLength(monomer.atom_details.length);
        expect(coords.every((c) => Number.isFinite(c.x) && Number.isFinite(c.y))).toBe(true);
        // A real layout, not everything stacked at the origin.
        const distinct = new Set(coords.map((c) => `${c.x.toFixed(3)},${c.y.toFixed(3)}`));
        expect(distinct.size).toBe(monomer.atom_details.length);
      } finally {
        mol.delete();
      }
    }, 30000);
  }
});
