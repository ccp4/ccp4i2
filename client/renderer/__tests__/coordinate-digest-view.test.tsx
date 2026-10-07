/**
 * #681: the coordinate-file Digest was raw JSON; it is now a chain table.
 * A chain's range is its polymer range where it has one, since the
 * whole-chain range runs over waters numbered after the polymer.
 */
import React from "react";
import { describe, it, expect } from "vitest";
import { render, screen, within } from "@testing-library/react";
import {
  chainDigestRows,
  CoordinateDigestView,
  isCoordinateDigest,
} from "../components/coordinate-digest-view";

// Chain details as the server reports them for 1jst (demo_data) plus a
// water-only chain, abbreviated.
const DIGEST = {
  spaceGroup: "P 21 21 21",
  cell: { a: 73.1, b: 134.5, c: 148.2, alpha: 90, beta: 90, gamma: 90 },
  composition: {
    nModels: 1,
    nAtoms: 9134,
    chainDetails: [
      {
        id: "B", type: "protein", nResidues: 282, firstRes: "175", lastRes: "97",
        nPolymer: 258, polymerFirstRes: "175", polymerLastRes: "432", nWater: 24, nOther: 0,
      },
      {
        id: "W", type: "solvent", nResidues: 134, firstRes: "1", lastRes: "134",
        nPolymer: 0, polymerFirstRes: "", polymerLastRes: "", nWater: 134, nOther: 0,
      },
      {
        id: "", type: "nucleic", nResidues: 12, firstRes: "1", lastRes: "12",
        nPolymer: 12, polymerFirstRes: "1", polymerLastRes: "12", nWater: 0, nOther: 0,
      },
    ],
    ligands: [{ chain: "A", name: "ATP", seqNum: 401, atomCount: 31 }],
    residueNameCounts: { HOH: 158, ALA: 20 },
  },
};

describe("chainDigestRows", () => {
  it("uses the polymer count and range, not the whole chain's", () => {
    const [b, w, blank] = chainDigestRows(DIGEST);
    expect(b).toEqual({
      chain: "B", type: "Protein", residues: 258, range: "175–432", waters: 24, other: 0,
    });
    expect(w).toMatchObject({ chain: "W", type: "Water", residues: 134, range: "1–134" });
    expect(blank).toMatchObject({ chain: "(blank)", type: "Nucleic acid", residues: 12 });
  });

  it("falls back to the whole-chain range for a digest without the split", () => {
    const [row] = chainDigestRows({
      composition: {
        chainDetails: [{ id: "A", type: "protein", nResidues: 300, firstRes: "1", lastRes: "300" }],
      },
    });
    expect(row).toMatchObject({ residues: 300, range: "1–300", waters: null });
  });

  it("recognises a coordinate digest only by its chain details", () => {
    expect(isCoordinateDigest(DIGEST)).toBe(true);
    expect(isCoordinateDigest({ cell: {} })).toBe(false);
    expect(isCoordinateDigest(null)).toBe(false);
  });
});

describe("CoordinateDigestView", () => {
  it("renders a chain table and the ligands", () => {
    render(<CoordinateDigestView digest={DIGEST} />);
    const chains = screen.getByRole("table", { name: "Chains" });
    expect(within(chains).getByText("175–432")).toBeTruthy();
    const ligands = screen.getByRole("table", { name: "Ligands" });
    expect(within(ligands).getByText("ATP")).toBeTruthy();
    expect(screen.getByText(/Space group P 21 21 21/)).toBeTruthy();
  });
});
