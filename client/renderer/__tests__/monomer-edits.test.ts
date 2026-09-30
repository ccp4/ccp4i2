/**
 * A monomer's edits as data, and the valence count that advises on them.
 *
 * The GLU and LYS shapes are the CCP4 monomer library's (heavy atoms, their
 * charges and hydrogen counts, Kekule bonds). The valence cases are the ones
 * AceDRG was seen to accept or refuse (ccp4-20260702).
 */
import { describe, it, expect } from "vitest";
import {
  NO_EDITS, isEmpty, toggleDelete, setBondOrder, setCharge, editedOrder, editedCharge,
  dictionaryOrder, valenceAdvisories,
} from "../lib/monomer-edits";
import type { MonomerAtomDetail, MonomerBond } from "../lib/monomer-molblock";

const atom = (name: string, element: string, charge = 0, hydrogens = 0): MonomerAtomDetail =>
  ({ name, element, charge, hydrogens });
const bond = (atom1: string, atom2: string, type = "single"): MonomerBond => ({ atom1, atom2, type });

const GLU = {
  atoms: [
    atom("N", "N", 1, 3), atom("CA", "C", 0, 1), atom("C", "C"), atom("O", "O"),
    atom("CB", "C", 0, 2), atom("CG", "C", 0, 2), atom("CD", "C"),
    atom("OE1", "O"), atom("OE2", "O", -1), atom("OXT", "O", -1),
  ],
  bonds: [
    bond("N", "CA"), bond("CA", "C"), bond("CA", "CB"), bond("C", "O", "double"),
    bond("C", "OXT"), bond("CB", "CG"), bond("CG", "CD"), bond("CD", "OE1", "double"),
    bond("CD", "OE2"),
  ],
};

const LYS = {
  atoms: [
    atom("N", "N", 1, 3), atom("CA", "C", 0, 1), atom("C", "C"), atom("O", "O"),
    atom("CB", "C", 0, 2), atom("CG", "C", 0, 2), atom("CD", "C", 0, 2),
    atom("CE", "C", 0, 2), atom("NZ", "N", 1, 3), atom("OXT", "O", -1),
  ],
  bonds: [
    bond("N", "CA"), bond("CA", "C"), bond("CA", "CB"), bond("C", "O", "double"),
    bond("C", "OXT"), bond("CB", "CG"), bond("CG", "CD"), bond("CD", "CE"), bond("CE", "NZ"),
  ],
};

describe("editing", () => {
  it("toggles a deletion on and off", () => {
    const once = toggleDelete(NO_EDITS, "OE2");
    expect(once.deletes).toEqual(["OE2"]);
    expect(toggleDelete(once, "OE2")).toEqual(NO_EDITS);
  });

  it("drops the edits that touch an atom when it is deleted", () => {
    let edits = setBondOrder(NO_EDITS, "C4A", "O4A", "SINGLE", "DOUBLE");
    edits = setCharge(edits, "O4A", -1, 0);
    edits = toggleDelete(edits, "O4A");
    expect(edits).toEqual({ deletes: ["O4A"], bondOrders: [], charges: [] });
  });

  it("keeps one order per bond, whichever way round it is named", () => {
    let edits = setBondOrder(NO_EDITS, "CD", "OE1", "SINGLE", "DOUBLE");
    edits = setBondOrder(edits, "OE1", "CD", "TRIPLE", "DOUBLE");
    expect(edits.bondOrders).toEqual([{ atom1: "OE1", atom2: "CD", order: "TRIPLE" }]);
    expect(editedOrder(edits, "CD", "OE1")).toBe("TRIPLE");
  });

  it("stores nothing for the dictionary's own order or charge", () => {
    const edits = setBondOrder(NO_EDITS, "CD", "OE1", "SINGLE", "DOUBLE");
    expect(isEmpty(setBondOrder(edits, "CD", "OE1", "DOUBLE", "DOUBLE"))).toBe(true);
    expect(isEmpty(setBondOrder(edits, "CD", "OE1", null, "DOUBLE"))).toBe(true);
    expect(isEmpty(setCharge(setCharge(NO_EDITS, "NZ", 0, 1), "NZ", 1, 1))).toBe(true);
  });

  it("keeps one charge per atom", () => {
    const edits = setCharge(setCharge(NO_EDITS, "NZ", 0, 1), "NZ", -1, 1);
    expect(edits.charges).toEqual([{ atom: "NZ", charge: -1 }]);
    expect(editedCharge(edits, "NZ")).toBe(-1);
  });

  it("names only the orders it can edit", () => {
    expect(dictionaryOrder("double")).toBe("DOUBLE");
    expect(dictionaryOrder("aromatic")).toBeNull();
  });
});

describe("valence advisories", () => {
  const advise = (monomer: typeof GLU, edits = NO_EDITS, link?: string) =>
    valenceAdvisories(monomer.atoms, monomer.bonds, edits, link);

  it("find nothing to say about the isopeptide AceDRG accepted", () => {
    // GLU CD linked, OE2 deleted: CG + =OE1 + link = 4.
    expect(advise(GLU, toggleDelete(NO_EDITS, "OE2"), "CD")).toEqual([]);
    // LYS NZ linked: AceDRG takes a hydrogen off it, so no complaint.
    expect(advise(LYS, NO_EDITS, "NZ")).toEqual([]);
  });

  it("flag the deletion AceDRG refused, where RDKit would add a hydrogen", () => {
    // Deleting OXT leaves C with CA and =O: valence 3.
    const advisories = advise(GLU, toggleDelete(NO_EDITS, "OXT"));
    expect(advisories).toHaveLength(1);
    expect(advisories[0].atom).toBe("C");
    expect(advisories[0].message).toBe("C would have valence 3; C usually has 4");
  });

  it("flag a linking atom that keeps all its bonds", () => {
    // Linking GLU CD without deleting OE2: CG + =OE1 + OE2 + link = 5.
    const [advisory] = advise(GLU, NO_EDITS, "CD");
    expect(advisory.message).toBe("CD would have valence 5 with the link; C usually has 4");
  });

  it("follow a changed bond order", () => {
    const edits = setBondOrder(NO_EDITS, "CD", "OE1", "TRIPLE", "DOUBLE");
    expect(advise(GLU, edits).map((a) => a.atom).sort()).toEqual(["CD", "OE1"]);
  });

  it("judge a changed charge against the new charge", () => {
    // O4A-like: O with one single bond and no H, made -1: valence 1 is right.
    const alcohol = { atoms: [atom("C1", "C", 0, 3), atom("O1", "O", 0, 0)], bonds: [bond("C1", "O1")] };
    expect(advise(alcohol, setCharge(NO_EDITS, "O1", -1, 0))).toEqual([]);
    expect(advise(alcohol, setCharge(NO_EDITS, "O1", 1, 0))).not.toEqual([]);
  });

  it("leave alone what it cannot count", () => {
    const ring = {
      atoms: [atom("C1", "C", 0, 1), atom("C2", "C", 0, 1), atom("O1", "O")],
      bonds: [bond("C1", "C2", "aromatic"), bond("C1", "O1")],
    };
    expect(advise(ring, toggleDelete(NO_EDITS, "O1"))).toEqual([]);
    const metal = { atoms: [atom("FE", "Fe"), atom("O1", "O")], bonds: [bond("FE", "O1")] };
    expect(advise(metal, toggleDelete(NO_EDITS, "O1"))).toEqual([]);
  });
});
