/**
 * The edits MakeLink makes to one monomer, as data.
 *
 * Mirrors the task's DELETE_ATOMS_n / BOND_ORDERS_n / CHARGES_n: what the
 * modified monomer should be, not the AceDRG commands that make it (the
 * server writes those, server/ccp4i2/pipelines/MakeLink/script/
 * link_instruction.py). Every function here returns a new MonomerEdits, and
 * keeps it free of the contradictions the server would refuse: deleting an
 * atom drops any edit that touches it, and an edit that restores the
 * dictionary's own value is removed rather than stored.
 *
 * valenceAdvisories() is a count, not chemistry. RDKit cannot do this job:
 * it rejects an over-valent atom without saying which, and quietly adds
 * hydrogens to an under-valent one -- the very case AceDRG refuses ("atom C
 * in monomer GLU has a total valence of 3"). AceDRG stays the authority, so
 * what this finds is advice, never a block.
 */
import type { MonomerAtomDetail, MonomerBond } from "./monomer-molblock";

export type BondOrder = "SINGLE" | "DOUBLE" | "TRIPLE";
export const BOND_ORDERS: BondOrder[] = ["SINGLE", "DOUBLE", "TRIPLE"];

export interface BondOrderEdit {
  atom1: string;
  atom2: string;
  order: BondOrder;
}

export interface ChargeEdit {
  atom: string;
  charge: number;
}

export interface MonomerEdits {
  deletes: string[];
  bondOrders: BondOrderEdit[];
  charges: ChargeEdit[];
}

export const NO_EDITS: MonomerEdits = { deletes: [], bondOrders: [], charges: [] };

export function isEmpty(edits: MonomerEdits): boolean {
  return !edits.deletes.length && !edits.bondOrders.length && !edits.charges.length;
}

const sameBond = (a1: string, a2: string, b1: string, b2: string) =>
  (a1 === b1 && a2 === b2) || (a1 === b2 && a2 === b1);

/** A dictionary bond's order as one of ours; aromatic and the like have none. */
export function dictionaryOrder(type: string | undefined): BondOrder | null {
  const upper = String(type ?? "").toUpperCase();
  return (BOND_ORDERS as string[]).includes(upper) ? (upper as BondOrder) : null;
}

/** Delete an atom, or restore it if it was deleted. */
export function toggleDelete(edits: MonomerEdits, atom: string): MonomerEdits {
  if (edits.deletes.includes(atom)) {
    return { ...edits, deletes: edits.deletes.filter((name) => name !== atom) };
  }
  // A deleted atom cannot also have its bonds or charge changed.
  return {
    deletes: [...edits.deletes, atom],
    bondOrders: edits.bondOrders.filter((b) => b.atom1 !== atom && b.atom2 !== atom),
    charges: edits.charges.filter((c) => c.atom !== atom),
  };
}

/**
 * Give a bond an order. `order` null, or equal to the dictionary's own, means
 * "as in the dictionary", which removes any edit rather than storing one.
 */
export function setBondOrder(
  edits: MonomerEdits,
  atom1: string,
  atom2: string,
  order: BondOrder | null,
  dictionary: BondOrder | null
): MonomerEdits {
  const others = edits.bondOrders.filter((b) => !sameBond(b.atom1, b.atom2, atom1, atom2));
  const bondOrders = order && order !== dictionary ? [...others, { atom1, atom2, order }] : others;
  return { ...edits, bondOrders };
}

/** Give an atom a formal charge; null, or the dictionary's own, removes the edit. */
export function setCharge(
  edits: MonomerEdits,
  atom: string,
  charge: number | null,
  dictionary: number
): MonomerEdits {
  const others = edits.charges.filter((c) => c.atom !== atom);
  const charges = charge !== null && charge !== dictionary ? [...others, { atom, charge }] : others;
  return { ...edits, charges };
}

export function editedOrder(edits: MonomerEdits, atom1: string, atom2: string): BondOrder | null {
  return edits.bondOrders.find((b) => sameBond(b.atom1, b.atom2, atom1, atom2))?.order ?? null;
}

export function editedCharge(edits: MonomerEdits, atom: string): number | null {
  return edits.charges.find((c) => c.atom === atom)?.charge ?? null;
}

// The usual valences by element and formal charge. Elements not listed
// (metals, selenium, ...) are not judged.
const VALENCES: Record<string, Record<number, number[]>> = {
  C: { 0: [4], 1: [3], [-1]: [3] },
  N: { 0: [3], 1: [4], [-1]: [2] },
  O: { 0: [2], 1: [3], [-1]: [1] },
  S: { 0: [2, 4, 6], 1: [3], [-1]: [1] },
  P: { 0: [3, 5], 1: [4] },
  B: { 0: [3], [-1]: [4] },
  F: { 0: [1] },
  CL: { 0: [1] },
  BR: { 0: [1] },
  I: { 0: [1] },
};

const ORDER_VALUE: Record<BondOrder, number> = { SINGLE: 1, DOUBLE: 2, TRIPLE: 3 };

export interface ValenceAdvisory {
  atom: string;
  message: string;
}

/**
 * Atoms whose valence would be unusual once the edits and the link are made.
 *
 * Only atoms an edit or the link touches are judged. Hydrogens are the
 * dictionary's; for the linking atom and any atom whose charge changes,
 * AceDRG may remove hydrogens itself, so there only too many bonds to heavy
 * atoms counts against it. An atom with a bond order this cannot count
 * (aromatic, delocalised, metal) is not judged at all.
 */
export function valenceAdvisories(
  atomDetails: MonomerAtomDetail[],
  bonds: MonomerBond[],
  edits: MonomerEdits,
  linkAtom?: string | null,
  linkOrder: BondOrder = "SINGLE"
): ValenceAdvisory[] {
  const deleted = new Set(edits.deletes);
  const touched = new Set<string>();
  if (linkAtom) touched.add(linkAtom);
  for (const bond of bonds) {
    if (deleted.has(bond.atom1) && !deleted.has(bond.atom2)) touched.add(bond.atom2);
    if (deleted.has(bond.atom2) && !deleted.has(bond.atom1)) touched.add(bond.atom1);
  }
  for (const b of edits.bondOrders) {
    touched.add(b.atom1);
    touched.add(b.atom2);
  }
  for (const c of edits.charges) touched.add(c.atom);

  const advisories: ValenceAdvisory[] = [];
  for (const atom of atomDetails) {
    if (!touched.has(atom.name) || deleted.has(atom.name)) continue;
    const element = String(atom.element).toUpperCase();
    const charge = editedCharge(edits, atom.name) ?? atom.charge;
    const allowed = VALENCES[element]?.[charge];
    if (!allowed) continue;

    let heavy = 0;
    let countable = true;
    for (const bond of bonds) {
      if (bond.atom1 !== atom.name && bond.atom2 !== atom.name) continue;
      const other = bond.atom1 === atom.name ? bond.atom2 : bond.atom1;
      if (deleted.has(other)) continue;
      const order = editedOrder(edits, bond.atom1, bond.atom2) ?? dictionaryOrder(bond.type);
      if (!order) {
        countable = false;
        break;
      }
      heavy += ORDER_VALUE[order];
    }
    if (!countable) continue;
    const isLink = atom.name === linkAtom;
    if (isLink) heavy += ORDER_VALUE[linkOrder];

    const hydrogens = atom.hydrogens ?? 0;
    const flexible = isLink || editedCharge(edits, atom.name) !== null;
    const low = flexible ? heavy : heavy + hydrogens;
    const high = heavy + hydrogens;
    if (allowed.some((v) => v >= low && v <= high)) continue;

    const valence = flexible && heavy > Math.max(...allowed) ? heavy : high;
    const usual = allowed.join(" or ");
    const chargeText = charge ? ` with charge ${charge > 0 ? "+" : ""}${charge}` : "";
    advisories.push({
      atom: atom.name,
      message: `${atom.name} would have valence ${valence}${isLink ? " with the link" : ""}; ${atom.element}${chargeText} usually has ${usual}`,
    });
  }
  return advisories;
}
