/**
 * A monomer dictionary's atoms and bonds as a V2000 molblock.
 *
 * This exists so a depiction can be keyed by dictionary atom name. RDKit is
 * asked only to lay the molecule out in 2D; the drawing is ours, and every
 * element it emits carries the name the AceDRG tasks speak ("NZ", "C4A").
 *
 * The contract that makes that work is atom order: molblock atom i is
 * `atoms[i]`, so coordinates come back in the order they went in and no
 * index has to be mapped back to a name. The server guarantees the same
 * order between `atoms` and `atom_details`.
 *
 * Note this deliberately does NOT go via SMILES. The server can already
 * produce molblocks that way (`generate_all_molblocks`), but SMILES loses
 * the atom names and reorders the atoms, which is exactly what we need kept.
 */

export interface MonomerBond {
  atom1: string;
  atom2: string;
  type: string;
}

export interface MonomerAtomDetail {
  name: string;
  element: string;
  charge: number;
  /** Hydrogens the dictionary bonds to this atom; not drawn, but counted. */
  hydrogens?: number;
}

/** Molfile bond-block codes. */
const BOND_ORDER: Record<string, number> = {
  single: 1,
  double: 2,
  triple: 3,
  // CCP4 dictionaries also use these two. Both mean "delocalised" to a
  // depiction, and 4 is the molfile's aromatic code; without this they would
  // fall back to single and the ring would be drawn wrong.
  aromatic: 4,
  deloc: 4,
  // A metal coordination bond is not a covalent bond. Drawn as single: the
  // molfile has no better code, and 0 makes RDKit reject the molecule.
  metal: 1,
};

/**
 * Molfiles encode charge as a code, not a number: 1 => +3, 3 => +1,
 * 5 => -1, 7 => -3, i.e. `4 - charge`, with 0 meaning uncharged and 4
 * reserved for a doublet radical. Anything outside +/-3 has no code and is
 * dropped here rather than mis-encoded -- the M CHG block below carries the
 * real value, and that is what RDKit reads.
 */
function chargeCode(charge: number): number {
  if (!Number.isFinite(charge) || charge === 0) return 0;
  const rounded = Math.round(charge);
  if (rounded > 3 || rounded < -3) return 0;
  return 4 - rounded;
}

/**
 * An element symbol as a molfile spells it: "Cl", "Br", where CIF
 * dictionaries write "CL", "BR". RDKit's molfile reader happens to accept
 * either (checked against the bundled MinimalLib), and the server already
 * sends gemmi's normalised symbol; this keeps the molblock canonical for any
 * reader that is stricter, and is the one place every depiction passes
 * through. (In SMILES the case is not cosmetic -- "CL" is not chlorine --
 * which is one more reason the picker never goes that way.)
 */
export function elementSymbol(element: string | undefined): string {
  const symbol = (element || "C").trim().slice(0, 3);
  return symbol.charAt(0).toUpperCase() + symbol.slice(1).toLowerCase();
}

function padLeft(text: string, width: number): string {
  return text.length >= width ? text : " ".repeat(width - text.length) + text;
}

function counts(value: number): string {
  return padLeft(String(value), 3);
}

/** Property-block fields (M CHG and friends) are four characters wide. */
function field4(value: number): string {
  return padLeft(String(value), 4);
}

function coordinate(value: number): string {
  return padLeft(value.toFixed(4), 10);
}

/**
 * Build a V2000 molblock. Coordinates are all zero: RDKit is asked to
 * generate them, and passing the dictionary's own x/y/z would be worse than
 * useless -- those are ideal *3D* coordinates, which overlap badly for rings
 * when flattened.
 */
export function buildMolblock(
  atomDetails: MonomerAtomDetail[],
  bonds: MonomerBond[],
  title = ""
): string {
  const index = new Map<string, number>();
  atomDetails.forEach((atom, i) => index.set(atom.name, i + 1)); // molfiles are 1-based

  // A bond naming an atom that is not in the list (a hydrogen the server
  // dropped, a typo in a hand-edited dictionary) is skipped, not fatal.
  const usable = bonds.filter((b) => index.has(b.atom1) && index.has(b.atom2));

  const lines: string[] = [
    title,
    "  ccp4i2  2D",
    "",
    `${counts(atomDetails.length)}${counts(usable.length)}  0  0  0  0  0  0  0  0999 V2000`,
  ];

  for (const atom of atomDetails) {
    const symbol = elementSymbol(atom.element);
    lines.push(
      `${coordinate(0)}${coordinate(0)}${coordinate(0)} ${symbol.padEnd(3)} 0` +
        `${counts(chargeCode(atom.charge ?? 0))}`
    );
  }

  for (const bond of usable) {
    const order = BOND_ORDER[String(bond.type).toLowerCase()] ?? 1;
    lines.push(
      `${counts(index.get(bond.atom1)!)}${counts(index.get(bond.atom2)!)}${counts(order)}  0`
    );
  }

  // An M CHG block: the atom-line charge column is legacy and RDKit prefers
  // this. Only non-zero charges appear, eight per line at most.
  //
  // Its entries are FOUR characters wide, not three like every count on the
  // lines above (the spec writes each as " aaa"). Getting that wrong does not
  // fail loudly -- RDKit rejects the whole molecule, and only for some inputs,
  // so LYS and ATP came back null while PLP and TRP parsed.
  const charged = atomDetails
    .map((atom, i) => [i + 1, Math.round(atom.charge ?? 0)] as const)
    .filter(([, charge]) => charge !== 0);
  for (let i = 0; i < charged.length; i += 8) {
    const chunk = charged.slice(i, i + 8);
    lines.push(
      `M  CHG${counts(chunk.length)}` +
        chunk.map(([idx, charge]) => `${field4(idx)}${field4(charge)}`).join("")
    );
  }

  lines.push("M  END");
  return lines.join("\n") + "\n";
}

/** 2D coordinates read back out of a molblock, in atom order. */
export function readCoordinates(molblock: string): { x: number; y: number }[] {
  const lines = molblock.split("\n");
  const countsLine = lines[3];
  if (!countsLine) return [];
  const atomCount = parseInt(countsLine.slice(0, 3), 10);
  if (!Number.isFinite(atomCount) || atomCount <= 0) return [];
  const out: { x: number; y: number }[] = [];
  for (let i = 0; i < atomCount; i += 1) {
    const line = lines[4 + i];
    if (!line) break;
    out.push({
      x: parseFloat(line.slice(0, 10)),
      y: parseFloat(line.slice(10, 20)),
    });
  }
  return out;
}
