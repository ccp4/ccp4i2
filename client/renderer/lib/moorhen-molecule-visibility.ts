/**
 * Hiding and re-showing a whole molecule from a side-panel row.
 *
 * Moorhen's own MoleculeCard restores a hidden molecule with
 * `rep.interfaceOption.visible`, and every panel in this app copied that.
 * But that flag is only ever set true by Moorhen's `RepresentationBuilder`
 * and by the representation chips in its Models drawer. Every representation
 * *this* app creates goes through the deprecated `molecule.addRepresentation`,
 * which leaves `interfaceOption.visible` as `undefined` — so the restore
 * branch never fired, and a molecule toggled off could not be toggled back
 * on. The eye icon worked once, in one direction.
 *
 * Showing everything instead is not the fix. The scene resolver deliberately
 * hides each molecule's loader-default ribbon so that an applied scene alone
 * owns the look; a blanket `show()` would bring one ghost ribbon per molecule
 * back the first time someone toggled a ligand off and on.
 *
 * So the toggle remembers: at the moment of hiding, which representations
 * were actually drawn, and on showing, exactly those come back.
 */
import { moorhen } from "moorhen/types/moorhen";

/**
 * What was drawn on each hidden molecule, keyed by molNo. Callers hold one of
 * these per panel (a `useRef`), so it lives as long as the loaded molecules
 * do and no longer.
 */
export type MoleculeVisibilityMemory = Map<number, Set<string>>;

/** Moorhen's representation objects, in the shape this module needs. */
type Representation = {
  uniqueId: string;
  visible: boolean;
  interfaceOption?: { visible?: boolean };
  show: () => void;
  hide: () => void;
};

function representationsOf(mol: moorhen.Molecule): Representation[] {
  return ((mol as unknown as { representations?: Representation[] })
    .representations ?? []);
}

/** Hide every representation, remembering which of them were drawn. */
export function hideMoleculeRepresentations(
  mol: moorhen.Molecule,
  memory: MoleculeVisibilityMemory,
): void {
  const drawn = new Set<string>();
  for (const rep of representationsOf(mol)) {
    if (rep.visible) drawn.add(rep.uniqueId);
    rep.hide();
  }
  if (mol.molNo != null) memory.set(mol.molNo, drawn);
}

/**
 * Put back exactly what `hideMoleculeRepresentations` took away.
 *
 * With no record of the hiding — a molecule hidden by some other control, or
 * by a Moorhen release whose own card did it — fall back to Moorhen's flag,
 * reading the `undefined` it leaves on our representations as "was drawn".
 * That is best-effort by construction: it is the state this module exists
 * because nobody recorded.
 */
export function showMoleculeRepresentations(
  mol: moorhen.Molecule,
  memory: MoleculeVisibilityMemory,
): void {
  const drawn = mol.molNo != null ? memory.get(mol.molNo) : undefined;
  for (const rep of representationsOf(mol)) {
    const wasDrawn = drawn
      ? drawn.has(rep.uniqueId)
      : rep.interfaceOption?.visible !== false;
    if (wasDrawn) rep.show();
  }
  if (mol.molNo != null) memory.delete(mol.molNo);
}
