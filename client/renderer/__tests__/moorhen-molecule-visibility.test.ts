/**
 * The eye icon in the Moorhen side panels: hide, then show, must put back
 * what was there.
 *
 * The bug this guards against was invisible in a type check and silent at
 * runtime. Moorhen's own MoleculeCard restores with
 * `rep.interfaceOption.visible`, which only Moorhen's own builder and its
 * representation chips ever set. Everything this app draws goes through the
 * deprecated `addRepresentation`, which leaves that flag `undefined`, so the
 * restore branch never ran: a molecule toggled off stayed off, with no error
 * anywhere. The two cases below are the two halves of getting it right —
 * scene-drawn representations come back, and the loader-default ribbon the
 * scene resolver hid stays hidden.
 */
import { describe, expect, it } from "vitest";

import {
  hideMoleculeRepresentations,
  showMoleculeRepresentations,
  type MoleculeVisibilityMemory,
} from "../lib/moorhen-molecule-visibility";

/** A representation as this module sees it, with show/hide that record. */
function makeRep(
  uniqueId: string,
  opts: { visible?: boolean; interfaceVisible?: boolean } = {},
) {
  return {
    uniqueId,
    visible: opts.visible ?? true,
    // `undefined` is what `molecule.addRepresentation` leaves here, and is
    // therefore the case that matters.
    interfaceOption: { visible: opts.interfaceVisible },
    show() {
      this.visible = true;
    },
    hide() {
      this.visible = false;
    },
  };
}

function makeMolecule(representations: ReturnType<typeof makeRep>[], molNo = 1) {
  return { molNo, representations } as never;
}

describe("molecule visibility toggle", () => {
  it("restores representations whose interfaceOption.visible was never set", () => {
    const sticks = makeRep("sticks");
    const mol = makeMolecule([sticks]);
    const memory: MoleculeVisibilityMemory = new Map();

    hideMoleculeRepresentations(mol, memory);
    expect(sticks.visible).toBe(false);

    showMoleculeRepresentations(mol, memory);
    expect(sticks.visible).toBe(true);
  });

  it("leaves representations that were already hidden hidden", () => {
    // What the scene resolver does: the loader's default ribbon is hidden so
    // the applied scene alone owns the look. Toggling the molecule off and on
    // must not resurrect it.
    const ghostRibbon = makeRep("loader-default-ribbon", { visible: false });
    const sceneSticks = makeRep("scene-sticks");
    const mol = makeMolecule([ghostRibbon, sceneSticks]);
    const memory: MoleculeVisibilityMemory = new Map();

    hideMoleculeRepresentations(mol, memory);
    showMoleculeRepresentations(mol, memory);

    expect(ghostRibbon.visible).toBe(false);
    expect(sceneSticks.visible).toBe(true);
  });

  it("keeps each molecule's record separate", () => {
    const a = makeRep("a");
    const bHidden = makeRep("b-hidden", { visible: false });
    const bShown = makeRep("b-shown");
    const molA = makeMolecule([a], 1);
    const molB = makeMolecule([bHidden, bShown], 2);
    const memory: MoleculeVisibilityMemory = new Map();

    hideMoleculeRepresentations(molA, memory);
    hideMoleculeRepresentations(molB, memory);
    showMoleculeRepresentations(molB, memory);

    expect(a.visible).toBe(false); // still hidden — nobody asked for it back
    expect(bHidden.visible).toBe(false);
    expect(bShown.visible).toBe(true);

    showMoleculeRepresentations(molA, memory);
    expect(a.visible).toBe(true);
  });

  it("survives a repeated toggle", () => {
    const sticks = makeRep("sticks");
    const mol = makeMolecule([sticks]);
    const memory: MoleculeVisibilityMemory = new Map();

    for (let i = 0; i < 3; i++) {
      hideMoleculeRepresentations(mol, memory);
      expect(sticks.visible).toBe(false);
      showMoleculeRepresentations(mol, memory);
      expect(sticks.visible).toBe(true);
    }
    // The record is spent once it has been used, so a stale entry cannot
    // outlive the molecule it described.
    expect(memory.size).toBe(0);
  });

  it("falls back to Moorhen's own flag when it did not do the hiding", () => {
    // A molecule hidden by some other control: no record to go on. Moorhen's
    // flag is honoured where it is set, and its `undefined` is read as "was
    // drawn" — the same reading that makes the normal path work.
    const explicitlyHidden = makeRep("chip-hidden", {
      visible: false,
      interfaceVisible: false,
    });
    const unset = makeRep("ours", { visible: false });
    const mol = makeMolecule([explicitlyHidden, unset]);

    showMoleculeRepresentations(mol, new Map());

    expect(explicitlyHidden.visible).toBe(false);
    expect(unset.visible).toBe(true);
  });
});
