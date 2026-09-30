/**
 * The monomer editor turns clicks into a description of the modified monomer.
 *
 * Every click hands the caller a whole new MonomerEdits; the editor stores
 * nothing itself. GLU is the CCP4 library's (heavy atoms, charges, hydrogen
 * counts, Kekule bonds).
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";

const mocks = vi.hoisted(() => ({ rdkit: { module: null as any } }));

vi.mock("../providers/rdkit-provider", () => ({
  useRDKit: () => ({ rdkitModule: mocks.rdkit.module, isLoading: false, error: null }),
}));

import { MonomerEditor } from "../components/monomer/monomer-editor";
import { NO_EDITS, type MonomerEdits } from "../lib/monomer-edits";

const atom = (name: string, element: string, charge = 0, hydrogens = 0) => ({ name, element, charge, hydrogens });
const bond = (atom1: string, atom2: string, type = "single") => ({ atom1, atom2, type });

const GLU = {
  atom_details: [
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
const MONOMER = { atoms: GLU.atom_details.map((a) => a.name), ...GLU };

// A deterministic stand-in for MinimalLib's layout (see monomer-picker.test).
function fakeRDKit() {
  return {
    get_mol: (molblock: string) => {
      const count = parseInt(molblock.split("\n")[3].slice(0, 3), 10);
      return {
        set_new_coords: () => {},
        get_molblock: () =>
          ["", "", "", molblock.split("\n")[3],
            ...Array.from({ length: count }, (_, i) =>
              `${String(i).padStart(10)}${String((i * 7) % 5).padStart(10)}    0.0000 C   0  0`)].join("\n"),
        delete: vi.fn(),
      };
    },
  };
}

beforeEach(() => {
  mocks.rdkit.module = fakeRDKit();
});

function setup(edits: MonomerEdits = NO_EDITS, linkAtom = "CD") {
  const onEdits = vi.fn();
  const onPickLink = vi.fn();
  const view = render(
    <MonomerEditor
      monomer={MONOMER}
      code="GLU"
      emptyMessage="none"
      linkAtom={linkAtom}
      linkOrder="SINGLE"
      edits={edits}
      onPickLink={onPickLink}
      onEdits={onEdits}
    />
  );
  const atomEl = (name: string) => view.container.querySelector(`[data-atom="${name}"]`)!;
  const bondEl = (label: string) => view.container.querySelector(`[data-bond="${label}"]`)!;
  return { ...view, onEdits, onPickLink, atomEl, bondEl };
}

const mode = (name: string) => fireEvent.click(screen.getByRole("button", { name }));

describe("MonomerEditor", () => {
  it("picks the linking atom by default", () => {
    const { atomEl, onPickLink } = setup();
    fireEvent.click(atomEl("OE1"));
    expect(onPickLink).toHaveBeenCalledWith("OE1");
  });

  it("deletes an atom, and restores it on a second click", () => {
    const { atomEl, onEdits, rerender } = setup();
    mode("Delete");
    fireEvent.click(atomEl("OE2"));
    expect(onEdits).toHaveBeenLastCalledWith({ deletes: ["OE2"], bondOrders: [], charges: [] });
    const deleted = { deletes: ["OE2"], bondOrders: [], charges: [] };
    rerender(
      <MonomerEditor monomer={MONOMER} code="GLU" emptyMessage="none" linkAtom="CD"
        linkOrder="SINGLE" edits={deleted} onPickLink={vi.fn()} onEdits={onEdits} />
    );
    fireEvent.click(atomEl("OE2"));
    expect(onEdits).toHaveBeenLastCalledWith(NO_EDITS);
  });

  it("will not delete the linking atom, and says why", () => {
    const { atomEl, onEdits } = setup();
    mode("Delete");
    fireEvent.click(atomEl("CD"));
    expect(onEdits).not.toHaveBeenCalled();
    expect(screen.getByText("CD is the linking atom, so it cannot be deleted")).toBeTruthy();
  });

  it("sets a bond order from the bond's own row", () => {
    const { bondEl, onEdits } = setup();
    mode("Bond order");
    fireEvent.click(bondEl("CD-OE1"));
    expect(screen.getByText("dictionary: double")).toBeTruthy();
    fireEvent.click(screen.getByRole("button", { name: "single" }));
    expect(onEdits).toHaveBeenLastCalledWith({
      deletes: [], bondOrders: [{ atom1: "CD", atom2: "OE1", order: "SINGLE" }], charges: [],
    });
  });

  it("steps a formal charge from the dictionary's", () => {
    const { atomEl, onEdits } = setup();
    mode("Charge");
    fireEvent.click(atomEl("OE2"));
    expect(screen.getByTestId("charge-value").textContent).toBe("-1");
    fireEvent.click(screen.getByRole("button", { name: "Raise the charge" }));
    expect(onEdits).toHaveBeenLastCalledWith({
      deletes: [], bondOrders: [], charges: [{ atom: "OE2", charge: 0 }],
    });
  });

  it("lists the edits, each removable", () => {
    const edits = {
      deletes: ["OE2"], bondOrders: [{ atom1: "CD", atom2: "OE1", order: "DOUBLE" as const }],
      charges: [{ atom: "N", charge: 0 }],
    };
    const { onEdits } = setup(edits);
    expect(screen.getByText("delete OE2")).toBeTruthy();
    expect(screen.getByText("CD–OE1 double")).toBeTruthy();
    expect(screen.getByText("N 0")).toBeTruthy();
    const chip = screen.getByText("delete OE2").closest(".MuiChip-root")!;
    fireEvent.click(chip.querySelector(".MuiChip-deleteIcon")!);
    expect(onEdits).toHaveBeenLastCalledWith({ ...edits, deletes: [] });
  });

  it("draws the edits over the dictionary", () => {
    const { atomEl, bondEl } = setup({
      deletes: ["OE2"], bondOrders: [{ atom1: "CD", atom2: "OE1", order: "SINGLE" }], charges: [],
    });
    expect(atomEl("OE2").getAttribute("data-deleted")).toBe("true");
    expect(bondEl("CD-OE2").getAttribute("data-deleted")).toBe("true");
    expect(bondEl("CD-OE1").getAttribute("data-order")).toBe("single");
  });

  it("advises when a deletion leaves an atom short, as AceDRG will refuse it", () => {
    setup({ deletes: ["OXT"], bondOrders: [], charges: [] }, "");
    expect(screen.getByText("C would have valence 3; C usually has 4")).toBeTruthy();
  });

  it("says nothing about the isopeptide AceDRG accepts", () => {
    setup({ deletes: ["OE2"], bondOrders: [], charges: [] });
    expect(screen.queryByRole("alert")).toBeNull();
  });
});
