/**
 * The picker renders dictionary atom names and reports them on a click.
 *
 * That is the whole contract: the task parameters hold strings like "NZ" and
 * "C4A", so the depiction has to be keyed by those names rather than by the
 * indices a toolkit-drawn picture would give.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";

const mocks = vi.hoisted(() => ({
  rdkit: {
    isLoading: false,
    module: null as any,
  },
}));

vi.mock("../providers/rdkit-provider", () => ({
  useRDKit: () => ({
    rdkitModule: mocks.rdkit.module,
    isLoading: mocks.rdkit.isLoading,
    error: null,
  }),
}));

import { MonomerPicker } from "../components/monomer/monomer-picker";
import { buildMolblock } from "../lib/monomer-molblock";

const ATOMS = [
  { name: "N", element: "N", charge: 1 },
  { name: "CA", element: "C", charge: 0 },
  { name: "C", element: "C", charge: 0 },
  { name: "O", element: "O", charge: 0 },
];
const BONDS = [
  { atom1: "N", atom2: "CA", type: "single" },
  { atom1: "CA", atom2: "C", type: "single" },
  { atom1: "C", atom2: "O", type: "double" },
];

/**
 * A stand-in for MinimalLib that lays atoms out on a diagonal. The real
 * library is exercised by monomer-molblock-rdkit.test.ts; here the point is
 * the rendering and the callbacks, so a deterministic layout is better than
 * a real one.
 */
function fakeRDKit({ reject = false } = {}) {
  return {
    get_mol: (molblock: string) => {
      if (reject) return null;
      const atomCount = parseInt(molblock.split("\n")[3].slice(0, 3), 10);
      return {
        set_new_coords: () => {},
        get_molblock: () => {
          const lines = ["", "", "", molblock.split("\n")[3]];
          for (let i = 0; i < atomCount; i += 1) {
            lines.push(
              `${String(i).padStart(10)}${String(i * 2).padStart(10)}` +
                "    0.0000 C   0  0"
            );
          }
          return lines.join("\n");
        },
        delete: vi.fn(),
      };
    },
  };
}

beforeEach(() => {
  mocks.rdkit.module = fakeRDKit();
  mocks.rdkit.isLoading = false;
});

const svgAtoms = (container: HTMLElement) =>
  Array.from(container.querySelectorAll("[data-atom]")).map(
    (el) => el.getAttribute("data-atom")!
  );

describe("MonomerPicker", () => {
  it("labels every atom with its dictionary name", () => {
    const { container } = render(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(svgAtoms(container).sort()).toEqual(["C", "CA", "N", "O"]);
  });

  it("draws every bond, keyed independently of atom order", () => {
    const { container } = render(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(container.querySelectorAll("[data-bond]")).toHaveLength(BONDS.length);
  });

  it("reports the atom name, not an index, when one is clicked", () => {
    const onPickAtom = vi.fn();
    const { container } = render(
      <MonomerPicker atomDetails={ATOMS} bonds={BONDS} mode="atom" onPickAtom={onPickAtom} />
    );
    fireEvent.click(container.querySelector('[data-atom="CA"]')!);
    expect(onPickAtom).toHaveBeenCalledWith("CA");
  });

  it("reports both atom names when a bond is clicked", () => {
    const onPickBond = vi.fn();
    const { container } = render(
      <MonomerPicker atomDetails={ATOMS} bonds={BONDS} mode="bond" onPickBond={onPickBond} />
    );
    fireEvent.click(container.querySelector('[data-bond="C-O"]')!);
    expect(onPickBond).toHaveBeenCalledWith("C", "O");
  });

  it("ignores clicks in the mode that is not active", () => {
    const onPickAtom = vi.fn();
    const onPickBond = vi.fn();
    const { container } = render(
      <MonomerPicker
        atomDetails={ATOMS}
        bonds={BONDS}
        mode="bond"
        onPickAtom={onPickAtom}
        onPickBond={onPickBond}
      />
    );
    fireEvent.click(container.querySelector('[data-atom="CA"]')!);
    expect(onPickAtom).not.toHaveBeenCalled();
  });

  it("is reachable from the keyboard when it is interactive", () => {
    const onPickAtom = vi.fn();
    render(
      <MonomerPicker atomDetails={ATOMS} bonds={BONDS} mode="atom" onPickAtom={onPickAtom} />
    );
    const button = screen.getByRole("button", { name: "Atom CA" });
    fireEvent.keyDown(button, { key: "Enter" });
    expect(onPickAtom).toHaveBeenCalledWith("CA");
  });

  it("marks the selection so the user can see what is chosen", () => {
    const { container } = render(
      <MonomerPicker atomDetails={ATOMS} bonds={BONDS} mode="atom" selectedAtom="CA" onPickAtom={vi.fn()} />
    );
    expect(container.querySelector('[data-atom="CA"]')).toHaveAttribute("aria-pressed", "true");
    expect(container.querySelector('[data-atom="N"]')).toHaveAttribute("aria-pressed", "false");
  });

  it("treats a selected bond the same whichever way round it is given", () => {
    const { container } = render(
      <MonomerPicker
        atomDetails={ATOMS}
        bonds={BONDS}
        mode="bond"
        selectedBond={{ atom1: "O", atom2: "C" }}
        onPickBond={vi.fn()}
      />
    );
    expect(container.querySelector('[data-bond="C-O"]')).toHaveAttribute("aria-pressed", "true");
  });

  it("offers no buttons at all when it is not a picker", () => {
    render(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} mode="none" />);
    expect(screen.queryAllByRole("button")).toHaveLength(0);
  });

  it("says why it is empty rather than rendering a blank box", () => {
    // Three different reasons, three different messages: a blank panel
    // leaves the user with nothing to act on.
    const { rerender } = render(<MonomerPicker atomDetails={[]} bonds={[]} />);
    expect(screen.getByText("No monomer selected")).toBeTruthy();

    mocks.rdkit.module = null;
    mocks.rdkit.isLoading = true;
    rerender(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(screen.getByText(/Preparing the chemistry toolkit/)).toBeTruthy();

    mocks.rdkit.isLoading = false;
    rerender(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(screen.getByText(/could not be loaded/)).toBeTruthy();

    mocks.rdkit.module = fakeRDKit({ reject: true });
    rerender(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(screen.getByText(/could not be drawn/)).toBeTruthy();
  });

  it("releases the WASM molecule it built", () => {
    // MinimalLib objects are not garbage collected; a picker that re-lays out
    // on every parameter edit would leak one each time.
    const deletes: any[] = [];
    mocks.rdkit.module = {
      get_mol: () => {
        const mol = {
          set_new_coords: () => {},
          get_molblock: () => buildMolblock(ATOMS, BONDS),
          delete: vi.fn(),
        };
        deletes.push(mol.delete);
        return mol;
      },
    };
    render(<MonomerPicker atomDetails={ATOMS} bonds={BONDS} />);
    expect(deletes.length).toBeGreaterThan(0);
    expect(deletes.every((d) => d.mock.calls.length === 1)).toBe(true);
  });
});
