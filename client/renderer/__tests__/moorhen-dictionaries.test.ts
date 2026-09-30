import { describe, expect, it } from "vitest";
import {
  buildableCompIds,
  extractDictCompIds,
  loadWithDictionaries,
  privateDictionaryNumber,
  provenanceOf,
} from "../lib/moorhen-dictionaries";

const CIF = "data_comp_list\nloop_\n_chem_comp.id\nLIG\ndata_comp_LIG\n#\ndata_comp_DRG\n";

class FakeMolecule {
  molNo: number | null = null;
  events: string[] = [];
  async addDict(text: string) {
    this.events.push(`addDict:${text.slice(0, 4)}`);
  }
  async loadMissingMonomers() {
    this.events.push("missing-monomers");
  }
  async loadToCoot() {
    this.molNo = 3;
    this.events.push("coords");
    await this.loadMissingMonomers(); // what Moorhen does inside the load
  }
}

describe("extractDictCompIds", () => {
  it("lists the components and skips the list header", () => {
    expect(extractDictCompIds(CIF)).toEqual(["LIG", "DRG"]);
  });
});

// Shaped like AcedrgLink's two dictionaries: the linked pair (a whole
// monomer with atoms) and the link itself (modifications and a link only).
const LINKED_PAIR = `data_comp_list
loop_
_chem_comp.id
LIG
data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
LIG N
`;
const LINK_ONLY = `data_program_info
data_mod_list
data_link_list
data_mod_LYSm1
loop_
_chem_mod_atom.mod_id
LYSm1
data_link_LYS-GLU
`;

describe("buildableCompIds", () => {
  it("lists the components that have atoms", () => {
    expect(buildableCompIds(LINKED_PAIR)).toEqual(["LIG"]);
  });

  it("finds nothing to build in a link dictionary", () => {
    expect(buildableCompIds(LINK_ONLY)).toEqual([]);
  });

  it("skips a component block with no atoms", () => {
    expect(buildableCompIds(CIF)).toEqual([]);
  });
});

describe("privateDictionaryNumber", () => {
  it("never repeats and stays clear of Coot's special numbers and molecule numbers", () => {
    const a = privateDictionaryNumber();
    const b = privateDictionaryNumber();
    expect(b).not.toBe(a);
    expect(Math.min(a, b)).toBeGreaterThanOrEqual(1_000_000);
  });
});

describe("loadWithDictionaries", () => {
  it("attaches the job's dictionaries before Moorhen looks for missing monomers", async () => {
    const mol = new FakeMolecule();
    await loadWithDictionaries(mol, [{ text: "AAAA" }, { text: "BBBB" }], () => mol.loadToCoot());
    expect(mol.events).toEqual(["coords", "addDict:AAAA", "addDict:BBBB", "missing-monomers"]);
  });

  it("restores the molecule's own method afterwards", async () => {
    const mol = new FakeMolecule();
    await loadWithDictionaries(mol, [{ text: "AAAA" }], () => mol.loadToCoot());
    expect(Object.prototype.hasOwnProperty.call(mol, "loadMissingMonomers")).toBe(false);
    mol.events = [];
    await mol.loadMissingMonomers();
    expect(mol.events).toEqual(["missing-monomers"]);
  });

  it("leaves a load with no dictionaries exactly as Moorhen does it", async () => {
    const mol = new FakeMolecule();
    await loadWithDictionaries(mol, [], () => mol.loadToCoot());
    expect(mol.events).toEqual(["coords", "missing-monomers"]);
  });

  it("restores the method even when the load throws, and attaches nothing", async () => {
    const mol = new FakeMolecule();
    await expect(
      loadWithDictionaries(mol, [{ text: "AAAA" }], async () => {
        throw new Error("bad coordinates");
      }),
    ).rejects.toThrow("bad coordinates");
    expect(mol.events).toEqual([]);
    expect(Object.prototype.hasOwnProperty.call(mol, "loadMissingMonomers")).toBe(false);
  });
});

describe("provenanceOf", () => {
  it("maps every component to the file that defines it, skipping unsourced text", () => {
    const map = provenanceOf([{ text: CIF, fileId: 7, projectId: "p" }, { text: "data_comp_XYZ\n" }]);
    expect(Array.from(map.entries())).toEqual([
      ["LIG", { fileId: 7, projectId: "p" }],
      ["DRG", { fileId: 7, projectId: "p" }],
    ]);
  });
});
