/**
 * Which code, into which molecule: the two questions "Add ligand here" has
 * to answer from CCP4i2's data rather than by guessing. The dictionary text
 * here follows real acedrg / refmac output (see demo_data/baz2b/*.cif).
 */
import { describe, expect, it, vi } from "vitest";
import { candidateLigandCodes, moleculeForFile, placeLigand } from "../lib/ligand-codes";

const compList = (...codes: string[]) =>
  "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.three_letter_code\n_chem_comp.name\n_chem_comp.group\n" +
  codes.map((c) => `${c} ${c} 'Unknown' non-polymer\n`).join("");

const compBlock = (code: string) =>
  `data_comp_${code}\n#\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n${code} C1\n${code} N1\n#\n`;

const dict = (...codes: string[]) =>
  "# acedrg 290\n" + compList(...codes) + "#\n" + codes.map(compBlock).join("");

describe("candidateLigandCodes", () => {
  it("reads the fragment out of an acedrg dictionary", () => {
    expect(candidateLigandCodes([dict("DRG")])).toEqual(["DRG"]);
  });

  it("ignores the standard monomers a refmac LIBOUT carries alongside the ligand", () => {
    expect(candidateLigandCodes([dict("ALA", "HOH", "GOL", "LIG")])).toEqual(["LIG"]);
  });

  it("offers both fragments of a two-fragment dictionary, sorted", () => {
    expect(candidateLigandCodes([dict("F02", "F01")])).toEqual(["F01", "F02"]);
  });

  it("merges several dictionaries without repeating a code", () => {
    expect(candidateLigandCodes([dict("DRG"), dict("DRG", "ALA")])).toEqual(["DRG"]);
  });
});

const fakeMolecule = (molNo: number, uniqueId?: string) => ({
  molNo,
  uniqueId,
  addLigandOfType: vi.fn(async () => ({})),
  redraw: vi.fn(async () => {}),
});

describe("moleculeForFile", () => {
  it("finds the molecule by the file it was loaded from, not by position", () => {
    const decoy = fakeMolecule(9, "/api/proxy/ccp4i2/files/77/download/");
    const target = fakeMolecule(7, "/api/proxy/ccp4i2/files/42/download/");
    // The decoy is last, i.e. the "active" molecule Moorhen would default to.
    expect(moleculeForFile([target, decoy], 42)).toBe(target);
    expect(moleculeForFile([target, decoy], 43)).toBeUndefined();
    expect(moleculeForFile([target, decoy], null)).toBeUndefined();
  });

  it("does not match a file id that is a suffix of another", () => {
    const other = fakeMolecule(1, "/api/proxy/ccp4i2/files/142/download/");
    expect(moleculeForFile([other], 42)).toBeUndefined();
  });
});

describe("placeLigand", () => {
  it("adds the ligand to the tracked molecule, looking the dictionary up on that same molecule", async () => {
    const decoy = fakeMolecule(9, "/api/proxy/ccp4i2/files/77/download/");
    const target = fakeMolecule(7, "/api/proxy/ccp4i2/files/42/download/");
    const placed = await placeLigand([target, decoy], 42, "DRG");
    expect(placed).toBe(target);
    expect(target.addLigandOfType).toHaveBeenCalledWith("DRG", 7);
    expect(target.redraw).toHaveBeenCalledTimes(1);
    expect(decoy.addLigandOfType).not.toHaveBeenCalled();
  });

  it("refuses when the tracked molecule is gone rather than picking another", async () => {
    const decoy = fakeMolecule(9, "/api/proxy/ccp4i2/files/77/download/");
    await expect(placeLigand([decoy], 42, "DRG")).rejects.toThrow(/no longer loaded/);
    expect(decoy.addLigandOfType).not.toHaveBeenCalled();
  });
});
