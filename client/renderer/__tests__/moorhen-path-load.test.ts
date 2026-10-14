import { describe, expect, it } from "vitest";
import { classifyCif, isLoadableFile, planPathLoad } from "../lib/moorhen-path-load";

describe("planPathLoad", () => {
  it("splits a PanDDA dataset folder into models, maps and CIFs to read", () => {
    const d = "CCP4_JOBS/job_3/pandda2_out/processed_datasets/xtal-0001";
    const plan = planPathLoad([
      `${d}/xtal-0001-z_map.native.ccp4`,
      `${d}/xtal-0001-pandda-input.pdb`,
      `${d}/xtal-0001-event_1_1-BDC_0.27_map.native.ccp4`,
      `${d}/events.yaml`,
      `${d}/LIG.cif`,
    ]);
    expect(plan.coordinates).toEqual([`${d}/xtal-0001-pandda-input.pdb`]);
    expect(plan.maps).toEqual([
      `${d}/xtal-0001-event_1_1-BDC_0.27_map.native.ccp4`,
      `${d}/xtal-0001-z_map.native.ccp4`,
    ]);
    expect(plan.cifs).toEqual([`${d}/LIG.cif`]);
  });

  it("takes paths without a Moorhen extension for folders", () => {
    expect(isLoadableFile("a/b/model.pdb")).toBe(true);
    expect(isLoadableFile("a/b/xtal-0001")).toBe(false);
  });
});

describe("classifyCif", () => {
  it("knows coordinates by their atom_site loop", () => {
    expect(classifyCif("data_x\nloop_\n_atom_site.group_PDB\n_atom_site.id\n")).toBe("coordinates");
  });

  it("knows a monomer dictionary", () => {
    expect(classifyCif("data_comp_list\nloop_\n_chem_comp.id\nLIG\ndata_comp_LIG\nloop_\n_chem_comp_atom.comp_id\n")).toBe("dictionary");
  });

  it("says nothing of other CIFs", () => {
    expect(classifyCif("data_x\n_cell.length_a 10\n")).toBeNull();
  });
});
