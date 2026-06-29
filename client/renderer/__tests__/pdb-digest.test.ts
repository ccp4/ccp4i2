/**
 * PDB-id extraction + PDBe digest formatting for the authoring prompt.
 * The NL→`pdb:` choice happens in the LLM, so the app digests PDB ids it can
 * see in the request text and embeds chains + ligand CIDs.
 */
import { describe, it, expect, vi } from "vitest";

import {
  extractPdbIds,
  formatPdbContents,
  fetchPdbContents,
} from "../lib/moorhen-scene-prompt";

describe("extractPdbIds", () => {
  it("finds the id in the guinea-pig request, no false positives", () => {
    expect(
      extractPdbIds(
        "Using PDB file 1OGU, could we centre and slab on centre of mass of chains A and B",
      ),
    ).toEqual(["1OGU"]);
  });

  it("uppercases + dedupes classic ids and matches extended ids", () => {
    expect(extractPdbIds("compare 1abc with 1ABC")).toEqual(["1ABC"]);
    expect(extractPdbIds("entry pdb_00001abc")).toEqual(["pdb_00001abc"]);
  });

  it("does not match plain words or chain selectors", () => {
    expect(extractPdbIds("show the dimer of chains A and B as ribbon")).toEqual([]);
  });
});

describe("formatPdbContents", () => {
  it("lists polymer chains (named) and ligand CIDs", () => {
    const block = formatPdbContents(
      "1OGU",
      [
        { molecule_type: "polypeptide(L)", molecule_name: ["Kinase"], in_chains: ["A", "B"] },
        { molecule_type: "water", in_chains: ["A"] },
      ],
      [{ chem_comp_id: "GOL", chain_id: "A", author_residue_number: 501 }],
    );
    expect(block).toContain("PDB 1OGU (from PDBe)");
    expect(block).toContain("chain A: Kinase");
    expect(block).toContain("chain B: Kinase");
    expect(block).not.toContain("water"); // non-polymer excluded from chains
    expect(block).toContain("GOL at //A/501");
  });
});

describe("fetchPdbContents", () => {
  const ok = (body: unknown) =>
    ({ ok: true, json: async () => body }) as unknown as Response;

  it("fetches molecules + ligand_monomers and formats", async () => {
    const fetchFn = vi.fn(async (url: string) =>
      url.includes("ligand_monomers")
        ? ok({ "1ogu": [{ chem_comp_id: "GOL", chain_id: "A", author_residue_number: 501 }] })
        : ok({ "1ogu": [{ molecule_type: "polypeptide(L)", in_chains: ["A"] }] }),
    ) as unknown as typeof fetch;
    const block = await fetchPdbContents("1OGU", fetchFn);
    expect(block).toContain("PDB 1OGU");
    expect(block).toContain("GOL at //A/501");
  });

  it("returns null when the entry 404s (over-matched id self-filters)", async () => {
    const fetchFn = vi.fn(
      async () => ({ ok: false, status: 404 }) as unknown as Response,
    ) as unknown as typeof fetch;
    expect(await fetchPdbContents("2XYZ", fetchFn)).toBeNull();
  });
});
