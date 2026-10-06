/**
 * ProvideAsuContents fills its table from three kinds of file. The digests
 * differ in shape; these helpers turn each into table rows and merge them
 * into what is already there without duplicating a sequence.
 */
import { describe, it, expect } from "vitest";
import {
  chainsFromDigest,
  classifyChain,
  deduplicateChains,
  mergeAsuEntries,
  polymerTypeFromMolecule,
} from "../components/task/task-elements/mmcif-sequence-parser";

describe("polymerTypeFromMolecule", () => {
  it("resolves CSequence's NUCLEIC to DNA or RNA by alphabet", () => {
    expect(polymerTypeFromMolecule("NUCLEIC", "ACGU")).toBe("RNA");
    expect(polymerTypeFromMolecule("NUCLEIC", "ACGT")).toBe("DNA");
    expect(polymerTypeFromMolecule("PROTEIN", "MKV")).toBe("PROTEIN");
  });
});

describe("classifyChain", () => {
  const composition = {
    peptides: ["A"],
    nucleics: ["B", "C"],
    chainDetails: [{ id: "D", type: "nucleic" }],
  };
  it("does not call a nucleic chain a protein", () => {
    expect(classifyChain("A", composition, "MKV")).toBe("PROTEIN");
    expect(classifyChain("B", composition, "ACGU")).toBe("RNA");
    expect(classifyChain("C", composition, "ACGT")).toBe("DNA");
    expect(classifyChain("D", composition, "ACGT")).toBe("DNA");
    expect(classifyChain("Z", composition, "")).toBe("OTHER");
  });
});

describe("chainsFromDigest", () => {
  it("reads a sequence-file digest as one chain", () => {
    const chains = chainsFromDigest(
      { name: "CDK2_HUMAN", moleculeType: "PROTEIN", sequence: "MEN FQK" },
      "cdk2.fasta",
    );
    expect(chains).toEqual([
      { chainId: "CDK2_HUMAN", sequence: "MENFQK", polymerType: "PROTEIN", length: 6, description: "cdk2.fasta" },
    ]);
  });

  it("reads a coordinate-file digest as its polymer chains, in composition order", () => {
    const chains = chainsFromDigest(
      {
        composition: { peptides: ["A", "C"], nucleics: ["B"] },
        sequences: { A: "MENFQK", B: "ACGU", C: "MENFQK", W: "x" },
      },
      "1cbs",
    );
    expect(chains.map((c) => c.chainId)).toEqual(["A", "C", "B"]);
    expect(chains[2].polymerType).toBe("RNA");
    expect(chains[0].description).toBe("1cbs chain A");
  });

  it("gives nothing for a digest without sequences", () => {
    expect(chainsFromDigest({ cell: {} }, "x.mtz")).toEqual([]);
    expect(chainsFromDigest(undefined, "x")).toEqual([]);
  });
});

describe("mergeAsuEntries", () => {
  const cdk2 = { name: "CDK2", sequence: "MENFQK", polymerType: "PROTEIN", description: "", nCopies: 1 };
  const cyclin = { name: "CCNA", sequence: "VPDYHE", polymerType: "PROTEIN", description: "", nCopies: 1 };

  it("appends new sequences and folds a repeated one into copies", () => {
    const incoming = deduplicateChains(
      chainsFromDigest(
        { composition: { peptides: ["A", "B", "C"] }, sequences: { A: "MENFQK", B: "VPDYHE", C: "MENFQK" } },
        "1fin",
      ),
    );
    const merged = mergeAsuEntries([cdk2], incoming);
    expect(merged).toHaveLength(2);
    expect(merged[0]).toMatchObject({ name: "CDK2", nCopies: 3 }); // 1 + chains A and C
    expect(merged[1]).toMatchObject({ sequence: "VPDYHE", nCopies: 1 });
  });

  it("leaves the existing rows untouched", () => {
    const existing = [cdk2];
    mergeAsuEntries(existing, [cyclin, { ...cdk2 }]);
    expect(existing[0].nCopies).toBe(1);
  });

  it("matches on polymer type as well as residues", () => {
    const dna = { ...cdk2, sequence: "ACGT", polymerType: "DNA" };
    const rna = { ...cdk2, sequence: "ACGT", polymerType: "RNA" };
    expect(mergeAsuEntries([dna], [rna])).toHaveLength(2);
  });
});
