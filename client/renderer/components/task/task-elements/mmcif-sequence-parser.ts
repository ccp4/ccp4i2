/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
/**
 * Utilities for parsing multi-chain FASTA responses from PDB services
 * and extracting chain sequence information.
 */

export interface ChainSequenceInfo {
  chainId: string;
  sequence: string;
  polymerType: "PROTEIN" | "DNA" | "RNA" | "OTHER";
  length: number;
  description: string;
}

/**
 * Classify polymer type from FASTA header description.
 * PDBe/RCSB headers typically contain "mol:protein", "mol:na" etc.
 */
function classifyFromHeader(header: string): "PROTEIN" | "DNA" | "RNA" | "OTHER" {
  const lower = header.toLowerCase();
  if (lower.includes("mol:protein") || lower.includes("polypeptide")) return "PROTEIN";
  if (lower.includes("mol:na")) {
    if (lower.includes("dna") || lower.includes("deoxy")) return "DNA";
    return "RNA";
  }
  // Heuristic: if sequence contains U, likely RNA; if mostly ACGT, DNA
  return "PROTEIN"; // Default assumption for PDB entries
}

/**
 * Parse a multi-entry FASTA string into individual chain sequences.
 *
 * Handles FASTA from PDBe (`/pdbe/entry/pdb/{id}/fasta`) and
 * RCSB (`/fasta/entry/{ID}`) which return all chains in one response.
 *
 * Header formats:
 *   PDBe: >pdb|1cbs|A Chain A, ...
 *   RCSB: >1CBS_1|Chain A|...
 */
/**
 * ASU content entry with copy count, produced by deduplicating chains.
 */
export interface AsuSequenceEntry {
  name: string;
  sequence: string;
  polymerType: string;
  description: string;
  nCopies: number;
}

/**
 * Deduplicate chains by identical sequence + polymerType.
 * Chains with the same sequence are collapsed into a single entry
 * with nCopies reflecting the count and a combined name.
 *
 * e.g. CDK2 chains A,C → { name: "Chain_A_C", nCopies: 2 }
 */
export function deduplicateChains(chains: ChainSequenceInfo[]): AsuSequenceEntry[] {
  const groups = new Map<string, { chains: ChainSequenceInfo[]; polymerType: string }>();
  for (const chain of chains) {
    // Key on sequence content + polymer type
    const key = `${chain.polymerType}:${chain.sequence}`;
    const existing = groups.get(key);
    if (existing) {
      existing.chains.push(chain);
    } else {
      groups.set(key, { chains: [chain], polymerType: chain.polymerType });
    }
  }

  const entries: AsuSequenceEntry[] = [];
  for (const group of groups.values()) {
    const chainIds = group.chains.map((c) => c.chainId);
    const first = group.chains[0];
    entries.push({
      name: `Chain_${chainIds.join("_")}`,
      sequence: first.sequence,
      polymerType: first.polymerType,
      description: first.description || chainIds.map((id) => `Chain ${id}`).join(", "),
      nCopies: group.chains.length,
    });
  }
  return entries;
}

export function parseMultiChainFasta(fastaText: string): ChainSequenceInfo[] {
  const entries: ChainSequenceInfo[] = [];
  const lines = fastaText.split("\n");

  let currentHeader = "";
  let currentSeq = "";

  const flush = () => {
    if (!currentHeader || !currentSeq) return;

    // Extract chain ID from header
    let chainId = "?";
    const description = currentHeader.substring(1).trim(); // Remove leading >

    // PDBe format: >pdb|1cbs|A ...
    const pdbeMatch = description.match(/^pdb\|\w+\|(\w+)/);
    // RCSB format: >1CBS_1|Chains A, B|... or >1CBS_A|...
    const rcsbMatch = description.match(/^\w+_(\w+)\|/);
    // Generic: look for "Chain X" or "Chains X, Y"
    const chainMatch = description.match(/[Cc]hains?\s+([A-Za-z0-9,\s]+)/);

    if (pdbeMatch) {
      chainId = pdbeMatch[1];
    } else if (chainMatch) {
      // May be "Chains A, B" — take first chain letter
      chainId = chainMatch[1].split(/[,\s]+/)[0].trim();
    } else if (rcsbMatch) {
      chainId = rcsbMatch[1];
    }

    entries.push({
      chainId,
      sequence: currentSeq,
      polymerType: classifyFromHeader(description),
      length: currentSeq.length,
      description,
    });
  };

  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed.startsWith(">")) {
      flush();
      currentHeader = trimmed;
      currentSeq = "";
    } else if (trimmed) {
      currentSeq += trimmed;
    }
  }
  flush(); // Don't forget last entry

  return entries;
}
