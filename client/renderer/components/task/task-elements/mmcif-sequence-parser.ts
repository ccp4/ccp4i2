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

// ---------------------------------------------------------------------------
// From a file digest to ASU rows. Two digest shapes carry sequences:
//   a sequence file:   { name, moleculeType: "PROTEIN" | "NUCLEIC", sequence }
//   a coordinate file: { composition: { peptides, nucleics, chainDetails },
//                        sequences: { [chainId]: oneLetterSequence } }
// (an ASU file's digest is already a list of rows and needs none of this).
// ---------------------------------------------------------------------------

export type PolymerType = ChainSequenceInfo["polymerType"];

/** DNA or RNA, told apart by the alphabet: RNA has U, DNA has T. */
function nucleicType(sequence: string): "DNA" | "RNA" {
  return /U/i.test(sequence) ? "RNA" : "DNA";
}

/**
 * CSequence says PROTEIN or NUCLEIC; a CAsuContentSeq row says PROTEIN, DNA or
 * RNA. Resolve NUCLEIC from the sequence rather than guessing DNA.
 */
export function polymerTypeFromMolecule(moleculeType: string, sequence: string): PolymerType {
  const upper = (moleculeType || "").toUpperCase();
  if (upper === "PROTEIN" || upper === "DNA" || upper === "RNA") return upper;
  if (upper === "NUCLEIC") return nucleicType(sequence);
  return "OTHER";
}

/** Classify one chain of a coordinate-file digest from its composition lists. */
export function classifyChain(chainId: string, composition: any, sequence = ""): PolymerType {
  if (composition?.peptides?.includes(chainId)) return "PROTEIN";
  const detail = composition?.chainDetails?.find((d: any) => d.id === chainId);
  if (detail?.type === "protein") return "PROTEIN";
  if (composition?.nucleics?.includes(chainId) || detail?.type === "nucleic") {
    return nucleicType(sequence);
  }
  return "OTHER";
}

/**
 * Every polymer chain a digest describes, in file order. `label` names the
 * source for the row descriptions ("1cbs", "my_model.pdb").
 */
export function chainsFromDigest(digest: any, label: string): ChainSequenceInfo[] {
  if (!digest) return [];
  if (digest.moleculeType && digest.sequence) {
    const sequence = String(digest.sequence).replace(/\s/g, "");
    return [{
      chainId: digest.name || label,
      sequence,
      polymerType: polymerTypeFromMolecule(digest.moleculeType, sequence),
      length: sequence.length,
      description: digest.description || label,
    }];
  }
  if (digest.composition && digest.sequences) {
    const ids: string[] = [
      ...(digest.composition.peptides || []),
      ...(digest.composition.nucleics || []),
    ];
    return ids
      .filter((chainId) => digest.sequences[chainId])
      .map((chainId) => {
        const sequence = String(digest.sequences[chainId]);
        return {
          chainId,
          sequence,
          polymerType: classifyChain(chainId, digest.composition, sequence),
          length: sequence.length,
          description: `${label} chain ${chainId}`,
        };
      });
  }
  return [];
}

/**
 * Append `incoming` to `existing`: a sequence already in the table (same
 * polymer type and residues) adds its copies to that row instead of
 * repeating it, so loading a second source never produces duplicates.
 */
export function mergeAsuEntries<T extends AsuSequenceEntry>(
  existing: T[],
  incoming: AsuSequenceEntry[],
): (T | AsuSequenceEntry)[] {
  const key = (e: AsuSequenceEntry) =>
    `${String(e.polymerType).toUpperCase()}:${String(e.sequence).replace(/\s/g, "")}`;
  const merged: (T | AsuSequenceEntry)[] = existing.map((e) => ({ ...e }));
  const byKey = new Map(merged.map((e) => [key(e), e]));
  for (const entry of incoming) {
    const hit = byKey.get(key(entry));
    if (hit) {
      hit.nCopies = (Number(hit.nCopies) || 0) + entry.nCopies;
    } else {
      merged.push(entry);
      byKey.set(key(entry), entry);
    }
  }
  return merged;
}
