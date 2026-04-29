/**
 * Client-side scaffold (substructure) membership matcher.
 *
 * The aggregation page already has SMILES per compound and an RDKit
 * context. We compute scaffold matches locally rather than asking the
 * server which compounds match each scaffold — sub-100ms for ~50
 * scaffolds × ~500 compounds, no per-compound round-trip.
 *
 * The `RDKitModule` type exposed by the core provider is a minimal
 * subset (just `get_mol(smiles)` for SVG rendering). MinimalLib's full
 * surface includes `get_qmol(smarts)` and `mol.get_substruct_match(qmol)`
 * — we widen via a local interface here rather than mutating the shared
 * type, so the broader app keeps its conservative typing.
 */

export interface ScaffoldDefinition {
  name: string;
  aliases: string[];
  smarts: string;
  source?: string;
}

export interface CompoundForMatching {
  formatted_id: string;
  smiles: string | null | undefined;
}

interface SubstructMatchResult {
  atoms?: number[];
  bonds?: number[];
}

interface RDKitMolWithMatch {
  get_substruct_match: (qmol: unknown) => string;
  /**
   * Canonical SMARTS string derived from this mol. Available on mols
   * parsed via `get_mol(smiles)` — RDKit assigns aromaticity during
   * the SMILES parse, so the resulting SMARTS uses lowercase aromatic
   * atoms even when the input was Kekulé. We rely on this to recover
   * a usable query mol when the chemist pastes a Kekulé-style SMILES
   * into the "+ SMARTS" input (e.g. `C12=C([N]C=C2)N=CN=C1` for
   * pyrrolopyrimidine — parses as SMARTS but the literal-character
   * query then matches zero aromatic compounds).
   */
  get_smarts: () => string;
  delete: () => void;
}

interface RDKitModuleWithSubstruct {
  get_mol: (smiles: string) => RDKitMolWithMatch | null;
  get_qmol: (smarts: string) => (RDKitMolWithMatch & { delete: () => void }) | null;
}

/**
 * Build a query mol from a chemist-supplied pattern. Tries SMILES-first
 * (which forces aromatic perception) and converts to canonical SMARTS;
 * falls back to raw SMARTS parsing for genuine SMARTS-only inputs that
 * SMILES parsing rejects (atom queries like `[#6X3]=[NX2]`, recursive
 * SMARTS, etc.).
 *
 * Without this dual strategy, Kekulé-style SMILES pasted as a SMARTS
 * (the natural mistake when a chemist copies from ChemDraw) produce a
 * literal-character query that matches nothing in an aromatic-perceived
 * compound database — silently zero, no error.
 */
function buildQueryMol(
  mod: RDKitModuleWithSubstruct,
  pattern: string,
): RDKitMolWithMatch | null {
  const tempMol = mod.get_mol(pattern);
  if (tempMol) {
    try {
      const aromaticSmarts = tempMol.get_smarts();
      if (aromaticSmarts) {
        const qmol = mod.get_qmol(aromaticSmarts);
        if (qmol) return qmol;
      }
    } catch {
      // get_smarts can throw on edge-case mols; fall through to the
      // raw-SMARTS path.
    } finally {
      tempMol.delete();
    }
  }
  return mod.get_qmol(pattern);
}

/**
 * For each scaffold, build the set of compound formatted_ids that
 * carry it. Scaffolds whose SMARTS fails to parse are skipped (a
 * curated catalog shouldn't have these, but `ScaffoldExtension` rows
 * are user-supplied — defence in depth).
 */
export function buildScaffoldMatches(
  rdkitModule: unknown,
  scaffolds: ReadonlyArray<ScaffoldDefinition>,
  compounds: ReadonlyArray<CompoundForMatching>,
): Map<string, Set<string>> {
  const out = new Map<string, Set<string>>();
  if (!rdkitModule) return out;
  const mod = rdkitModule as RDKitModuleWithSubstruct;
  if (typeof mod.get_qmol !== 'function') return out;

  // Cache parsed query mols across all scaffolds for the call.
  // buildQueryMol handles both real SMARTS and Kekulé SMILES (via
  // SMILES → canonical aromatic SMARTS round-trip).
  const queryMols: Array<{ name: string; qmol: RDKitMolWithMatch } | null> = scaffolds.map((s) => {
    const qmol = buildQueryMol(mod, s.smarts);
    if (!qmol) return null;
    return { name: s.name, qmol };
  });

  try {
    for (const c of compounds) {
      if (!c.smiles) continue;
      const mol = mod.get_mol(c.smiles);
      if (!mol) continue;
      try {
        for (const entry of queryMols) {
          if (!entry) continue;
          const raw = mol.get_substruct_match(entry.qmol);
          // get_substruct_match returns "{}" (JSON) when there's no
          // match — atoms array is absent / empty.
          let parsed: SubstructMatchResult;
          try {
            parsed = JSON.parse(raw || '{}');
          } catch {
            continue;
          }
          if (parsed.atoms && parsed.atoms.length > 0) {
            let bucket = out.get(entry.name);
            if (!bucket) {
              bucket = new Set();
              out.set(entry.name, bucket);
            }
            bucket.add(c.formatted_id);
          }
        }
      } finally {
        mol.delete();
      }
    }
  } finally {
    for (const entry of queryMols) {
      if (entry) entry.qmol.delete();
    }
  }
  return out;
}
