// Build the "Copy prompt" scaffold for authoring a Moorhen scene with an LLM.
//
// The prompt has four parts: instructions, the scene grammar (embedded verbatim
// from types/moorhen-scene.md — the same contract the validator enforces), a
// manifest of the project's referenceable files, and a ground-truth summary of
// the loaded structure's chains and ligands. The user appends their request and
// pastes the whole thing into any chatbot; the returned YAML goes back into the
// scene editor, which validates and applies it. No API key, provider-free.

import sceneGrammar from "../types/moorhen-scene.md";

// ── Contents summary (from the coordinate `digest` endpoint) ────────────────

export interface DigestLigand {
  chainId: string;
  name: string;
  seqNum: number;
  cid: string;
  nAtoms: number;
}

export interface FileDigest {
  sequences?: Record<string, string>;
  ligands?: DigestLigand[];
}

/** Format a coordinate file's digest (per-chain sequences + ligand CIDs) into a
 *  contents block. Gives the author chain ids and exact ligand CIDs so it never
 *  has to guess a selection. */
export function buildContentsBlock(digest: FileDigest, label: string): string {
  const lines: string[] = [];
  const seqs = digest.sequences ?? {};
  const chains = Object.keys(seqs);
  if (chains.length) {
    lines.push("Polymer chains (reference as e.g. //A or selection ranges):");
    for (const c of chains) {
      lines.push(`  - chain ${c}: ${seqs[c].length} residues`);
    }
  }
  const ligs = digest.ligands ?? [];
  if (ligs.length) {
    lines.push("Ligands / ions (reference by the exact CID shown):");
    for (const l of ligs) {
      lines.push(`  - ${l.name} at ${l.cid}  (${l.nAtoms} atoms)`);
    }
  }
  if (!lines.length) return `Loaded structure "${label}": no chains or ligands reported.`;
  return `Loaded structure "${label}":\n${lines.join("\n")}`;
}

// ── Project manifest (from jobs + files REST) ───────────────────────────────

export interface ManifestJob {
  id: number;
  number?: string;
  task_name?: string;
  title?: string;
}

export interface ManifestFile {
  job: number;
  job_param_name?: string;
  type?: string;
  sub_type?: number | null;
  annotation?: string;
  name?: string;
  directory: number;
}

// Renderable output types worth advertising to the author. Other outputs (logs,
// dictionaries, etc.) are omitted to keep the manifest focused.
const RENDERABLE_TYPES = new Set([
  "chemical/x-pdb",
  "chemical/x-mmcif",
  "chemical/x-cif",
  "application/CCP4-map",
  "application/CCP4-mtz-map",
  "application/CCP4-mtz",
]);

/** Format the project's job outputs into a manifest the author references via
 *  `{ job: <number>, param: <PARAM>, projectId }`. One group per job, listing
 *  its renderable output params. */
export function buildManifestBlock(
  jobs: ManifestJob[],
  files: ManifestFile[],
  projectId: string | undefined,
): string {
  const jobByPk = new Map<number, ManifestJob>(jobs.map((j) => [j.id, j]));
  // group renderable output files by their producing job
  const byJob = new Map<number, ManifestFile[]>();
  for (const f of files) {
    if (f.directory !== 1) continue; // job outputs only
    if (!f.job_param_name) continue;
    if (!f.type || !RENDERABLE_TYPES.has(f.type)) continue;
    const arr = byJob.get(f.job) ?? [];
    arr.push(f);
    byJob.set(f.job, arr);
  }
  if (byJob.size === 0) return "No referenceable job outputs in this project yet.";

  const lines: string[] = [];
  if (projectId) lines.push(`projectId for every job/param reference: ${projectId}`);
  lines.push("Referenceable files (use { job: <number>, param: <PARAM>, projectId }):");
  // stable order by job number (numeric-ish)
  const jobPks = [...byJob.keys()].sort((a, b) => {
    const na = parseFloat(jobByPk.get(a)?.number ?? "0");
    const nb = parseFloat(jobByPk.get(b)?.number ?? "0");
    return na - nb;
  });
  for (const pk of jobPks) {
    const job = jobByPk.get(pk);
    const header = `  job ${job?.number ?? "?"}${job?.task_name ? ` (${job.task_name})` : ""}:`;
    lines.push(header);
    for (const f of byJob.get(pk)!) {
      const mask = f.sub_type === 4 ? ", mask" : "";
      lines.push(`    - param ${f.job_param_name} -> ${f.type}${mask}`);
    }
  }
  return lines.join("\n");
}

// ── Whole prompt ────────────────────────────────────────────────────────────

export function buildAuthoringPrompt(opts: {
  contents: string;
  manifest: string;
}): string {
  return [
    "You are drafting a Moorhen \"scene\": a YAML document describing a molecular",
    "view. Follow the grammar below EXACTLY and output ONLY the YAML (no prose, no",
    "code fences). Reference project files with { job, param, projectId } using the",
    "manifest; use pdb:/url: only for structures not in the project; never use path:.",
    "Use the exact ligand CIDs from the contents summary for ligand selections.",
    "",
    "=== SCENE GRAMMAR ===",
    sceneGrammar,
    "",
    "=== REFERENCEABLE PROJECT FILES ===",
    opts.manifest,
    "",
    "=== LOADED STRUCTURE CONTENTS ===",
    opts.contents,
    "",
    "=== YOUR REQUEST ===",
    "Describe the view you want, then produce the scene YAML:",
    "",
  ].join("\n");
}
