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
): string {
  const jobByPk = new Map<number, ManifestJob>(jobs.map((j) => [j.id, j]));
  // group renderable referenceable files by their producing job. Include both
  // job-output-dir (1) and import-dir (2) files: imported coordinates carry a
  // job + param too and are valid job+param targets.
  const byJob = new Map<number, ManifestFile[]>();
  for (const f of files) {
    if (!f.job_param_name) continue;
    if (!f.type || !RENDERABLE_TYPES.has(f.type)) continue;
    const arr = byJob.get(f.job) ?? [];
    arr.push(f);
    byJob.set(f.job, arr);
  }
  if (byJob.size === 0) return "No referenceable files in this project yet.";

  const lines: string[] = [];
  lines.push("Reference files by { job: <number>, param: <PARAM>, projectId }.");
  lines.push('A quoted name (e.g. — "CDK2") is a human label to help you choose the');
  lines.push("right file; it is NOT a reference key — always reference by job + param.");
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
      // The annotation is a human label (e.g. "CDK2-Cyclin A") to help pick the
      // right file — it is NOT how the file is referenced (use job + param).
      const ann = f.annotation ? `  — "${f.annotation}"` : "";
      lines.push(`    - param ${f.job_param_name} -> ${f.type}${mask}${ann}`);
    }
  }
  return lines.join("\n");
}

// ── Interpretation guidance ─────────────────────────────────────────────────

// Static guidance teaching the model to turn everyday phrasing into concrete
// selections against the contents it was given. The high-value, reliably
// meetable cases are back-references ("the dimer" = the chains just named) and
// "the ligand" = the single ligand present; genuinely ambiguous cases (which
// dimer in a tetramer) are handled by "pick the most likely reading and proceed".
const INTERPRETATION_GUIDANCE = [
  "=== INTERPRETING THE REQUEST ===",
  "Turn everyday phrasing into concrete selections using the CONTENTS and PROJECT",
  "above. Never invent chains, ligands, or residues that are not listed.",
  "- Back-references point to what the request already named: \"the dimer\", \"the",
  "  complex\", \"it\", \"both\", \"that\" refer to the chains/entities mentioned earlier",
  "  in the SAME request. After \"chains A and B ...\", \"the dimer\" means chains A and B.",
  "- A representation draws its own `selection` (the WHOLE molecule if you omit it),",
  "  so to show only certain chains/residues you MUST set that selection. Colour does",
  "  NOT limit what is drawn: \"chains A and B as ribbon\" needs a CRs representation",
  "  with selection //A||//B — an unscoped one draws every chain (C and D too), even",
  "  if the colours only cover A and B.",
  "- When a representation is requested \"of\"/\"for\" an entity, scope THAT",
  "  representation's selection to the same entity. \"a surface of the dimer\" after",
  "  \"chains A and B\" => a MolecularSurface whose selection covers chains A and B.",
  "- \"the protein\" / \"the model\" / \"everything\" => all polymer chains.",
  "- \"monomer\" => one polymer chain; \"dimer\"/\"trimer\"/... => that many polymer",
  "  chains (prefer the chains the request names; otherwise the polymer chains",
  "  present, in order).",
  "- A ligand named or given by 3-letter code (\"ATP\", \"the ligand\", \"the inhibitor\",",
  "  \"the drug\") => the matching ligand CID from the contents. If exactly one ligand",
  "  is present, \"the ligand\" means it.",
  "- \"active site\" / \"binding site\" / \"around the ligand\" => the ligand and its",
  "  surroundings (a neighbourhood selection if the grammar supports one, else the",
  "  ligand's chain/residue).",
  "- Colour cues: \"by domain\" => the domains block; \"by chain\" / \"rainbow\" /",
  "  \"spectrum\" => the matching colour scheme in the grammar.",
  "- To cover several chains or residue ranges in ONE selection, join them with",
  "  \"||\" (e.g. //A||//B). This works in representation selections AND in view",
  "  directives (centre, slab) — use the same form throughout.",
  "- Camera vs depth are SEPARATE directives. `view.centre` moves the camera to a",
  "  selection; `view.slab` only sets the clip DEPTH (it does not move the camera).",
  "  So \"centre on chains A and B AND slab to contain them\" needs BOTH:",
  "  view.centre.selection = //A||//B  and  view.slab.selection = //A||//B. Emitting",
  "  slab alone will not centre the view.",
  "For minor ambiguity, choose the most likely reading and proceed. Ask a concise",
  "clarifying question ONLY when the ambiguity would materially change the scene and",
  "no reasonable default exists (e.g. which two chains form \"the dimer\" in a",
  "tetramer) — ask once, then return the YAML code block once answered.",
].join("\n");

// ── Whole prompt ────────────────────────────────────────────────────────────

export function buildAuthoringPrompt(opts: {
  project?: { name?: string; id?: string };
  contents: string;
  manifest: string;
  request: string;
}): string {
  const projName = opts.project?.name;
  const projId = opts.project?.id;
  // Pin the project identity hard and up front: the single most important
  // constraint, because a chat model carrying context from an earlier project
  // will otherwise reference the wrong one. Every job/param ref must carry these.
  const projectLines: string[] = ["=== CURRENT PROJECT (use ONLY this one) ==="];
  if (projName) projectLines.push(`projectName: ${projName}`);
  if (projId) projectLines.push(`projectId: ${projId}`);
  if (!projName && !projId) {
    projectLines.push("(project identity unavailable — ask the user before referencing job/param files)");
  }
  projectLines.push(
    "Every { job, param } reference MUST include the projectId above (or projectName).",
    "Ignore any project mentioned earlier in this conversation; use ONLY the project above.",
  );

  return [
    "You are drafting a Moorhen \"scene\": a YAML document describing a molecular",
    "view. Follow the grammar below EXACTLY. Return the scene as a SINGLE fenced YAML",
    "code block (```yaml ... ```) and nothing else outside it — YAML is",
    "whitespace-sensitive, and a code block lets the chat UI show a copy button that",
    "preserves the exact indentation (selecting the text by hand does not). If — and",
    "only if — the request is materially ambiguous and no reasonable default exists,",
    "you may FIRST ask one concise clarifying question in plain text and wait for the",
    "reply before producing the code block.",
    "Reference project files with { job, param, projectId } using the manifest; use",
    "pdb:/url: only for structures not in the project; never use path:.",
    "Use the exact ligand CIDs from the contents summary for ligand selections.",
    "",
    projectLines.join("\n"),
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
    INTERPRETATION_GUIDANCE,
    "",
    "=== REQUEST ===",
    "Produce a scene YAML that satisfies the following request. Use the project,",
    "files, and CIDs above. For minor ambiguity choose sensible defaults; only for a",
    "material ambiguity with no reasonable default, ask one concise question first.",
    "When you produce the scene, return it as a SINGLE ```yaml fenced code block and",
    "nothing else outside it. The request:",
    "",
    opts.request.trim() || "(no request given — produce a clear default overview of the loaded structure)",
    "",
  ].join("\n");
}
