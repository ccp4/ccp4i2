// Build the "Copy prompt" scaffold for authoring a Moorhen scene with an LLM.
//
// The prompt has four parts: instructions, the scene grammar (a compact brief
// GENERATED from the committed JSON Schema contract, so it can never drift from
// what the validator enforces), a manifest of the project's referenceable
// files, and a ground-truth summary of the loaded structure's chains and
// ligands. The user appends their request and pastes the whole thing into any
// chatbot; the returned YAML goes back into the scene editor, which validates
// and applies it. No API key, provider-free.

import { buildSceneBrief } from "./scene/brief";

// Conventions the schema can't express (CID syntax, by-domain, colour forms,
// ||-join) — the curated semantic supplement to the generated grammar.
const SCENE_CONVENTIONS = [
  "=== CONVENTIONS ===",
  "- Selections are Coot CIDs: //A (whole chain A), //A/703-740 (residue range),",
  "  //*/LIG (residues named LIG), //A/750/CA (one atom). Join several with || :",
  "  //A||//B. Same syntax in representation `selection` and in view.centre/slab.",
  "- A representation draws its own `selection` (the WHOLE molecule if omitted);",
  "  colour does NOT limit what is drawn — scope the selection to limit it.",
  "- colour is a hex \"#rrggbb\"; OR a named scheme (by-domain, b-factor, af2-plddt,",
  "  secondary-structure, jones-rainbow, mol-symm); OR a per-selection list",
  "  [{selection, colour}]. `by-domain` colours by the top-level `domains:` block:",
  "  define domains (name + selection + color) and set `colour: by-domain` on reps.",
  "- `element.colour` sets a molecule-wide default (any colour form above) inherited by",
  "  every representation of that file; a representation's own `colour` overrides it.",
  "- geometry dimensions are Ångström. `hints` (lighting/effects) are ADVISORY —",
  "  a viewer may ignore them; never rely on them for what must be visible.",
].join("\n");

const WORKED_EXAMPLE = [
  "=== EXAMPLE (shape only — use the PROJECT/CONTENTS below for real refs) ===",
  "```yaml",
  "scene: example",
  "version: 1",
  "files:",
  "  - { name: prot, pdb: 1ABC }",
  "domains:",
  "  - { name: nterm, selection: \"//A/1-100\",   color: \"#4b8bbe\" }",
  "  - { name: cterm, selection: \"//A/101-200\", color: \"#e74c3c\" }",
  "elements:",
  "  - file: prot",
  "    representations:",
  "      - { style: CRs, selection: \"//A\", colour: by-domain }",
  "view:",
  "  centre: { file: prot, selection: \"//A\" }",
  "```",
].join("\n");

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

// ── PDB references in the request (digest from PDBe) ────────────────────────

// The NL→`pdb:` decision happens inside the LLM, which has no access to our
// endpoints — so a `pdb:` ref can't be digested lazily. Instead the APP scans
// the request TEXT for PDB ids it can see ("…from pdb 1OGU…"), fetches their
// chains + ligand CIDs from PDBe (via the COEP proxy), and embeds them, exactly
// like a loaded-structure digest. Only fires when an id is literally present;
// a purely NL reference ("CDK2") degrades gracefully to no digest.

/** Extract candidate PDB ids from free text. Classic ids start with a digit
 *  (1-9) + 3 alphanumerics; extended ids are `pdb_########`. Over-matching is
 *  harmless — a non-existent id 404s at PDBe and is dropped. */
export function extractPdbIds(text: string): string[] {
  const ids = new Set<string>();
  for (const m of text.matchAll(/\b[1-9][a-zA-Z0-9]{3}\b/g)) ids.add(m[0].toUpperCase());
  for (const m of text.matchAll(/\bpdb_[0-9a-z]{8}\b/gi)) ids.add(m[0].toLowerCase());
  return [...ids];
}

interface PdbeMolecule {
  molecule_type?: string;
  molecule_name?: string[] | string;
  in_chains?: string[];
}
interface PdbeLigand {
  chem_comp_id?: string;
  chain_id?: string;
  author_residue_number?: number;
}

const MAX_LIGAND_LINES = 30;

/** Format PDBe `molecules` + `ligand_monomers` arrays into a contents block,
 *  matching buildContentsBlock so the LLM reads PDB and project refs alike. */
export function formatPdbContents(
  id: string,
  molecules: PdbeMolecule[],
  ligands: PdbeLigand[],
): string {
  const isPolymer = (t?: string) =>
    !!t && /polypeptide|polyribo|polydeoxyribo|polynucleotide|nucleotide|peptide|saccharide|carbohydrate/i.test(t);
  const lines: string[] = [];

  const chainLines: string[] = [];
  for (const m of molecules) {
    if (!isPolymer(m.molecule_type)) continue;
    const name = Array.isArray(m.molecule_name) ? m.molecule_name[0] : m.molecule_name;
    for (const c of m.in_chains ?? []) chainLines.push(`  - chain ${c}${name ? `: ${name}` : ""}`);
  }
  if (chainLines.length) {
    lines.push("Polymer chains (reference as e.g. //A or selection ranges):");
    lines.push(...chainLines);
  }

  const ligLines = ligands
    .filter((l) => l.chem_comp_id)
    .map((l) => `  - ${l.chem_comp_id} at //${l.chain_id}/${l.author_residue_number}`);
  if (ligLines.length) {
    lines.push("Ligands / ions (reference by the exact CID shown):");
    lines.push(...ligLines.slice(0, MAX_LIGAND_LINES));
    if (ligLines.length > MAX_LIGAND_LINES)
      lines.push(`  - … and ${ligLines.length - MAX_LIGAND_LINES} more`);
  }

  if (!lines.length) return `PDB ${id} (from PDBe): no chains or ligands reported.`;
  return `PDB ${id} (from PDBe):\n${lines.join("\n")}`;
}

/** Fetch + format a PDB id's contents from PDBe via the proxy. Returns null on
 *  404/error (so over-matched ids self-filter). `fetchFn`/`base` injectable for tests. */
export async function fetchPdbContents(
  id: string,
  fetchFn: typeof fetch = fetch,
  base = "/api/proxy/pdbe/api/pdb/entry",
): Promise<string | null> {
  const lid = id.toLowerCase();
  try {
    const [mr, lr] = await Promise.all([
      fetchFn(`${base}/molecules/${lid}`),
      fetchFn(`${base}/ligand_monomers/${lid}`),
    ]);
    if (!mr.ok) return null;
    const mj = (await mr.json()) as Record<string, PdbeMolecule[]>;
    const lj = lr.ok ? ((await lr.json()) as Record<string, PdbeLigand[]>) : {};
    const molecules = mj[lid] ?? [];
    if (!molecules.length) return null;
    return formatPdbContents(id, molecules, lj[lid] ?? []);
  } catch {
    return null;
  }
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

// ── System prompt (static, cacheable) ────────────────────────────────────────

/**
 * The static, query-general half of the authoring prompt: instructions +
 * grammar + conventions + worked example + interpretation guidance. It contains
 * NOTHING that varies per project or request, so it is byte-identical across
 * every call and every user.
 *
 * That invariance is the whole point. It is (a) the cacheable prefix for a chat
 * model's / Azure OpenAI's prefix cache — it must sit ahead of the first
 * variable byte or that variable's churn re-invalidates it; and (b) the system
 * message Materia's `nlp_scene` endpoint holds server-side. It is committed
 * verbatim as the drift-tested artifact `lib/scene/moorhen-scene.system-prompt.v1.md`
 * so a pinned submodule has a stable system-prompt source. Keep this function
 * free of per-call inputs.
 */
export function buildSceneSystemPrompt(): string {
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
    "pdb:/url: only for structures not in the project; never use relativeUrl:.",
    "Use the exact ligand CIDs from the contents summary for ligand selections.",
    "",
    "=== SCENE GRAMMAR ===",
    buildSceneBrief(),
    "",
    SCENE_CONVENTIONS,
    "",
    WORKED_EXAMPLE,
    "",
    INTERPRETATION_GUIDANCE,
  ].join("\n");
}

// ── Whole prompt ────────────────────────────────────────────────────────────

export function buildAuthoringPrompt(opts: {
  project?: { name?: string; id?: string };
  contents: string;
  manifest: string;
  request: string;
}): string {
  const projName = opts.project?.name;
  const projId = opts.project?.id;
  // Pin the project identity hard: a chat model carrying context from an earlier
  // project will otherwise reference the wrong one. Every job/param ref must
  // carry these. This is the FIRST variable element — everything above it (the
  // system prompt) is static and cacheable; everything from here down varies per
  // call, ordered most-stable-to-least (project → manifest → contents → request).
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
    buildSceneSystemPrompt(),
    "",
    projectLines.join("\n"),
    "",
    "=== REFERENCEABLE PROJECT FILES ===",
    opts.manifest,
    "",
    "=== LOADED STRUCTURE CONTENTS ===",
    opts.contents,
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
