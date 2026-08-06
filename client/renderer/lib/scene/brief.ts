/**
 * Compact LLM authoring brief, generated from the committed JSON Schema
 * contract (the CI-gated artefact) — so it can never drift from the format,
 * and is ~10x lighter than embedding the full prose grammar.
 *
 * The brief conveys SHAPE + enums + terse per-field docs. Worked examples and
 * prose conventions (CID syntax, ||-join, by-domain) are NOT here — they stay
 * curated in the prompt, because the schema can't express them.
 *
 * See MOORHEN_SCENES_SCHEMA_V1_DESIGN.md §7a.
 */
import schema from "./moorhen-scene.ccp4i2.v1.json";

interface Node {
  type?: string | string[];
  enum?: unknown[];
  const?: unknown;
  anyOf?: Node[];
  oneOf?: Node[];
  properties?: Record<string, Node>;
  required?: string[];
  items?: Node | Node[];
  prefixItems?: Node[];
  description?: string;
  format?: string;
}

const INLINE_MAX = 72; // inline an object/array on one line if it fits

function enumStr(values: unknown[]): string {
  return values.map((v) => JSON.stringify(v)).join("|");
}

function typeStr(node: Node, indent: string): string {
  if (node.const !== undefined) return JSON.stringify(node.const);
  if (node.enum) return enumStr(node.enum);
  const variants = node.anyOf ?? node.oneOf;
  if (variants) return variants.map((v) => typeStr(v, indent)).join(" | ");
  if (node.prefixItems) {
    return `[${node.prefixItems.map((n) => typeStr(n, indent)).join(", ")}]`;
  }
  const t = Array.isArray(node.type) ? node.type[0] : node.type;
  switch (t) {
    case "string":
    case "number":
    case "integer":
    case "boolean":
      return t === "integer" ? "number" : t;
    case "array": {
      const it = Array.isArray(node.items) ? node.items[0] : node.items;
      if (!it) return "any[]";
      const inner = typeStr(it, indent);
      return /\s\|\s|\|/.test(inner) ? `(${inner})[]` : `${inner}[]`;
    }
    case "object":
      return objStr(node, indent);
    default:
      return "any";
  }
}

/** Render an object — inline if short, else one field per line. */
function objStr(node: Node, indent: string): string {
  const props = node.properties ?? {};
  const required = new Set(node.required ?? []);
  const keys = Object.keys(props);
  if (keys.length === 0) return "{}";

  const inline = `{ ${keys
    .map((k) => `${k}${required.has(k) ? "" : "?"}: ${typeStr(props[k], indent)}`)
    .join(", ")} }`;
  if (inline.length <= INLINE_MAX && !inline.includes("\n")) return inline;

  const inner = indent + "  ";
  const lines = keys.map((k) => {
    const child = props[k];
    const opt = required.has(k) ? "" : "?";
    const doc = child.description ? `  # ${child.description}` : "";
    return `${inner}${k}${opt}: ${typeStr(child, inner)}${doc}`;
  });
  return `{\n${lines.join("\n")}\n${indent}}`;
}

/**
 * Human-readable Markdown docs for the scene format. The Grammar section is
 * generated from the contract (always in sync); the surrounding prose +
 * examples are curated here. Regenerated + drift-checked in the test suite, so
 * it can never go stale the way the hand-maintained doc did.
 */
export function buildSceneMarkdown(): string {
  return `<!-- GENERATED from the Zod schema via buildSceneMarkdown(); do not edit by hand.
     The Grammar section below is authoritative. Edit lib/scene/core.ts (the Zod
     source) and regenerate (UPDATE_SCHEMA=1 vitest run scene-schema). -->
# Moorhen Scene format

A Moorhen **scene** is a portable YAML description of how to render one or more
structures: which files to load, which domains to recognise, what
representations and colours to apply, and where the camera sits. Scenes validate
against the published contract (\`moorhen-scene.core.v1.json\`, or
\`moorhen-scene.ccp4i2.v1.json\` for project references) and can be re-applied
across different structures of the same protein.

## Grammar

Fields marked \`?\` are optional. \`A | B\` is a choice; \`T[]\` is a list of \`T\`.

\`\`\`yaml
${buildSceneBrief()}
\`\`\`

## File references

Each \`files[]\` entry sets **exactly one** source:

| Ref | Meaning | Portable? |
|-----|---------|-----------|
| \`pdb\` | PDB id, fetched via the PDBe proxy | yes (core) |
| \`url\` | absolute URL | yes (core) |
| \`bundle\` | asset inside a \`.scene.zip\` | yes (core) |
| \`cifText\` | inline dictionary CIF (\`kind: dictionary\` only) | yes (core) |
| \`relativeUrl\` | origin-relative URL (\`/api/…\`) | within one deployment |
| \`fileId\` (+\`projectId\`) | ccp4i2 project file | within one deployment |
| \`job\`+\`param\` (+\`projectId\`) | ccp4i2 job output | within one deployment |

The resolver matches a ref against an already-loaded molecule by its loader URL
(or by \`fileId\` extracted from it); a \`pdb\`/\`url\`/\`bundle\`/\`job+param\` ref is
fetched if not present. A molecule's identity is its loader URL, never a bare
path — so \`relativeUrl\` matches loaded molecules, it does not read local files.

## Selections (Coot CIDs)

Selections everywhere (representation \`selection\`, \`view.centre\`/\`slab\`,
domains) are Coot CIDs:

- \`//A\` — the whole of chain A
- \`//A/703-740\` — residues 703–740 of chain A
- \`//*/LIG\` — every residue named LIG
- \`//A/750/CA\` — one atom
- join several with \`||\`: \`//A||//B\`

A representation draws **its own** \`selection\` (the whole molecule if omitted);
colour does **not** limit what is drawn — scope the selection to limit it.

## Colours

A \`colour\` is one of:

1. a hex string — \`"#4b8bbe"\` (or \`#rrggbbaa\`);
2. a named scheme — \`by-domain\`, \`b-factor\`, \`b-factor-norm\`, \`af2-plddt\`,
   \`secondary-structure\`, \`jones-rainbow\`, \`mol-symm\`;
3. a per-selection list — \`[{ selection: "//A", colour: "#a08766" }, …]\`;
4. a raw Coot rule (escape hatch).

\`by-domain\` colours by the top-level \`domains:\` block: define domains
(\`name\` + \`selection\` + \`color\`) and set \`colour: by-domain\` on the reps.

Colour cascades over two levels (matching Moorhen's molecule-level colour):
\`element.colour\` is the molecule-wide default inherited by every representation
of that file, and a representation's own \`colour\` overrides it for that rep.

## Geometry and hints

- \`geometry\` (per representation) is **honoured** — physical dimensions in
  Ångström a renderer must reproduce (bond/ball radii, ribbon widths, …).
- \`hints\` (scene-level lighting + effects) are **advisory** — a renderer may
  ignore them and still produce a correct image. \`hints.lighting.direction\` is
  a direction vector (+x right, +y up, +z toward the viewer); effects are
  scene-authoritative (an \`effects\` block fully determines effect state).

## Examples

Minimal portable scene:

\`\`\`yaml
scene: minimal
version: 1
files:
  - { name: prot, pdb: 1ABC }
elements:
  - file: prot
    representations:
      - { style: CRs, selection: "//A", colour: secondary-structure }
\`\`\`

Domains, honoured geometry, advisory hints (validates under the core schema):

\`\`\`yaml
scene: gamma-demo
version: 1
files:
  - { name: gamma, pdb: 1B9K }
domains:
  - { name: head,   selection: "//A/703-740",  color: "#4b8bbe" }
  - { name: tail,   selection: "//A/900-1000", color: "#e74c3c" }  # clamps to 938
elements:
  - file: gamma
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }
      - { style: CBs, selection: "//A/750", geometry: { bondRadius: 0.18 } }
hints:
  lighting: { direction: [1, 1, 1], shininess: 24 }
  effects:  { ssao: { enabled: true }, edgeDetect: { enabled: true } }
resolver:
  onMissingResidues: clamp-and-log
\`\`\`
`;
}

/** Build the compact grammar brief from the committed contract. */
export function buildSceneBrief(): string {
  const root = schema as unknown as Node;
  const props = root.properties ?? {};
  const required = new Set(root.required ?? []);
  const lines: string[] = [
    "Moorhen scene — YAML grammar (fields, ? = optional):",
    "",
  ];
  for (const key of Object.keys(props)) {
    const child = props[key];
    const opt = required.has(key) ? "" : "?";
    const doc = child.description ? `  # ${child.description}` : "";
    lines.push(`${key}${opt}: ${typeStr(child, "")}${doc}`);
  }
  return lines.join("\n");
}
