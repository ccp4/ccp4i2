# Moorhen Scene format — v1

A **scene** is a portable, human-editable YAML description of how to
render one or more structures in Moorhen. It captures *intent*
(domains, colour rules, representations, superpositions, camera)
separately from any specific PDB file, so the same visual treatment
can be re-applied across different structures of the same protein.

This document is the authoring reference — written for a person, a
script, or a **language model** building a scene from outside Moorhen.
It is the hand-written companion to the two generated contracts in this
directory: `moorhen-scene.system-prompt.v1.md` (the grammar, generated
from the Zod schema) and `moorhen-scene.structured.v1.json` (the strict
Structured-Outputs profile). Those two state what is *well-formed*;
this one carries what is *correct in practice* — the pitfalls that
produce a scene which validates cleanly and still renders nothing.

Read it alongside them, not instead of them.

**Using this with any model.** Everything here is provider-neutral. The
files ship as `ccp4i2` package data, so an adopter can load them
directly:

```python
from importlib.resources import files
contracts = files("ccp4i2") / "scene_contracts"
grammar   = (contracts / "moorhen-scene.system-prompt.v1.md").read_text()
authoring = (contracts / "moorhen-scene-authoring.v1.md").read_text()
# Concatenate as the system message for whichever model you use;
# constrain the response with moorhen-scene.structured.v1.json where
# the provider supports JSON-Schema-constrained output.
```

The canonical implementation is the Zod schema in the CCP4i2 frontend
(`client/renderer/lib/scene/core.ts`), with TypeScript types in
`client/renderer/types/moorhen-scene.ts`.

## File extensions

| Suffix             | Purpose                                                              |
| ------------------ | -------------------------------------------------------------------- |
| `*.scene.yaml`     | The scene as a plain YAML file. Use when all file refs are portable. |
| `*.scene.zip`      | Bundle: `scene.yaml` at root + `assets/` directory of attached data. |
| `*.session.json`   | Moorhen-native cache, regenerable. Not for hand-authoring.           |

The bundle is the right shape when the scene references local files
(coords or dictionaries) that aren't otherwise reachable by URL or PDB
ID. The yaml inside uses `bundle: <relpath>` to point at zipped
assets. Both forms parse identically — the yaml is the source of truth.

## Top-level grammar

```yaml
scene: <string>                      # required: human-readable identifier
version: 1                           # required: schema version

authoredIn:                          # optional: provenance (never resolved)
  projectId: <uuid>                  # optional
  projectName: <name>                # optional
  createdAt: <iso-8601>              # optional
  createdBy: <author>                # optional
  ccp4i2Version: <string>            # optional

files: [ ... ]                       # optional: named coord + dict refs
superpose: [ ... ]                   # optional: alignments, applied in order
globalDictionaries: [ <name>, ... ]  # optional: dicts loaded on every molecule
domains: [ ... ]                     # optional: reusable named residue blocks
elements: [ ... ]                    # optional: per-file rendering instructions
view: { ... }                        # optional: camera, clip, fog, background

resolver:                            # optional: apply-time policy
  onMissingResidues: clamp-and-log   # | strict
```

Order of evaluation at apply-time:

1. **Fetch dictionaries** (load each globally, so coords parse).
2. **Fetch coordinates** (in declared order).
3. **Scope dictionaries** per element (re-associate to that molecule's molNo).
4. **Run superpositions** (mutate moving structures' transforms).
5. **Apply representations** per element.
6. **Set camera**.

## `files`

A list of named file references. The name (e.g. `protein`, `apo`,
`x0034`) is what other parts of the scene refer to; it's local to this
file and doesn't have to mean anything outside.

Each entry has:

- `name`: **required**, unique within the block.
- `kind`: optional, `"coordinates"` (default) or `"dictionary"`.
- Exactly one resolution method:

  | Field                              | Use when                                                                                                                 |
  | ---------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
  | `pdb: <id>`                        | Deposited structure. Fetched via PDBe. Most portable. Coordinates only.                                                  |
  | `url: <https://...>`               | Coord or dict published at a CORS-friendly URL.                                                                          |
  | `fileId: <int>` + `projectId: <uuid>` | Project-internal ref to a ccp4i2 file.                                                                                |
  | `job: <int>` + `param: <str>` + `projectId` | Project-internal ref to a job's output parameter.                                                              |
  | `bundle: <relpath>`                | Asset packaged inside a `.scene.zip` (e.g. `assets/coords/x0034.cif`).                                                   |
  | `cifText: \|<inline cif>`          | Dictionary text inlined directly. Only valid for `kind: dictionary`.                                                     |
  | `path: <abspath>`                  | Local-only. Won't resolve on another machine. Mostly used as a marker.                                                   |

### Examples

```yaml
files:
  # A deposited structure — most portable.
  - { name: ref, pdb: 1m17 }

  # Project-internal ccp4i2 file.
  - name: refined
    projectId: 3f8a-aaaa-bbbb-cccc-uuid
    fileId: 482

  # Multiple coord files via a single bundle.
  - { name: x0034, bundle: assets/coords/x0034.cif }
  - { name: x0092, bundle: assets/coords/x0092.cif }

  # A scoped dictionary, bundled alongside its coords.
  - { name: lig-A, kind: dictionary, bundle: assets/dict/lig-A.cif }

  # Inline dict (cifText) — verbose but always works on any machine.
  - name: minimal-dict
    kind: dictionary
    cifText: |
      data_comp_LIG
      _chem_comp.id LIG
      ...
```

## `domains`

Reusable named residue blocks, referenced from `colour: by-domain`
inside an element. Hoisted to the top level so a multi-structure scene
doesn't duplicate them per element.

```yaml
domains:
  - { name: CARD, selection: "//A/1-92",    color: "#4b8bbe" }
  - { name: NBD,  selection: "//A/104-260", color: "#2ecc71" }
  - { name: HD1,  selection: "//A/261-330", color: "#9b59b6" }
```

Fields:

- `name`: required.
- `selection`: required (see the deprecated alternative below). Any
  valid Coot CID — the same grammar as a representation `selection`:
  - `"//A/1-92"` — residues 1–92 of chain A.
  - `"//A"` — the whole of chain A.
  - `"//*/32-64"` — that range on every chain (useful for symmetric
    assemblies like the apoptosome heptamer).
  - `"//A/(ALA,GLY)"`, `"//A/55/CA[C]"` — residue names, atoms, and
    other things the old chain+range form could not express.
- `color`: required hex `#rrggbb` or `#rrggbbaa`.

When the CID has the plain `//chain/start-end` shape, the resolver
clamps the range to the residues actually present in the loaded
structure and logs what it trimmed (see `resolver.onMissingResidues`);
any richer CID is passed straight to Coot. A domain whose range lies
entirely outside the model is skipped with a log entry, not an error.

### Deprecated: `chain` + `range`

Older scenes wrote a domain as `chain:` plus `range:` instead of a
`selection:`. That form still validates and still resolves, but it is
deprecated — **write new scenes with `selection:`**. Set one form or the
other, never both; a domain with neither is rejected.

```yaml
domains:
  - { name: NBD, chain: A, range: 104-260, color: "#2ecc71" }  # deprecated
  - { name: NBD, selection: "//A/104-260", color: "#2ecc71" }  # preferred
```

In the deprecated form `chain` accepts `"A"`, `"*"`, or `["A","B"]`, and
an omitted `range` means the whole chain.

Note that `superpose` also has `chain` and `range` fields; those are
**not** deprecated — see [`superpose`](#superpose).

## `elements`

Per-file rendering instructions: which representations to draw, with
what colour, on which selection.

```yaml
elements:
  - file: protein                        # name from the files: block
    dictionaries: [lig-A, cofactor]      # optional: per-molecule dict scoping
    representations:
      - { style: CRs,      selection: "//A",       colour: by-domain }
      - { style: ligands,  selection: "//*/(LIG)", colour: "#2ecc71" }
      # Multi-residue highlight: || joins single-residue / range CIDs.
      - { style: CBs,      selection: "//A/115||//A/116||//A/121-122", colour: "#e74c3c" }
```

### `style`

A Moorhen representation style string. The values below are the
authoritative set, taken from Moorhen's
[`baby-gru/src/utils/enums.ts`](https://github.com/moorhen-coot/Moorhen/blob/main/baby-gru/src/utils/enums.ts)
(`representationLabelMapping`) and the dispatch switch in
[`MoorhenMoleculeRepresentation.ts`](https://github.com/moorhen-coot/Moorhen/blob/main/baby-gru/src/utils/MoorhenMoleculeRepresentation.ts).
**Case matters** — the resolver passes the string through verbatim, so
`MetaBalls` works and `metaballs` does not.

#### Everyday styles — what you'll author 95% of the time

| Style              | Label (in the UI)  | What it draws                                     |
| ------------------ | ------------------ | ------------------------------------------------- |
| `CRs`              | Ribbons            | Cartoon ribbons. The default for showing a fold.  |
| `CBs`              | Bonds              | All-atom sticks (carbon-bond style).              |
| `CAs`              | C-Alpha            | Cα-only trace.                                     |
| `MolecularSurface` | Mol Surface        | Smooth solvent-excluded surface. Shows pockets.   |
| `gaussian`         | Gaussian Surface   | Softer, blobbier surface than `MolecularSurface`. |
| `VdwSpheres`       | Spheres            | Van-der-Waals spheres (CPK).                       |
| `ligands`          | Ligands            | Sticks restricted to HET groups. **Always pair with an explicit `//*/(COMPID)` selection** (see pitfalls). |
| `MetaBalls`        | MetaBalls          | Smooth fused-blob isosurface around the selected atoms. Reads as a solid "plug" — excellent for a ligand sitting in a pocket where sticks look thin. |
| `DishyBases`       | Bases              | Cartoon nucleotide bases (DNA/RNA).               |
| `StickBases`       | —                  | Stick-style nucleotide bases (alternative to `DishyBases`). |
| `allHBonds`        | H-Bonds            | Hydrogen bonds as dashed lines.                   |
| `glycoBlocks`      | Glyco-Blocks       | SNFG-style block cartoons for glycans/sugars.     |

#### Choosing between the surfaces and the blobby styles

These four all produce "solid" geometry; pick by intent:

- **`MolecularSurface`** — the protein's outer surface; use it to reveal
  a *cleft or pocket* (e.g. show a surface on the receptor, sticks/
  metaballs on the ligand inside it).
- **`gaussian`** — a smoother, lower-detail surface; good for a clean
  poster silhouette where atomic bumpiness is distracting.
- **`MetaBalls`** — fuses selected atoms into rounded merged blobs. Best
  on a *small* selection (a ligand, a handful of residues). On a whole
  protein it's heavy and loses the fold. This is the style to reach for
  when a ligand drawn as `ligands` (thin sticks) doesn't read clearly.
- **`VdwSpheres`** — hard CPK spheres; more literal/atomic than
  `MetaBalls`, no fusing between atoms.

`MetaBalls` is a *Coot bond representation* (same family as `CBs`,
`VdwSpheres`, `CAs`, `ligands`): it is selected by CID and coloured by
the element's `colour:` rule exactly like those. So swapping
`style: ligands` → `style: MetaBalls` on the same selection/colour is a
safe, drop-in change. Its internal grid/radius parameters are not
exposed in the scene format — you get Moorhen's defaults.

#### Validation / analysis styles — valid, but situational

These render real geometry and are legal in a scene; you'd use them for
a validation or analysis view rather than a presentation figure:

| Style              | Label             | Draws                                          |
| ------------------ | ----------------- | ---------------------------------------------- |
| `CDs`              | Contact dots      | All-atom contact dots (clashes/contacts).      |
| `rama`             | Ramachandran Balls | Per-residue Ramachandran validation markers.  |
| `rotamer`          | Rotamer Dodec.    | Rotamer-quality dodecahedra.                   |
| `restraints`       | Restraints        | Geometry restraints.                           |
| `adaptativeBonds`  | Adaptive Bonds    | Bonds whose detail adapts to zoom (note the spelling — `adaptative`, from the source). |

#### Do NOT author these — internal / interaction states

The style union in the source also contains values that are driven by
Moorhen's UI at runtime (hovering, selecting, transforming) or need a
context the scene format can't supply. They are **not** meant to be set
from a scene and will render nothing useful (or nothing at all):
`hover`, `residueSelection`, `transformation`, `unitCell`,
`environment`, `ligand_environment`, `residue_environment`,
`contact_dots`, `chemical_features`, `ligand_validation`, `VdWSurface`,
`Calpha`. (Note `VdWSurface`/`Calpha` are the *M2T-internal* spellings;
for scenes use `MolecularSurface` and `CAs` instead.)

If you genuinely need one of these, confirm against the current Moorhen
source before authoring — don't infer behaviour from the name.

### `selection`

A Coot CID string (gemmi-parsed selection syntax). Format:
`/<model>/<chain>/<residue>/<atom>`. Wildcards: `*` matches any one
component; omitting trailing components defaults them to `*`.

The residue field accepts **either a single residue number, or a
single `start-end` range** — *not* a comma-separated list. Comma lists
ARE valid for residue *names* (`(ALA,GLY,SER)`) and *elements*
(`[N,O,S]`), but not for residue numbers.

| Pattern                    | Meaning                                            |
| -------------------------- | -------------------------------------------------- |
| `/*/*/*/*`                 | Every atom (the default if `selection` is omitted) |
| `//A`                      | Chain A, all residues                              |
| `//A/115`                  | Residue 115 of chain A                             |
| `//A/115-200`              | Residues 115–200 inclusive                         |
| `//*/(LIG)`                | Every chain, residue name LIG (parens are required)|
| `//A/(ALA,GLY)`            | Chain A, residue name ALA or GLY                   |
| `//A/115/CA[C]`            | The Cα atom of residue 115, element carbon         |

**Comma-separated residue numbers are NOT valid.** Coot's parser
rejects `//A/115,116,121` with `Invalid selection syntax`.

### Multiple disjoint residues: `||`-joined multi-CID

To select several non-contiguous residues on the same chain (a common
need for highlighting catalytic / binding / mutation hotspots), join
single-residue or single-range CIDs with `||`:

```yaml
# Correct: || is the only way to express "these N disjoint residues"
- style: CBs
  selection: "//A/115||//A/122||//A/155||//A/210||//A/221||//A/233-234||//A/246"
  colour: "#e74c3c"
```

The resolver splits on `||` and emits **one representation per chunk**,
each sharing the same style and colour. This is correct, but every
chunk becomes a row in Moorhen's Models drawer. For seven highlights
you get seven rows; for fifty, fifty.

**There is currently no compact form that produces one row.** The
mmdb/gemmi CID grammar simply doesn't have one. If panel-row count
matters for you, group highlights into contiguous ranges where you
can, and accept the row-per-chunk cost where you can't. Improving
this needs a Moorhen API change — see the upstream issue tracker.

## Common authoring pitfalls

A small set of CID and rep-shape mistakes account for almost every
"why isn't anything rendering?" problem. If you're writing a converter
that produces scenes from another tool (PyMOL `.pse`, ChimeraX
session, etc.), check these first.

### Never emit `selection: //`

`//` is **not** a valid CID. Coot reads it as "model nothing, chain
nothing" and the selection comes back empty, so the representation
renders no atoms. The viewer will look broken or near-empty.

To mean "the whole molecule", do one of:

```yaml
# Best: omit `selection:` entirely. The resolver defaults to /*/*/*/*.
- { style: CRs, colour: by-domain }

# Equivalent: explicit wildcard.
- { style: CRs, selection: "/*/*/*/*", colour: by-domain }

# Single chain.
- { style: CRs, selection: "//A", colour: by-domain }
```

**Never** emit `selection: //`, `selection: ""`, `selection: "/"`, or
any partially-empty CID. They all evaluate to "select nothing" and the
rep silently draws nothing.

### `style: ligands` — always emit an explicit selection

Moorhen's `ligands` style does its own non-polymer atom-discovery
inside Coot, which depends on entity types and parsing quirks of the
loaded cif. In practice this misbehaves on cifs from many sources:
it can pick up amino-acid residues (especially glycines) as
"ligands", miss the actual fragment ligand entirely, or — worst —
share an auto-discovered list across all molecules in a multi-load
session, so every fragment in a campaign renders the *same* ligand
neighbourhood instead of its own.

**Always specify the ligand explicitly** by residue name. The CID
form is `//*/(COMPID)`:

```yaml
- file: x0682
  dictionaries: [x0682_LIG]
  representations:
    - style: ligands
      selection: "//*/(LIG)"          # only LIG, by residue name
```

`//*/(LIG)` reads as: any chain, any residue with name `LIG`. The
parentheses are required — they tell the gemmi CID parser this is a
residue-name selector, not a sequence number.

If the fragment comp_id varies per structure (typical in a fragment
campaign — each soaked compound has its own three-letter code), the
converter should read the comp_id out of the dict file's
`data_comp_<X>` line and emit the right name per element. **The
converter knows what the ligand is — there's a dict for it.** Anything
the converter knows should be expressed explicitly; don't trust
Moorhen's auto-discovery for anything load-bearing.

For multiple distinct ligands on one molecule, `||`-join them as
usual:

```yaml
- style: ligands
  selection: "//*/(LIG)||//*/(ADP)"
```

**Why the doc previously said the opposite:** earlier versions of this
file recommended bare `style: ligands` with no selection, relying on
Moorhen's auto-discovery. That advice produced wrong renderings on
real cifs and has been retracted.

### Don't lift `view.clipStart` / `clipEnd` from PyMOL

PyMOL's view state uses massive clip-plane values (often ±10000+ Å)
because PyMOL works in a much larger conceptual volume than Moorhen.
Importing those verbatim into a scene's `view:` block will push
Moorhen's depth-test out so far that fog culling / depth buffer
precision break down and the structure may render as a blank canvas
or with severe z-fighting.

**Either omit `clipStart` / `clipEnd` entirely** (the resolver leaves
Moorhen's defaults of 0 / 1000) **or clamp at author time** to
Moorhen-sensible values. The same applies to `fogStart` / `fogEnd`.

```yaml
# Good: omit, let Moorhen pick defaults
view:
  origin: [2.5, 16.6, -8.6]
  quat: [0.64, 0.67, -0.01, 0.37]
  zoom: 2.5

# Bad: PyMOL-derived values that don't translate
view:
  clipStart: -12586.35       # nope
  clipEnd: 13099.16          # nope
```

`origin`, `quat`, `zoom`, and `background` translate cleanly between
viewers and are safe to copy verbatim.

### Highlights produce one drawer row per CID chunk

The Coot CID grammar does not have a compact form for "these N
disjoint residue numbers" — comma lists of residue numbers are
rejected. The only way to express disjoint residues is a `||`-joined
multi-CID, and the resolver renders each chunk as a separate
representation, each of which becomes a row in Moorhen's Models drawer.

There's nothing the author can do to collapse N highlights into one
row today. The practical advice:

- Group highlights into contiguous ranges where you can.
  `//A/115-117` is one rep / one row; `//A/115||//A/116||//A/117`
  is three.
- Accept the row-count cost otherwise. Twelve catalytic residues
  produce twelve rows.
- If you want a *single* colour applied across a whole multi-highlight
  set, all the chunks share the same colour anyway — the picture
  reads as one logical thing even if the drawer has many rows.

### `colour`

Four forms:

```yaml
colour: "#4b8bbe"             # 1. Hex literal
colour: by-domain             # 2. Compile from the domains: block
colour: b-factor              # 3. Named Moorhen scheme
colour:                       # 4. Raw escape hatch (lossless, ugly)
  raw:
    ruleType: <string>
    args: [...]
    isMultiColourRule: <bool>
```

Named schemes available out of the box:
`by-domain`, `b-factor`, `b-factor-norm`, `af2-plddt`,
`secondary-structure`, `jones-rainbow`, `mol-symm`.

`by-domain` compiles to a single multi-residue colour rule built from
the top-level `domains:` block. So one element with `colour: by-domain`
on a ribbon gives you the whole domain-coloured protein in one rep.

#### Prefer `by-domain` over a bare hex on ribbons (`CRs`)

Observed in practice: a **bare per-element hex** on a `CRs` ribbon does
**not reliably survive an apply / re-export round-trip** — the ribbon
can come back in Moorhen's default chain colours instead of the hex you
asked for. Driving the same colour through the `domains:` block +
`colour: by-domain` is robust, because it compiles to a proper
`isMultiColourRule` colour rule that is applied as one unit.

Rule of thumb:

- **Ribbons / cartoons (`CRs`) you want a specific colour on → use
  `domains:` + `colour: by-domain`,** not a bare hex. To colour a whole
  chain one colour, declare a domain spanning it (an over-wide range like
  `1-1500` is fine — `clamp-and-log` trims it to the residues present):

  ```yaml
  domains:
    - { name: chainA, selection: "//A/1-1500", color: "#1f5fa8" }
    - { name: chainB, selection: "//B/1-1500", color: "#2aa39b" }
  elements:
    - file: complex
      representations:
        - { style: CRs, selection: "//A", colour: by-domain }
        - { style: CRs, selection: "//B", colour: by-domain }
  ```

- **Surfaces (`MolecularSurface`) and stick reps (`CBs`) → a bare hex
  `colour:` took reliably** on re-render. Use hex for highlighted residue
  sticks and surfaces where a domain block would be overkill.
- **`MetaBalls` (and likely `ligands`/`VdwSpheres`) → a bare hex did NOT
  reliably take; the rep fell back to per-atom CPK colouring** (orange C,
  red O, blue N, etc.) on re-render. This is often *desirable* for a
  ligand — CPK reads as "a real molecule with real atoms" — so it may be
  a feature, not a bug. But do not assume a uniform hex will hold on a
  metaball ligand. If you genuinely need a single flat colour on a
  metaball/ligand, expect to verify it in Moorhen and possibly fall back
  to a `colour: { raw: ... }` rule.

When in doubt for a presentation figure: colour structure *chains/domains*
via `by-domain` (reliable), colour *highlighted residue sticks and
surfaces* via hex (reliable), and accept that *metaball/ligand* groups
tend to render **CPK** regardless of the hex you set — usually fine.

### `dictionaries`

A list of dictionary file names (from `files:` with `kind: dictionary`).
The resolver loads each globally first (so the coord parses), then
re-associates with the specific molecule's molNo. This is the key
mechanism for fragment-campaign work: two molecules with same-named
ligands can carry different chemistry because each gets its own scoped
dictionary.

```yaml
elements:
  - file: x0034
    dictionaries: [x0034_LIG, x0034_ADP]      # scoped to this molecule only
    representations: [ ... ]
```

For dictionaries that should apply to *every* molecule (cofactors,
common ions), use the top-level `globalDictionaries:` block instead.

#### Round-trip artifact: drop auto-captured `dictionaries:` for PDB ligands

When you load a deposited entry containing a ligand (e.g. ANP, UUB),
Moorhen auto-acquires a chemical dictionary for that ligand and
associates it with the molecule. On **export**, it faithfully serialises
that loaded state, so the captured scene gains a `dictionaries:` entry
(and sometimes a matching `files:` line) like `dict-complex-ANP` that
**you did not author**.

For scenes whose coordinates come from `pdb:` (or any source where the
ligand ships *with* the structure), **delete these captured
`dictionaries:` references.** Moorhen re-acquires the dictionary on the
next load, so the entry is redundant at best; at worst it's a dangling
name with no backing `files:` entry. Stripping it keeps the scene
portable and matches the "prefer PDB fetch over bundled assets" rule.

**Keep a `dictionaries:` entry only** when the ligand's chemistry is
*not* derivable from the loaded coordinates — the fragment-campaign case
above, where each soaked compound has a custom restraint dict you bundle
in a `.scene.zip`. That's the one situation where the dictionary is
load-bearing rather than a capture artifact.

## `superpose`

Alignments to apply after fetching but before rendering. Each entry
mutates the *moving* structure's display transform in place; the
reference is untouched.

```yaml
superpose:
  # SSM — secondary-structure matching. Cheapest default; needs one
  # chain id on each side. Good for homologues / different conformations
  # of the same chain.
  - { method: ssm, move: holo, onto: apo, movChain: A, refChain: A }

  # LSQ on explicit residue ranges. Use when you know correspondences
  # or want to control which region drives the fit (e.g. align on a
  # binding-site loop only).
  - method: lsq
    move: nod-x
    onto: closed
    matches:
      - { refChain: A, refRange: 104-260, movChain: A, movRange: 104-260 }
      - { refChain: B, refRange: 50-150,  movChain: B, movRange: 60-160 }
    matchType: main      # one of: all | main | ca   (default: main)

  # LSQ convenience shorthand for the common "same chain + same range
  # on both sides" case. Mutually exclusive with `matches`.
  - { method: lsq, move: holo, onto: apo, chain: A, range: 104-260 }
```

The LSQ `matchType` controls which atoms are fitted:

- `all` — every atom in the named residue ranges.
- `main` — main-chain atoms only (default; usually what you want).
- `ca` — Cα atoms only (fastest, most tolerant of side-chain differences).

`gesamt` is not currently supported because the Moorhen build doesn't
expose it; the schema may grow to include it in a later version.

## `view`

The portable subset of Moorhen's view state. Anything you omit is left
alone at apply-time.

```yaml
view:
  origin: [10.5, 20.0, -5.25]            # camera origin
  quat:   [0, 0, 0, -1]                  # camera quaternion
  zoom:   1.5
  clipStart: 0
  clipEnd:   1000
  fogStart:  250
  fogEnd:    1250
  background: "#ffffff"                  # hex
```

## `resolver`

Apply-time policy. Currently one option:

```yaml
resolver:
  onMissingResidues: clamp-and-log   # default
  # onMissingResidues: strict        # fail loudly instead of silently clamping
```

`clamp-and-log` is the right default. The resolver clamps domain
ranges to the residues actually present in the loaded structure, splits
across internal gaps, and writes a sidecar log noting what was modified.
`strict` raises an error when any clamping or splitting would happen
— useful when you're authoring against a structure you know exactly.

## Authoring tips

### Keep the model panel manageable

Every representation entry — and every `||`-split chunk inside one —
becomes a row in Moorhen's Models drawer. For a scene with many
molecules each carrying many highlighted residues, the drawer can
become unwieldy. To keep it tractable:

- **Group highlights into ranges where you can.** `//A/115-117` is
  one chunk / one row; `//A/115||//A/116||//A/117` is three.
- **Use `by-domain` colouring instead of one rep per coloured block.**
  A ribbon with `colour: by-domain` is one row that colours the whole
  protein by the domain map; doing the same with N hex-coloured reps
  produces N rows.
- For repeated patterns across molecules (a fragment campaign), the
  same selection on N molecules costs N rows total, not N×M; you only
  pay for repeats *within* a single molecule.

### Naming bundle assets

Inside the `bundle:` value, use a sensible directory structure:

```
my-scene.scene.zip
├── scene.yaml
└── assets/
    ├── coords/
    │   ├── x0034.cif
    │   ├── x0092.cif
    │   └── ...
    └── dict/
        ├── x0034_LIG.cif
        ├── x0092_LIG.cif
        └── ...
```

The path inside `bundle:` matches the path inside the zip — pick
something stable and human-readable. `assets/` is convention; the
schema doesn't enforce it.

### Multi-comp dictionaries

A single `.cif` file can declare several `data_comp_*` blocks. Coot's
`read_dictionary_string` parses all of them in one call, so the
resolver doesn't need to split. One file ref = all blocks loaded.

### Bare-int ranges

In a domain `selection`, a single residue is just the bare number in the
CID: `"//A/245"` rather than `"//A/245-245"`. Both parse to the same
internal form.

```yaml
domains:
  - { name: catalytic, selection: "//A/245", color: "#e74c3c" }   # one residue
  - { name: helix-A,   selection: "//A/100-130", color: "#4b8bbe" }  # a range
```

The same shorthand applies to the `range:` fields that remain current —
LSQ matches and the `superpose` chain+range form — where YAML's bare
integer (`range: 115`) is accepted and normalised to `"115-115"`:

```yaml
superpose:
  - { method: lsq, move: holo, onto: apo, chain: A, range: 115 }
```

### Disjoint residues on one chain need `||`

The Coot CID grammar does **not** accept comma-separated residue
numbers — `//A/115,116,121` is a parse error. The only valid form for
disjoint residues is `||`-joined single-residue or single-range CIDs:

```yaml
# Correct: || joins single-residue / single-range CIDs
- style: CBs
  selection: "//A/115||//A/116||//A/121-122||//A/146||//A/155||//A/192||//A/198||//A/201||//A/204-205||//A/208||//A/213"
  colour: "#e74c3c"
```

Each `||` chunk becomes its own row in the Models drawer (the
resolver splits and emits one rep per chunk). Group adjacent
residues into ranges where you can to keep the count down.

### Multi-structure scenes share the same domains and view

When you're comparing 25 fragment-soak structures, write the domains
block *once* and one element per molecule that references it:

```yaml
domains:
  - { name: NBD, selection: "//A/104-260", color: "#2ecc71" }
  - { name: HD1, selection: "//A/261-330", color: "#9b59b6" }
  ...

elements:
  - file: x0034
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }
  - file: x0092
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }
  ...
```

For comparing conformations, add a `superpose:` block so the camera
stays meaningful across all of them.

## Worked example

A minimal scene that fetches two PDB entries, aligns them on the NBD,
and renders both with domain colouring:

```yaml
scene: apaf1-nod-aligned
version: 1

files:
  - { name: closed, pdb: 1z6t }
  - { name: nod-x,  pdb: 3sfz }

superpose:
  - { method: lsq, move: nod-x, onto: closed, chain: A, range: 104-260 }

domains:
  - { name: CARD, selection: "//A/1-92",    color: "#4b8bbe" }
  - { name: NBD,  selection: "//A/104-260", color: "#2ecc71" }
  - { name: HD1,  selection: "//A/261-330", color: "#9b59b6" }
  - { name: WHD,  selection: "//A/331-415", color: "#f1c40f" }
  - { name: HD2,  selection: "//A/416-586", color: "#e67e22" }

elements:
  - file: closed
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }
  - file: nod-x
    representations:
      - { style: CRs, selection: "//A", colour: by-domain }

resolver:
  onMissingResidues: clamp-and-log
```

Drop that into the Scenes editor, click Apply, and you have two
APAF-1 structures aligned on their NBDs with matching colour schemes
— ready to compare conformations.
