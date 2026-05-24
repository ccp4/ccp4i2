# Moorhen Scene format — v1

A **scene** is a portable, human-editable YAML description of how to
render one or more structures in Moorhen. It captures *intent*
(domains, colour rules, representations, superpositions, camera)
separately from any specific PDB file, so the same visual treatment
can be re-applied across different structures of the same protein.

This document is the authoring reference — written for a person (or a
script) building a scene from outside Moorhen. The TypeScript types
live in [`moorhen-scene.ts`](./moorhen-scene.ts); parse/serialise
logic lives in [`../lib/moorhen-scene.ts`](../lib/moorhen-scene.ts).

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
  - { name: CARD, chain: A, range: 1-92,    color: "#4b8bbe" }
  - { name: NBD,  chain: A, range: 104-260, color: "#2ecc71" }
  - { name: HD1,  chain: A, range: 261-330, color: "#9b59b6" }
```

Fields:

- `name`: required.
- `chain`: required. Three forms:
  - `"A"` — a single chain.
  - `"*"` — every chain present in the structure (useful for symmetric
    assemblies like the apoptosome heptamer).
  - `["A", "B", "C"]` — explicit list, applied to each chain in turn.
- `range`: required. Two forms:
  - `"start-end"` (string) — inclusive range, e.g. `"100-200"`.
  - `<integer>` (bare int) — single residue; the validator normalises
    `115` to `"115-115"` so downstream code only sees one shape.
- `color`: required hex `#rrggbb`.

The resolver clamps each range to the residues actually present in the
loaded structure (see `resolver.onMissingResidues`).

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

Any [Moorhen representation
style](https://github.com/moorhen-coot/Moorhen/) (27 supported as of
v1). The common ones:

| Style              | What it draws                       |
| ------------------ | ----------------------------------- |
| `CRs`              | Cartoon ribbons.                    |
| `CBs`              | All-atom sticks (carbon bonds).     |
| `CAs`              | Cα-only trace.                      |
| `MolecularSurface` | Smooth molecular surface.           |
| `VdwSpheres`       | Van-der-Waals spheres.              |
| `ligands`          | Stick representation of HET atoms.  |
| `DishyBases`       | Cartoon-style nucleotide bases.     |
| `allHBonds`        | Hydrogen bonds.                     |

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

Use bare integers for single residues. `range: 115` is more readable
than `range: "115-115"`, and parses to the same internal form. Same
applies inside LSQ matches:

```yaml
domains:
  - { name: catalytic, chain: A, range: 245, color: "#e74c3c" }   # bare int
  - { name: helix-A,   chain: A, range: 100-130, color: "#4b8bbe" }  # range
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
  - { name: NBD, chain: A, range: 104-260, color: "#2ecc71" }
  - { name: HD1, chain: A, range: 261-330, color: "#9b59b6" }
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
  - { name: CARD, chain: A, range: 1-92,    color: "#4b8bbe" }
  - { name: NBD,  chain: A, range: 104-260, color: "#2ecc71" }
  - { name: HD1,  chain: A, range: 261-330, color: "#9b59b6" }
  - { name: WHD,  chain: A, range: 331-415, color: "#f1c40f" }
  - { name: HD2,  chain: A, range: 416-586, color: "#e67e22" }

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
