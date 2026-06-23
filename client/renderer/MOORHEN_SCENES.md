# Moorhen Scenes — implementer notes

A working summary of the `moorhen-scenes` feature — what it does,
which files implement it, the design decisions worth knowing about,
and what's deliberately not done yet.

For the **authoring** reference (how to write a `.scene.yaml`), see
[`types/moorhen-scene.md`](renderer/types/moorhen-scene.md). This
document is the implementer's-eye view.

**Status:** merged to `django` via PR #164 and live in the app. Runs
against **Moorhen 1.0.0-alpha.3** (the `extraSidePanels` callback wiring
in `moorhen-wrapper.tsx` tracks that version's API). MTZ map support
landed after the initial merge; the lifter/resolver capture and re-apply
maps as well as coordinates.

## What it is

A YAML-based scene format for Moorhen. A scene describes *intent*
(domains, colour rules, representations, superpositions, camera) and
optionally bundles its own coordinate and dictionary data, so a
single `.scene.yaml` or `.scene.zip` is enough to reproduce a view on
any machine — including ones that have never seen the underlying
project files.

```
my-view.scene.zip                     ┌─────────────────────────┐
├── scene.yaml         drop on ─────► │  Moorhen Scenes panel   │
└── assets/                            │  • Open / Apply / Save  │
    ├── coords/x0034.cif               │  • Live Monaco editor   │
    ├── dict/x0034_LIG.cif             │  • Validation messages  │
    └── mtz/x0034.mtz                  └─────────────────────────┘
```

The format supports:

- **PDB IDs** (fetched via PDBe proxy), **ccp4i2 file IDs** (per
  project), **URLs**, and **bundled assets** (inside a `.scene.zip`).
- **Per-molecule scoped dictionaries** — the fragment-campaign case
  where two molecules carry same-named ligands with different chemistry.
- **Domain definitions** with by-domain colouring.
- **SSM and LSQ superpositions** (with a `chain`+`range` shorthand).
- **Electron-density maps** from MTZ refs (`kind: mtz` + a `maps:`
  block), including difference maps and an `activeMap:` for refinement.
- **Real-space CCP4 maps and masks** (`kind: map` + a `maps:` entry,
  no columns); masks (`isMask: true`, `File.sub_type ==
  CMapDataFile.SUBTYPE_MASK`) default to a translucent solid surface.
- **Inline dict text** (`cifText:`) and **bundle: assets/...**
  references for sources without a stable URL.

## File map

### Schema and parsing

| File | Role |
|------|------|
| `renderer/types/moorhen-scene.ts`         | TypeScript types for the scene model |
| `renderer/types/moorhen-scene.md`         | Authoring reference (the user-facing doc) |
| `renderer/lib/moorhen-scene.ts`           | YAML parser, validator, serialiser (`parseScene`, `serialiseScene`, `serialiseSceneWithComments`) |

### Resolver (scene → Moorhen state)

| File | Role |
|------|------|
| `renderer/lib/moorhen-scene-resolver.ts`  | The apply-time engine. Fetches files, runs superpositions, applies representations, loads/contours maps, sets the camera. Exposes a couple of small pure helpers (`splitMultiCid`, `clampRangeToPresent`, `expandLsqMatches`, `isFetchable`) for unit testing. |

### Lifter (Moorhen state → scene)

| File | Role |
|------|------|
| `renderer/lib/moorhen-scene-lifter.ts`    | The capture-time engine. `liftScene` writes a plain YAML; `liftSceneToBundle` inlines local coord/dict/MTZ data into a returned asset map suitable for zipping. Captures loaded maps (columns + render state) alongside coordinates. Conservative — only lifts things it recognises. |

### UI

| File | Role |
|------|------|
| `renderer/components/moorhen/moorhen-scenes-panel.tsx`        | The Scenes side-panel (Monaco editor + toolbar + live-apply + validation log + .scene.zip open/save). |
| `renderer/components/moorhen/moorhen-ccp4i2-tabbed-panel.tsx` | Tab container holding Controls + Scenes sub-panels under one Moorhen side-panel registration. (Workaround for a Moorhen quirk where registering two `extraSidePanels` only ever surfaces one tab.) |
| `renderer/components/moorhen/moorhen-wrapper.tsx`             | Wires the resolver + lifter to Moorhen via callbacks: `handleApplyScene`, `handleCaptureScene`, `handleFetchSceneFile`, `handleFetchSceneDictionary`, `handleLoadSceneDictionary`. |
| `renderer/components/moorhen/moorhen-control-panel.tsx`       | Existing panel; only touched to allow co-existence with the Scenes panel. |

### Supporting infrastructure

| File | Role |
|------|------|
| `renderer/app/api/proxy/pdbe/[...path]/route.ts` | Local proxy for PDBe fetches. Adds `Cross-Origin-Resource-Policy: cross-origin` so the COEP-enabled Moorhen page can fetch them. Sibling of the existing `proxy/uniprot/` route. |
| `renderer/components/task/task-elements/fetch-file-for-param.tsx` | Updated to route PDBe fetches through the new proxy (drive-by fix for an unrelated COEP/CORP issue). |
| `client/package.json`                           | New dependencies: `yaml` (parser), `jszip` (bundle round-trip). |

### Tests

| File | Tests |
|------|-------|
| `renderer/__tests__/moorhen-scene.test.ts`          | Schema, validator, serialiser, round-trip stability — including the `maps:` / `activeMap:` block (74 cases) |
| `renderer/__tests__/moorhen-scene-resolver.test.ts` | Pure resolver helpers — clamping, multi-CID splitting, LSQ expansion, isFetchable, chain selectors (28 cases) |
| `renderer/__tests__/moorhen-scene-lifter.test.ts`   | Lifter behaviour, dict enumeration, map capture, round-trip via fake molecules (20 cases) |

**~120 tests**, all passing. Run with `cd client && npx vitest run` for the live count.

The resolver itself (the bit that drives Moorhen) is exercised by
manual verification in the browser — too much runtime state to mock.

## Design decisions worth knowing about

### Two-phase model: format vs runtime

The format itself (parser, validator, serialiser) is pure TypeScript
and has no Moorhen dependency. The resolver and lifter call into
Moorhen but accept callbacks for the parts that need wrapper state
(fetchers, dictionary loader). This kept the test surface tractable
and means the YAML format could be lifted out into a standalone
library if it ever needs to be.

### Bundle (`.scene.zip`) is the portable shape

A `.scene.yaml` works fine when every file ref is portable (PDB ID,
public URL, ccp4i2 fileId in a project the recipient has). For
fragment campaigns and other "local data" cases, the
`.scene.zip` bundles `scene.yaml` + `assets/` so the whole view is
one droppable file. The format uses `bundle: assets/foo.cif` refs;
the resolver looks them up in an in-memory asset map.

### Lifter behaviour: conservative

The lifter recognises a few common cases (named colour schemes,
single-hex colours, by-domain pipe-delimited args, multi-domain
dict ownership) and falls back to a `{raw: ...}` escape hatch for
anything it doesn't understand. The escape hatch is lossless but
ugly; the recognised cases produce clean editable YAML. Round-trip
through `serialiseSceneWithComments(scene, hints)` adds a comment
above each file entry naming the original source.

### Resolver's "scene owns the look"

When the scene has any element for a molecule, the resolver hides
**every** loader-default non-custom representation on that molecule
— not just ones whose style the scene also adds. Without this, a
fragment-campaign scene that has 25 elements with only `style: ligands`
would leave 25 default ribbons visible underneath. The user can
toggle them back via the Models panel if they want.

### One Moorhen side panel, our own tabs inside

Moorhen's `extraSidePanels` API has an upstream quirk: registering
two panels only ever surfaces one tab in the switcher. We worked
around this by registering one `ccp4i2Controls` panel whose content
is a MUI Tabs container hosting the existing Controls UI and the new
Scenes UI as sub-tabs. Both subtrees stay mounted (one is
`display: none`) so Monaco's editor state and the open .zip's
asset map survive tab switches.

### `||`-split for multi-residue selections

Coot's CID grammar (via gemmi) accepts at most a single residue
number or a `start-end` range in the residue field. Disjoint
residues need `||`-joined multi-CIDs. The resolver splits on `||`
and emits one representation per chunk — each becomes a row in
Moorhen's Models drawer. There's no compact form; this is a
fundamental gemmi grammar constraint, not something to fix in the
resolver.

### `style: ligands` is broken for explicit selections

Moorhen's `ligands` representation style overwrites the rep's CID
with its own auto-discovered ligand list on every apply, throwing
away whatever `selection:` the author specified. On real cifs the
auto-discovery misfires (treats glycines as ligands, shares the
same CID across all molecules in a multi-load session). The
authoring doc tells the converter to use `style: CBs` with an
explicit residue-name CID like `//*/(LIG)` for ligand sticks.
This avoids the auto-discovery path entirely.

### Maps are best-effort, never fatal

A `maps:` entry needs the wrapper's map-fetcher to be wired in and an
MTZ ref it can actually fetch. If either is missing, the resolver
**drops the map with a log entry and carries on** — coordinates,
representations, and camera still apply. A scene that can't find its
density never renders as a hard failure. The lifter only emits
render-state fields (`contourLevel`, `style`, colours, …) when the
captured value differs from Moorhen's load-time default, so captured
map blocks stay small. v1 handles MTZ refs only; pre-computed CCP4
`.map` files and PDB-fetched maps are deferred.

### Camera handling: portable subset only

The `view:` block captures origin, quat, zoom, optionally clip and
fog. Lighting / SSAO / shadow / spec power are deliberately omitted
because they don't translate meaningfully across structures and they
bloat the YAML. PyMOL-derived clip values (often ±10000+ Å) should
be stripped at author time — the doc warns about this.

## Known limitations

These are things the format or implementation explicitly does not
do today. Listed in case they come up in conversation:

- **`gesamt` superposition is not exposed.** The WASM module is
  built into Moorhen but the JS binding is commented out. SSM and
  LSQ are available; gesamt is not. Worth filing upstream.
- **No UniProt / SIFTS domain import.** The format supports rich
  domain definitions; there's no automation to populate them from
  UniProt or PDBe SIFTS yet. Hand-author, or generate from your
  own script.
- **Multi-CID highlights produce one drawer row per chunk.** There
  is no compact CID form for disjoint residue numbers, so a
  highlight of 12 disjoint residues becomes 12 panel rows. Could be
  collapsed in a Moorhen UI patch upstream; can't be fixed in our
  format.
- **The lifter doesn't capture the loader's default reps.** It
  captures custom reps the user added via the Scenes apply path,
  but if you've fiddled with the Models drawer toggles, that state
  doesn't round-trip into the saved scene.
- **No project-attached scene storage.** Today scenes live in
  `.yaml` / `.zip` files on the user's local filesystem (or
  download folder). A future iteration could add a `scenes/`
  directory under each ccp4i2 project, exposed via REST.

## End-to-end smoke test

Without the test data we use internally, the cheapest way to verify
the round-trip works on a fresh install:

1. Open Moorhen in the Electron app or web build.
2. Pick a structure (any PDB ID — `3jbt` is small).
3. Open the Scenes panel (CCP4i2 side panel → Scenes tab).
4. Paste this and click Apply:

   ```yaml
   scene: smoke-test
   version: 1
   files:
     - { name: thing, pdb: 3jbt }
   domains:
     - { name: a, chain: A, range: 1-50,   color: "#4b8bbe" }
     - { name: b, chain: A, range: 51-100, color: "#e74c3c" }
   elements:
     - file: thing
       representations:
         - { style: CRs, selection: "//A", colour: by-domain }
   ```

5. Should see a structure load from PDBe and render with two
   coloured blocks on chain A.

## How it landed, and the Materia pin

The feature merged to `django` via **PR #164**
(`moorhen-scenes-into-django`). Follow-up commits added hex-alpha
colours, an ATP-class dict skip, and MTZ map capture/apply. The
Electron menu addition in `client/main/ccp4i2-menu.ts` is the only
main-process touch; everything else is renderer-side.

Two existing files are touched non-trivially:

- `client/renderer/components/moorhen/moorhen-wrapper.tsx` —
  the scene callbacks and the bundle/dict/map fetcher wiring.
- `client/renderer/components/task/task-elements/fetch-file-for-param.tsx` —
  routes PDBe fetches through the `/api/proxy/pdbe` route (not
  scenes-specific; a drive-by COEP/CORP fix).

Materia tracks scene work through its submodule pin — see
[`packages/ccp4i2-api/RELEASING.md`](../../packages/ccp4i2-api/RELEASING.md)
for the cherry-pick/branch convention when shipping scene changes
that Materia needs ahead of a `django` merge.
