---
name: moorhen-scenes
description: Author, validate, apply, lift, or change the schema of Moorhen scenes, the YAML descriptions of a molecular view used by CCP4i2's embedded Moorhen viewer. Use when writing a scene for a report, campaign summary, test fixture or Moorhen page; when editing the scene Zod schema or its generated contracts; or when a scene does not apply as expected.
---

# Moorhen scenes

A scene is a portable YAML document: which files to load, which domains to
recognise, what representations and colours to apply, which maps to contour,
and where the camera sits. The viewer's CCP4i2 side panel has a Scenes tab that
applies a scene to the live view and lifts the live view back into one.

This skill routes to the sources of truth and holds only the judgement the
schema cannot express. Do not copy the grammar into new documents: it is
generated and CI-gated, and copies drift.

## Where the truth lives

| What | Where | Status |
|---|---|---|
| The format (Zod) | `client/renderer/lib/scene/core.ts` (portable core), `dialect.ts` (ccp4i2 refs) | **source of truth** |
| Public API | `client/renderer/lib/scene/index.ts`: `parseScene`, `validateScene`, `serialiseScene`, `buildJsonSchemas` | hand-written |
| JSON Schemas | `client/renderer/lib/scene/moorhen-scene.{core,ccp4i2,structured}.v1.json` | generated, committed |
| Grammar document | `client/renderer/types/moorhen-scene.md` | generated |
| LLM system prompt | `client/renderer/lib/scene/moorhen-scene.system-prompt.v1.md` | curated prose + generated grammar |
| Compact brief for prompts | `client/renderer/lib/scene/brief.ts` (`buildSceneBrief`) | generated from the schema |
| Server mirror (slim, frontend-free deployments) | `server/ccp4i2/scene_contracts/` | generated copy; never edit |
| Decision record | `client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` | read before changing shape |

Read the grammar document first when authoring; read `core.ts` when changing
the format.

## Authoring a scene

1. Get the contents right before the style: chains, ligand residue names and
   CIDs come from the files, never from memory. In the app the prompt builder
   (`lib/moorhen-scene-prompt.ts`: `buildAuthoringPrompt`) assembles a contents
   summary and a project manifest for exactly this reason.
2. Write the YAML against the grammar document. `version: 1`.
3. Validate: `validateScene(yamlText)` from `lib/scene`; in the Scenes tab the
   editor shows the same errors as markers.
4. Apply in the viewer (Scenes tab) and check the picture, then lift the live
   view back to YAML to see what the viewer actually honoured.
5. Save as `.scene.yaml`, or `.scene.zip` when files are bundled. The panel
   downloads to the browser; scenes are not yet persisted server-side.

### Rules the schema cannot express

- **Selections are Coot CIDs**: `//A`, `//A/703-740`, `//*/(LIG)`,
  `//A/750/CA`. Join several with `||` (`//A||//B`), in representation
  selections and in `view.centre` / `view.slab` alike.
  **A residue *name* needs parentheses** — `//*/LIG` is a parse error, and a
  bad selection draws nothing rather than complaining, so the mistake is
  silent. The asymmetry that invites it: a comma list is valid *inside* the
  parens (`//*/(AY7,LIG)`) but not for residue numbers (`//A/115,116` is a
  parse error). Verified against gemmi.
- **A representation draws its own `selection`**, the whole molecule if
  omitted. Colour never limits what is drawn; scope the selection instead.
- **Colour forms**: hex `#rrggbb`; a named scheme (`by-domain`, `b-factor`,
  `af2-plddt`, `secondary-structure`, `jones-rainbow`, `mol-symm`); or a
  per-selection list. `by-domain` needs a top-level `domains:` block.
  `elements[].colour` is a molecule-wide default that a representation's own
  `colour` overrides.
- **Camera and depth are separate**: `view.centre` moves the camera,
  `view.slab` only sets the clip depth. To centre on and slab to a selection,
  give both.
- **MTZ maps need `columns`**; a bare `files[]` mtz entry has no columns and
  will not load. Real-space maps take no columns. `isDifference` selects the
  two-colour contour; `isMask` marks a mask, which is the same file type as a
  map distinguished only by sub-type.
- **Masking a map** (`maskMaps[]`) takes its source from a `maps[]` entry name,
  because that carries the columns, and renders through a second `maps[]`
  entry with `from:`. `keep: inside` (default) keeps density near the
  selection. Only reach for masking when the request is about carving density
  to a region.
- **Dictionaries belong to a molecule, and the molecule's job decides which.**
  The rule: a molecule's dictionaries are the ones its *job* took as input or
  wrote as output, and nothing else. A project routinely holds several
  ligands called LIG or DRG, so associating by residue name, or by "any
  dictionary in the project", is a guess that is wrong as soon as there are
  two. The server answers it from the database (`jobs/{id}/dictionaries/`,
  `files/{id}/companion_dictionaries/`); do not reconstruct it on the client.
  In a scene, put each dictionary in `files[]` with `kind: dictionary` and
  list it under its molecule's `elements[].dictionaries`. On apply it is
  attached to that molecule as it loads and is **never** loaded into Coot's
  global store, because a global entry is inherited by every molecule that
  has none of its own. `globalDictionaries` is for a monomer you really do
  want shared by every molecule. A dictionary listed in `files[]` but
  attached nowhere still goes global, which older hand-written scenes relied
  on; do not author new scenes that way.
- **Honoured versus hint.** Anything with a physical unit, or that governs
  visibility (`alpha`, clip planes, geometry radii, camera, background) is
  honoured: a renderer must reproduce it. `hints` (lighting, SSAO, edge
  detect, depth blur, shadows) are advisory; never rely on them for what must
  be visible. Test for a new field: if a renderer ignored it, would the image
  be *wrong* or merely *plainer*? Wrong means core with a unit; plainer means
  a hint.

### File references

| Ref | Meaning | Portable? |
|---|---|---|
| `pdb` | PDB id, fetched through the proxy | yes |
| `url` | absolute URL | yes |
| `bundle` | asset inside a `.scene.zip` | yes |
| `cifText` | inline dictionary CIF | yes |
| `fileId` | ccp4i2 file id; fetches with no project qualifier | this deployment only |
| `job` + `param` (+ `projectId` or `projectName`) | a job's output by role | this deployment only |
| `relativeUrl` | origin-relative loader URL | this deployment only; never author it |

Prefer `job` + `param` or `fileId` inside a project: they survive project
moves. A strict-portable export lowers deployment refs to `bundle` or `url`.

## Changing the format

1. Edit `core.ts` for anything portable and Moorhen-meaningful, `dialect.ts`
   for ccp4i2-only references. Decide the conformance class of every new
   field with the honoured-versus-hint test above and put it in the right
   layer.
2. Regenerate every derived artefact in one go, from `client/`:

   ```bash
   UPDATE_SCHEMA=1 npx vitest run renderer/__tests__/scene-schema.test.ts
   ```

   This rewrites the three JSON Schemas, the grammar document, and the server
   mirror. Without `UPDATE_SCHEMA` the same test fails if any of them is
   stale, which is the CI gate. Commit the regenerated files with the change.
3. Run the scene tests: `npx vitest run renderer/__tests__/moorhen-scene*`,
   plus `scene-schema.test.ts` and `renderer/__tests__/fixtures/demo.scene.yaml`,
   which every fenced YAML example must keep parsing.
4. Update the lifter (`lib/moorhen-scene-lifter.ts`) and resolver
   (`lib/moorhen-scene-resolver.ts`, `applyScene`) together: a field the
   resolver honours that the lifter cannot capture breaks the
   lift-then-apply round trip.
5. Record a decision that changes the shape in the design note. Keep
   `version: 1` unless the change is breaking.

## Producers and consumers

- **Client**: `lib/moorhen-scene-resolver.ts` applies (`planDictionaryScopes`
  decides what may go global), `lib/moorhen-scene-lifter.ts` captures,
  `lib/moorhen-dictionaries.ts` attaches dictionaries per molecule, `components/moorhen/moorhen-scenes-panel.tsx` is the UI,
  `lib/moorhen-scene-prompt.ts` and `components/moorhen/use-scene-nl-capability.ts`
  drive natural-language generation when a deployment provides it.
- **Server**: `lib/campaign_scene.py` (`build_summary_scene`) and the
  `summary_scene` action on `ProjectGroupViewSet` emit ccp4i2-dialect scenes
  with `fileId` references. A frontend-free server reads the contracts as
  package data from `scene_contracts/`.
- **Not this format**: the report `Picture` element (`report/pictures.py`)
  and `docs/scene-files.md` describe the older CCP4mg scene XML and its
  migration path.

## Traps

- A captured scene whose ligands draw with the wrong bonds after a paste
  almost always means an element with no `dictionaries:` beside another
  file's ligand of the same residue name. Every loader records where a
  molecule's dictionaries came from (`lib/moorhen-dictionaries.ts`), and
  capture writes them from that; a new load path that attaches dictionaries
  any other way, or reads one into the global store, breaks the round trip.

- Loading a P1 map with 90° angles from MTZ coefficients marks it EM in
  Moorhen; the viewer needs `primeEmMapHeaderInfo` (in `lib/moorhen-map-file.ts`)
  before the map is added or Moorhen's origin listener crashes the page.
  Every MTZ load site already calls it; keep that when adding one.
- `mol.uniqueId` is always the loader URL; the lifter derives file references
  from it, so a molecule loaded by a route it does not recognise lifts as a
  `relativeUrl`.
- The scenes panel's Generate tier only appears when the deployment reports
  the capability; on the desktop it is absent by design.
