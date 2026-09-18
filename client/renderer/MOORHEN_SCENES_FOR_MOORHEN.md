# Moorhen scenes — technical briefing for Moorhen developers

A *scene* is a portable, validated YAML description of a Moorhen view (which
files to load, what to draw, how to colour it, lighting, where the camera sits).
ccp4i2 uses it in production. This note is **informational** — it shows the data
model, how it's built and kept honest, and (importantly) that every part is
grounded in Moorhen's *own* API, so you can see where it would sit. It is not an
ask; the narrative case is in `MOORHEN_SCENES_UPSTREAM_PROPOSAL.md`.

## The package (artefacts to read, in order)

| Artefact | What it is |
|---|---|
| `lib/scene/moorhen-scene.core.v1.json` | The **published contract** — a JSON Schema, no ccp4i2 concepts. The thing to validate against / target. |
| `types/moorhen-scene.md` | Human-readable grammar + examples (generated from the contract). |
| `lib/scene/core.ts` | The Zod **source of truth** the contract is generated from. |
| `__tests__/fixtures/demo.scene.yaml` | A runnable example (`pdb: 1B9K`) — validates under the core schema, no project needed. |
| `MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` | The design rationale (layering, conformance model, the grounding findings below). |

`lib/scene/moorhen-scene.ccp4i2.v1.json` is the same core **plus** ccp4i2's
deployment-specific reference kinds — included only to show the extension seam.

## The data model

Three layers, by who they concern and how strictly a renderer must honour them:

| Layer | Contents | Conformance | In the core (upstreamable) schema? |
|---|---|---|---|
| **Core** | portable file refs (`pdb`/`url`/`bundle`/`cifText`), domains, selections, representations + **honoured geometry**, maps, camera/view, superpose | MUST honour | ✅ |
| **Hints** | scene-global **lighting** + perceptual **effects** | MAY ignore (advisory) | ✅ (optional) |
| **ccp4i2 dialect** | `fileId`/`job`+`param`/`relativeUrl` (deployment-coupled) | — | ❌ — injected extension |

The honoured-vs-advisory split is deliberate: a headless or constrained renderer
honours geometry/colours/camera and may drop lighting/effects, still producing a
*correct* image. "Honoured" fields carry physical units (Å) or govern visibility
(opacity); "advisory" fields have no ground truth (lighting direction, SSAO).

## It's grounded in Moorhen's API (not invented)

Every part maps onto Moorhen's existing surface — the resolver drives these
directly; the lifter reads them back:

| Scene field | Moorhen API |
|---|---|
| `representations[].style` | `RepresentationStyles` union (`MoorhenMoleculeRepresentation.d.ts`) — the schema enumerates all 29 |
| `representations[].geometry` (bond/ball/probe radii, ribbon widths, Å) | `m2tParameters` (`setM2tParams` + `useDefaultM2tParams`) |
| `representations[].colour` | representation colour rules (`ColourRule` objects assigned to `representation.colourRules`) + the named multi-colour rule types, whose data comes from `getMultiColourRuleArgs` |
| `hints.lighting` (direction, ambient/diffuse/specular, shininess) | `sceneSettingsSlice` (`glRefSlice` before 1.0.1): `lightPosition`, `ambient`/`diffuse`/`specular` (`[r,g,b,a]`), `specularPower` |
| `hints.effects` (ssao, edgeDetect, shadows, depthBlur, perspective) | `setDoSSAO`/`setDoEdgeDetect`/`setDoShadow`/`setUseOffScreenBuffers`/`setDoPerspectiveProjection` |
| `view` (origin, quat, zoom, clip/fog, background) | `sceneSettingsSlice` for origin, zoom and clip/fog, `glRefSlice` for `quat` (all of them `glRefSlice` before 1.0.1) + `setBackgroundColor` |

## How it's built and kept honest (the mechanism)

One source, many drift-free artefacts:

```
lib/scene/core.ts (+ dialect.ts)   ← Zod, the single source of truth
        │  z.toJSONSchema / generators (CI-gated: regenerate must match)
        ├─► moorhen-scene.{core,ccp4i2}.v1.json   the published contract
        ├─► types/moorhen-scene.md                human grammar + examples
        └─► a compact LLM authoring brief (~1.3k tokens)
```

- **Validation** is the schema (`parseScene`/`safeParse`); errors carry a path,
  usable for a generate→check→repair loop.
- The **resolver** (scene → Moorhen) takes the host's fetchers as callbacks and
  otherwise touches only molecule/map/representation/scene-settings APIs — no CCP4
  dependency. The **lifter** is the inverse (Moorhen state → scene), emitting
  only non-default values.
- **Apply policy**: `onMissingResidues: clamp-and-log` clamps a slightly-wrong
  residue range to what's present and logs it, instead of failing. Effects are
  *scene-authoritative* (an `effects` block fully determines effect state).

## Findings in the Moorhen codebase (may be useful to you)

While grounding the format we noticed a few things, shared in case they're news:

- **`outline` (`doStenciling`) is a hover/selection highlight, not a model
  silhouette.** `setOutlinesOn` only outlines buffers with `doStencil` set, which
  is the `isHoverBuffer` perfect-spheres path
  (`MoorhenMoleculeRepresentation.ts`). So it does nothing on a static view with
  nothing picked — we excluded it from the scene format for that reason. The
  prominent line effect is `doEdgeDetect`.
- **`useOffScreenBuffers` is the depth-blur enable** ("historically named", per
  the comment in `mgWebGL.tsx`) — not a general off-screen-buffer gate.
- **`setDoOutline` updates the store but the view doesn't re-render** when driven
  programmatically (we confirmed `sceneSettings.doOutline` tracks the dispatch,
  but the display only changed via a manual toggle). The `[doOutline]` effect in
  `MoorhenWebMG.tsx` is structurally identical to the working `[doSSAO]` one, so
  this looks like a render-path quirk worth a glance if you care about
  programmatic outline.
- **`lightPosition` is a position** (default `[25,25,50,1]`, |·|≈61), not a unit
  direction — the resolver maps a conceptual `direction` to a position along it.
- **Colour rules have two levels, and since 1.0.1 they are properly separate.**
  A representation uses the parent molecule's `defaultColourRules` until it is
  given rules of its own; assigning `representation.colourRules` copies the list
  and turns the defaults off, so per-representation colours no longer leak. (Up
  to 1.0.0-beta.1 the list was shared by reference and `addColourRule` pushed
  onto it molecule-wide; we nulled `colourRules` first to get a private list.)
  The scene models the same two levels (molecule + per-rep). One gap remains:
  `representation.addColourRule` has no parameter for `multiColourData`, so a
  multi-colour rule cannot be made through it. We build `ColourRule` directly
  (seventh constructor argument) and assign the list.
- **SSM/LSQ superpose moves coordinates inside coot** (`changesMolecules` +
  `setAtomsDirty`), with no display transform — so superposition can't be
  captured from a molecule's state, only re-applied from a remembered directive.

## What carving out a core would involve

The ccp4i2-specific surface is narrow: a couple of reference *kinds*
(`fileId`+`projectId`, `job`+`param`, `relativeUrl`), the fetcher
implementations, and the routine that recovers a file id from a loader URL. A
core would keep the portable kinds (`pdb`/`url`/`bundle`/`cifText`) and turn the
deployment-specific ones into injected extension points — the same seam where an
application plugs in its data backend. ccp4i2 would then consume the core and
register its own kinds, proving the extension API holds.
