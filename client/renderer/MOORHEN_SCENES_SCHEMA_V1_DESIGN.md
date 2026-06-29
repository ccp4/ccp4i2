# Moorhen Scene format — v1 design decisions

Status: **agreed, pre-implementation.** This records the decisions for hardening
the Moorhen scene format into a strict, published, versioned schema. The schema
is intended to change *infrequently*, so the rationale is captured here, not just
the result. Companion docs: `types/moorhen-scene.md` (grammar), `MOORHEN_SCENES.md`
(implementer notes), `MOORHEN_SCENES_UPSTREAM_PROPOSAL.md` (upstreaming case).

## 0. Context that frames every decision

- The format is meant to be **published and eventually upstreamed to Moorhen**, and
  consumed by editors and LLMs — so the contract must be a machine-readable,
  language-agnostic artifact, not TS-internal validation logic.
- **Backward compatibility is a non-issue.** The feature has had a single user
  (the author) and there are no persisted `.scene.yaml` / `.scene.zip` artefacts
  worth preserving. Therefore **all breaking changes land now, before any freeze** —
  the cost of breaking changes only rises once there is a second user or a frozen v1.
- "No data to migrate" is *not* "code may be inconsistent": the lifter, resolver and
  LLM prompt all reference these fields and must change in lockstep.

## 1. Source of truth: Zod → generated, committed JSON Schema

Chosen over (a) generating JSON Schema from the TS interfaces and (b) hand-authoring
the JSON Schema. Rationale: a strict published format's enemy is **drift** between
its representations (types ↔ validator ↔ schema ↔ docs). Zod minimises that to a
**single** hand-maintained source:

- TS types become `z.infer<typeof Scene>` — compiler-proven identical to validation.
- Runtime validation *is* the schema (`.parse`/`.safeParse`) — the validator cannot
  drift from the types because it is the same object. The hand-written
  `validateScene` in `lib/moorhen-scene.ts` is **deleted**.
- Cross-reference semantics (referential integrity no JSON Schema can express —
  `activeMap` names a `maps:` entry; superpose `move`/`onto`/`file` name `files:`
  entries; `globalDictionaries` name `kind: dictionary` refs) live in
  `.superRefine()`, *attached to the schema*.
- The JSON Schema is **generated and committed** to the repo (Zod 4 `z.toJSONSchema()`
  if available, else `zod-to-json-schema`), with a CI test that regenerates and
  asserts byte-equality. The committed `moorhen-scene.v1.json` is the published
  contract: stable `$id` encoding the version, referenced by `.scene.yaml` via
  `# yaml-language-server: $schema=…` for live editor validation, and handed to LLMs
  as a structured-output constraint.

Flip condition (recorded, not chosen): if the format's centre of gravity becomes the
external spec that third parties implement against and its exact wording is the
deliverable, move to schema-first (hand-authored JSON Schema → TS types generated via
`json-schema-to-typescript`). Zod-as-source captures ~95% of that benefit at far
lower maintenance cost.

## 2. Validation strictness

With compat off the table, the parser **rejects unknown keys** (`.strict()`), where
today's hand validator silently ignores them. For a strict format this is the better
default — it catches typos (`colur:`) at author time. Two **profiles**:

| Profile | Use | Behaviour |
|---|---|---|
| `permissive` | default read/apply | structural validation only |
| `strict-portable` | export / publish / share | also rejects `relativeUrl` and any unlowered deployment ref |

Portability is thus a *checkable property*, not a per-field guess.

## 3. Three-layer stratification

The lever that keeps the **core** stable while volatile concerns move freely.

| Layer | Contents | Conformance | Upstream core? |
|---|---|---|---|
| **Strict core** | portable file refs, domains, selections, representations + honoured geometry, maps, view geometry, superpose | MUST honour; validated | ✅ |
| **`hints`** | the lighting model (principal directional light + per-rep material) and perceptual post-processing — SSAO/ambient occlusion, edge detection, depth-of-field blur, depth cue, shadows | MAY ignore; advisory | optional |
| **ccp4i2 dialect** | symbolic `fileId`/`job`+`param`/`projectId`; `relativeUrl` | resolve / lower on export | ❌ (separate `$id`, `allOf`-imports core) |

In Zod terms: one `CoreScene` object; the dialect is `CoreScene.extend({…})` in a
separate module; the published core JSON Schema stays uncontaminated by ccp4i2 concepts.

## 4. The conformance rule (honoured vs hint)

The dividing axis is **conformance**, not aesthetics. Test, to be written into the spec
so future additions classify consistently:

> If a conforming renderer ignored this property, would the image be **wrong** or merely
> **plainer**? Wrong → honoured core (and give it a unit). Plainer → hint.

Corollary: **a property carrying a physical unit (Å, …) is must-honour by definition** —
units are the tell.

- **Honoured (core):** bond/atom radius, ribbon & arrow dimensions, clip/fog planes (Å),
  background, camera, **opacity (`alpha`)**. These have objective ground truth, *or*
  govern what is **visible**; a triangle generator must undertake to honour them. (Bond
  width is *not* a hint — it has a real radius in Å.)
- **Hints:** the lighting model (light direction/colour/intensity, per-rep material
  reflectance) and perceptual post-processing (SSAO/AO, edge detect, DoF blur, depth cue,
  shadows) — no physical ground truth; ignoring them yields a correct-but-plainer render.

### 4a. Three conformance classes (documentation, not a schema layer)

Light direction exposed a gap in the two-bucket model. It has no ground truth (so it is
not honoured) yet has no meaningful "off" state (every renderer always has *a* light), so
"MAY ignore" is the wrong verb — a renderer that can't honour it **substitutes its
default**, it doesn't omit it. Three classes, recorded as guidance:

| Class | Ground truth? | Unsupported renderer… | Examples |
|---|---|---|---|
| **Honoured** | yes (units), or governs visibility | must reproduce | bond/atom radius, clip planes, `alpha` |
| **Substituted** | no, but always-present | falls back to its own default | light direction, lighting intensities |
| **Additive hint** | no, has an "off" | renders without it | SSAO, edge detect, DoF, shadows |

Both *substituted* and *additive* live in the `hints` layer; the distinction is
documented, not structural. Lighting is flagged as **the first hint a faithful renderer
should honour**.

### 4b. Lighting and materials

**Grounding finding (overrides the earlier per-rep-material plan).** Inspecting Moorhen's
`glRefSlice`: the scene has a single light with `lightPosition: [n,n,n,n]` and
`ambient` / `diffuse` / `specular` as **`[r,g,b,a]` colours of that light**, plus
`specularPower`. Moorhen has **no per-object material concept at all** — diffuse/specular
are properties of the *light*, not of a surface. The earlier "per-rep material
(diffuse/specular/shininess)" plan assumed those were material reflectances; they are not.

**v1 decision (implemented):** model lighting faithfully as Moorhen's single scene-global
light — `hints.lighting` = `{ direction, ambient, diffuse, specular, shininess }`
(direction → `lightPosition`, shininess → `specularPower`; colours as hex). Zero fiction;
every field maps to a real setter. `direction` is the *substituted* class.

**Per-representation material is DEFERRED** (not in v1). The glossy-surface-over-matte-cartoon
goal is real, but with no Moorhen material concept it would be a `hints` field that maps to
nothing — the fiction the hints layer exists to avoid. Sequencing: add a material concept to
Moorhen upstream first (the author has commit access), then spec `SceneRepresentation.material`
against a real target.

**Opacity stays core** regardless: `alpha` governs what is *visible* (rendering a translucent
surface opaque *occludes* the model it was meant to reveal — closer to *wrong* than *plainer*),
unlike shading, which is advisory.

**`direction` → `lightPosition` convention (implemented + render-confirmed 2026-06-29).**
Moorhen's `lightPosition` is a *position* (default `[25,25,50,1]`, |·|≈61, w=1), not a unit
vector. The resolver (`directionToLightPosition`) normalises the scene's conceptual
`direction` and places the light at distance 60 along it, w=1 — so authors think in
*direction*, not scene-unit coordinates. Axes are view-space: **+x = right, +y = up,
+z = toward the viewer**. Confirmed in a real render: `[0,0,1]`→head-on, `[1,1,1]`→upper-right-front.
SSAO reads best on `MolecularSurface` representations (per Martin), so surface reps make the
`effects.ssao` hint visible.

**Scope:** one principal directional light + a per-rep Blinn-Phong material now. **Defer**
multiple lights and full PBR (roughness/metalness) — no authoring need yet, and
must-honour-grade lighting that nothing renders is exactly the fiction the `hints` layer
exists to avoid.

## 5. File reference model

`mol.uniqueId` is always the loader URL. The format's refs are:

| Ref | Meaning | Portability | Layer |
|---|---|---|---|
| `pdb` | PDB id, fetched via proxy | portable | core |
| `url` | absolute URL (`https://…`) | portable | core |
| `relativeUrl` *(was `path`)* | origin-relative URL (`/api/…`) | within-deployment only | dialect |
| `bundle` | asset inside a `.scene.zip` | portable | core |
| `cifText` | inline dictionary CIF | portable | core (dictionaries) |
| `fileId` (+`projectId`) | symbolic project file | within-deployment only | dialect |
| `job`+`param` (+`projectId`) | symbolic job output | within-deployment only | dialect |

- **Rename `path` → `relativeUrl`.** It holds the origin-relative loader URLs (e.g.
  `/api/proxy/pdbe/…`, `/api/ccp4i2/download?…`) — the ones that fail the `^https?://`
  test and fall through today (lifter:678/719). The lifter's `url`-vs-`path` split is
  literally absolute-vs-relative URL; the rename makes that branch self-documenting.
  `path` was actively misleading (reads as filesystem in JS/web).
- **Lowering on export:** a `strict-portable` export rewrites deployment refs
  (`relativeUrl`, `fileId`, `job`+`param`) into portable refs (`bundle` asset, or
  absolute `url` if the server is publicly reachable). Symbolic dialect refs are
  "recipes for producing a core ref."

## 6. Lockstep code changes (no data migration, but code must stay coherent)

- `lib/moorhen-scene.ts` — replace hand validator with Zod parse/validate/serialise.
- `lib/moorhen-scene-lifter.ts` — emit `relativeUrl` not `path` (lifter:678,719); the
  bundler's `delete ref.path` (426,551) → `delete ref.relativeUrl`.
- `lib/moorhen-scene-resolver.ts` — consume `relativeUrl` where it consumes `path`.
- `lib/moorhen-scene-prompt.ts` — update the "never use path:" guidance to the new name.
- `types/moorhen-scene.ts` — becomes `z.infer` re-exports (or is regenerated).

## 7. Testing strategy — two tiers, one fixture corpus

The corpus is the multiplier for "many many tests"; it is *not* a backward-compat net
(there's nothing to preserve) — it is the forward conformance + doc-lint suite.

- **Tier 1 — format (vitest, fast, no browser):** a `fixtures/valid/*.scene.yaml` and
  `fixtures/invalid/*.scene.yaml` corpus (each invalid paired with expected
  path/message), driven table-style against the Zod schema. Plus **property tests**:
  parse→serialise→parse is stable; lift→resolve→lift reaches a fixed point.
- **Tier 2 — behaviour (Playwright, slow):** one parameterised spec iterating the same
  corpus, applying each scene in the real Moorhen runtime and asserting on
  **serialisable Moorhen state** (molecule count, representation list, active map,
  applied superpose transforms) rather than pixels. Visual-regression screenshots for a
  *handful* of canonical scenes only — headless WebGL pixel diffs are the flaky tier.
- **Docs consistency:** a doc-lint test asserts every fenced ```yaml block in
  `types/moorhen-scene.md` parses and validates against the schema, so docs cannot drift
  from the format. **Superposition** (currently under-exemplified vs its rich type) gets
  first-class coverage across docs + fixtures + both tiers.

## 7a. AI authoring & prompt compactness (front-and-central goal)

AI authoring is a primary use case; the authoring prompt should be compact enough to
*aspire* to fit an on-device / OS model (Apple Intelligence, Windows Copilot). Key
correction: **the JSON Schema is token-heavy, not compact** — generated JSON Schema (with
its `$defs`/`$ref`/boilerplate) is typically *larger* than a terse grammar, so it must not
be used as prompt text. Separate the two artefacts by audience:

| Artefact | Form | Audience | Optimised for |
|---|---|---|---|
| `moorhen-scene.v1.json` | verbose JSON Schema | machines / **constrained decoders** | precision |
| **LLM authoring brief** *(generated)* | terse TS-signature + curated examples | model **context** | token economy |
| `types/moorhen-scene.md` | prose | humans | learnability |

All three generate from the **one Zod source** → drift-free.

Two delivery paths:

1. **Capable models** — use the JSON Schema for **constrained decoding** (fed to the
   decoder, not the prompt): structure guaranteed, ~zero prompt tokens spent on shape,
   tokens spent on *intent*. Format is YAML but constrained decoding is JSON-shaped, so the
   model emits **JSON → convert to YAML** on save.
2. **On-device / small models** — no custom-schema decoding and tiny context, so spec must
   live in-prompt; compactness is everything. **Tiered / progressive prompting**, which the
   §3 stratification gives for free: a minimal **core mini-grammar** (files + a couple of
   representations + colour + view; target ≲1k tokens) for on-device, with
   superpose/maps/hints **escalating** to a larger model on demand.

Format-level consequences (compactness as a design constraint): terse guessable keys,
strong defaults so most fields are omitted (mirrors the lifter's emit-only-non-defaults
rule), shallow nesting, and an explicit small **authoring-core subset** that hides escape
hatches (colour `raw`, etc.) a small model would fumble. The generated brief's token count
is **measured**, so "fits on-device" is a number, not a hope. `lib/moorhen-scene-prompt.ts`
becomes a generated compact brief, not hand-maintained prose.

## 8. Recommended implementation order

1. **Zod structural schema** — core (portable) + ccp4i2 dialect; generate + commit both
   `moorhen-scene.core.v1.json` (upstreamable) and `moorhen-scene.ccp4i2.v1.json`; CI drift
   test. *(This increment.)*
2. **Additive layers:** honoured-geometry fields (with units), `hints` block (lighting +
   per-rep material + effects), profiles (permissive default).
3. **Replace** the hand validator; update lifter/resolver/prompt in lockstep (§6).
4. **Generated artefacts:** the compact **LLM authoring brief** (§7a) + token measurement;
   constrained-decoding (JSON→YAML) path.
5. **Fixture corpus + doc-lint** (Tier 1), then **Playwright-on-state** (Tier 2).

## 9. Deferred / open

- **Per-representation material** — deferred until Moorhen gains a material concept (§4b).
- Honoured-geometry field set — **confirmed** against Moorhen `m2tParameters` and shipped
  (bond/ball/probe radii, ribbon widths, vdwScale). Resolver mapping lands in step 3.
- `hints` field set — **confirmed** against Moorhen `glRefSlice` (lighting) and the SSAO/
  edge-detect/outline/shadow/depth-blur/perspective setters (effects); shipped.
- `$id` URL scheme and where the published `moorhen-scene.{core,ccp4i2}.v1.json` are hosted.
- Resolver/lifter wiring for the new fields (step 3) — lifter should emit-only-non-default.
