# Moorhen Scene authoring UI — natural-language & capability-dependent design

Status: **design, pre-implementation.** Records how the scene panel offers
natural-language authoring, and how that UI adapts to what's actually available
in the current deployment (integrated LLM endpoint vs copy-paste vs neither).
Companion docs: `MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` (the format + §7a tiered
prompting + §12 the `nlp_scene` endpoint contract), `types/moorhen-scene.md`
(grammar). Materia-side notes live in `../../../materia/notes-for-ccp4i2/`.

## 0. Goal & principle

Authoring a scene ranges from hand-writing YAML to "describe it in a sentence and
get a live view". What's *available* differs per deployment: Materia-web has an
integrated Azure OpenAI endpoint; desktop Electron has no server-side LLM; any
context can reach an external chatbot by hand. The UI must **detect capability at
runtime and present the highest-value affordance that works here**, while never
removing the floor.

Guiding principle — **progressive disclosure, never hide the floor:**

> The YAML editor (with schema validation) is always present. Everything else is
> an *enhancement layered on top*, shown only when its prerequisite exists. A
> missing LLM degrades the experience; it never blocks authoring.

## 1. Three capability tiers

| Tier | Prerequisite | Affordance | Always on? |
|---|---|---|---|
| **1 · Hand-edit** | nothing | YAML editor + JSON-Schema markers/autocomplete + Apply | ✅ the floor |
| **2 · Copy-paste** | a clipboard | "Copy as prompt" → paste into any chatbot → paste result back → Apply | ✅ universal, provider-free |
| **3 · Integrated** | a live NL→scene endpoint | in-app **Generate** → endpoint → fills editor | only where detected |

Tiers 1 and 2 already exist in `components/moorhen/moorhen-scenes-panel.tsx` (the
editor, and the "Generate a scene-authoring prompt" modal). Tier 3 is the new
work — and it is a **thin cap on infrastructure that already exists**, not a
parallel UI (see §3).

Tier 2 is deliberately kept even when tier 3 is live: it is the escape hatch to a
**better/other model** than the deployment's, and the fallback when the endpoint
is rate-limited or down (§5). Copy-paste is not a legacy path; it's the
provider-agnostic ceiling.

## 2. Capability detection — a provider abstraction, not a Materia dependency

The frontend must not hard-depend on Materia. It asks a single hook "is a
scene-NL provider available here, and how do I call it?", with pluggable
providers:

```ts
useSceneNlCapability(): {
  available: boolean;
  provider: "none" | "materia" | /* future */ "desktop-byo" | "on-device";
  responseFormat?: "strict" | "raw";   // hint for ingest (see §4)
  model?: string;                       // for display
  quota?: { limit?: number; remaining?: number };
  generate(request: string, grounding?: string): Promise<string>; // raw model text
}
```

- **`materia` provider** — probes `GET /api/proxy/compounds/nlp/status` (see §6),
  and on Generate POSTs `/api/proxy/compounds/nlp/scene`.
- **`desktop-byo` / `on-device`** *(future, not built)* — a user-supplied API key
  called from the Electron main process, or a bundled small model (the
  `MOORHEN_SCENES_SCHEMA_V1_DESIGN.md` §7a on-device aspiration). Adding one is a
  new provider behind this same hook — a drop-in, not a rework.

The UI consumes only `available` + `generate()`; it never branches on which
provider answered.

### Desktop degrades correctly for free

Desktop Electron runs the local ccp4i2 Django, which does **not** serve the
compounds app — so the `materia` provider's status probe 404s, `available` is
`false`, and tier 3 is simply absent. Desktop users get hand-edit + copy-paste,
which is the *right* outcome (copy-paste into their preferred external model needs
nothing). **We do not build a server endpoint into desktop** to make degradation
work; absence is the signal. In-app generation on desktop, if ever wanted, comes
from a `desktop-byo`/`on-device` provider — our side, no server.

## 3. One input, context-dependent action

The existing modal ("Describe the view you want…", `promptRequest` state) is
already the NL input surface. Tier 3 does **not** add a second box — it changes
the modal's **primary action** based on `useSceneNlCapability().available`:

| Capability | Primary button | Secondary |
|---|---|---|
| tier 3 available | **Generate** (calls `generate()`, fills the editor) | Copy as prompt |
| tier 3 absent | **Copy as prompt** (today's behaviour) | — |

This mirrors the panel's existing capability-gating idiom: `onBuildAuthoringPrompt`
is an *optional* prop, and the prompt button is hidden when it's absent. Tier 3
adds a sibling optional prop (e.g. `onGenerateScene?`), injected by the wrapper
only where the capability hook reports available. The panel stays dumb about
providers.

Everything downstream of "we now have scene text" is **shared** across tiers 2
and 3: the text lands in the editor → JSON-Schema markers show validity → the
user Applies. One review surface, one apply path, one renderer.

## 4. The Generate flow (review-then-apply)

```
NL request ──▶ generate() ──▶ raw model text ──▶ ingest ──▶ editor (markers) ──▶ [user] Apply ──▶ render
                (provider)      YAML or JSON                  validate
```

- **No auto-apply** (decision, §7). Generation is imperfect and strict mode may be
  off, so the result is *previewed* in the editor with validation markers; the
  human Applies. This reuses the existing trust surface — tier 3 introduces no new
  one.
- **Ingest handles both shapes.** `parseScene` already strips a ` ```yaml `/
  ` ```json ` fence and runs `YAML.parse` (JSON ⊂ YAML), so raw-YAML (spike) and
  strict-JSON (production) both parse. **One addition for strict:** OpenAI strict
  mode emits explicit `null`s for absent optionals, and Zod treats *optional* ≠
  *null* — so **strip nulls before `parseScene`**. This is the single genuinely
  new bit of client logic tier 3 needs. `responseFormat: "strict"` from the
  capability hint tells us when it's needed, but stripping nulls unconditionally
  is safe (a hand-authored scene has none).
- **Grounding** is the same project/manifest/contents block `buildAuthoringPrompt`
  already assembles (the wrapper's `handleBuildAuthoringPrompt`, incl. PDB digest
  folding); for tier 3 it becomes the endpoint's `grounding` field rather than
  inline prompt text.

## 5. Fallback & error states

Capability is not a one-shot boolean at load — live failures must degrade in
place. The `nlp_scene` POST error kinds map to UI:

| Kind | UI response |
|---|---|
| `rate_limited` (daily cap) | Inline notice + **fall back to "Copy as prompt"** ("quota reached — paste into your own model instead") |
| `disabled` | Same fallback; note the feature is off |
| `llm_error` (upstream 502) | Retry offered, then copy-paste fallback |
| `bad_request` | Surface the message (missing request) |

So `nlp/status` drives the *default* UI; the POST error kinds drive *live*
fallback if state changes between the status check and the call. Copy-paste is the
common safety net for all of them — another reason it stays visible.

## 6. The `nlp/status` capability endpoint (implemented, Materia)

Detection uses a cheap, read-only, no-quota probe. **Implemented by Materia**
(commit `6095575`; ask in
`../../../materia/notes-for-ccp4i2/ccp4i2-ask-nlp-status-endpoint.md`):

`GET /api/proxy/compounds/nlp/status` (no trailing slash; `IsAuthenticated`; no
model call, no cap spend), returning per-feature capability:

```jsonc
{ "scene": { "enabled": true, "responseFormat": "raw",   // "raw" live today
             "model": "gpt-4o", "dailyLimit": 500, "dailyRemaining": 500 },
  "query": { "enabled": true } }
```

The two semantics that drive the UI are honoured:

- **The status endpoint is NOT gated by the feature flag.** So where the compounds
  app exists it always answers `200` and states `enabled` explicitly:
  `enabled:true` (on) or `enabled:false` (present-but-off, → "disabled on this
  deployment"). A **`404`** therefore means only "no NLP app here" (desktop
  Electron), never "disabled". We do not lean on 404 to detect the off-state — it
  is purely the absent-app signal, which is the honest truth on desktop. This is
  co-location: the endpoint reporting the capability lives in the same app that
  provides it, so absent app ⇒ 404, by construction — not a new backend
  obligation, and deliberately *not* added to ccp4i2 core (which must not know
  about the optional compounds app).
- **`responseFormat` reflects what `generate_scene` will actually use** — `strict`
  only when `MOORHEN_SCENE_STRICT` is on *and* the schema loads, else `raw` (the
  runtime fallback). Lets the client anticipate JSON (strip nulls) vs YAML;
  stripping nulls unconditionally is safe regardless.

Client side: SWR-cached; `404`/unreachable → `available:false` (quiet, expected);
transient `5xx`/timeout → retry, still safe-fallback to copy-paste. Only
`scene.enabled` is load-bearing; the rest is UX polish.

## 7. Decisions

**Locked:**
- **Review-then-apply**, no auto-apply — generated scenes preview in the editor.
- **Copy-paste stays** as a secondary action even when integrated is live (better
  model / offline / quota fallback).
- **Strict is proven** end-to-end on DDU (raw + strict both green), so the client
  targets **strict JSON as the production shape** → null-strip before `parseScene`.
- **Capability is a provider abstraction**; desktop degrades via absence, no
  desktop server endpoint required.

**Deferred:**
- `desktop-byo` (BYO-API-key from Electron) and `on-device` providers.
- Showing remaining-quota / model name in the button (needs the optional
  `nlp/status` fields; nice-to-have).
- Streaming the generation into the editor (v1 waits for the full document).

## 8. Implementation sketch (sequencing)

1. **`useSceneNlCapability()` hook** with the `materia` provider only: SWR probe of
   `nlp/status`; `generate()` = POST `nlp/scene` (no trailing slash) + fence/null
   strip.
2. **Wrapper wiring** — inject `onGenerateScene?` into the panel only when
   `available`; reuse the existing grounding assembly for the `grounding` field.
3. **Panel** — swap the modal's primary action on `available`; keep "Copy as
   prompt" as secondary; route both results into the editor (no auto-apply).
4. **Fallback** — map POST error kinds to in-place copy-paste (§5).
5. **Later** — additional providers; quota/model display; streaming.

§6 (`nlp/status`) is **live** (Materia commit `6095575`), so the hook can consume a
real capability signal today rather than defaulting to unavailable — tier 3 lights
up wherever the endpoint answers `enabled:true`.
