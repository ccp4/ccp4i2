# Moorhen scenes: authoring by an LLM — scoping notes

Status: scoping, not a plan. This records what is true today, what the real
choices are, and a suggested order — so we can decide before building. It
assumes the scene format and resolver described in `MOORHEN_SCENES.md` and the
grammar in `types/moorhen-scene.md`, and it builds on the case made in
`MOORHEN_SCENES_UPSTREAM_PROPOSAL.md` for why a validated, declarative scene is a
good thing for a program — including a language model — to produce.

## The premise, in one paragraph

A scene is data, not code. A model that drafts one cannot invent a method name
or leave the viewer in a broken state; the worst it can do is emit YAML that
fails `parseScene`, which is a usable error to repair against, or a selection
that is slightly wrong, which the clamp-and-log resolver absorbs. So the
interesting questions are not "can a model do this" but the practical ones around
it: how does it name the files to load, how does it learn what is actually inside
those files, and how do we deliver the feature across desktop and deployed
without making either depend on the other.

## 1. File references: what actually works, and where

The grammar offers several ways to name a file. They are not equivalent — they
differ in where they can load from. Measured against the apply-time fetcher
(`handleFetchSceneFile` in `components/moorhen/moorhen-wrapper.tsx`):

| Ref form | Loads from scratch? | Electron | Deployed (browser) |
|----------|--------------------|----------|--------------------|
| `pdb: 1jst` | yes — PDBe proxy | yes | yes |
| `url: https://…` | yes — direct fetch | yes | yes (any URL the browser can reach) |
| `fileId` + `projectId` | yes — ccp4i2 proxy | yes | yes |
| `job` + `param` | **not wired** (type stub only) | no | no |
| `bundle:` (`.scene.zip`) | yes — in-memory asset | yes | yes |
| `cifText` (dictionaries) | yes — inline | yes | yes |
| `path:` (local filesystem) | **not implemented** | no | no |

Two of these are traps worth calling out:

- **`path:` does not load anything.** It is documented as an absolute local
  path, but the fetcher has no branch for it. It only *matches an already-loaded*
  molecule (by comparing the path against the open molecule's `uniqueId` in
  `matchOneFile`). A scene that names a fresh local path you have not already
  opened silently fails to load it — including in Electron. Do not advertise
  `path:` to an authoring model.

- **`job` + `param` is a stub.** The fields exist on `SceneFileRef` ("job
  number" + "job parameter name, e.g. XYZOUT") but neither fetcher resolves them.
  This is the one ccp4i2-specific gap worth closing — see section 4.

### Local files are a separate feature, not a ref form

If we want a browser-deployed Moorhen to open a file off the user's disk at all,
that is its own feature, decoupled from the scene grammar:

- **Browser:** the File System Access API (`showOpenFilePicker`) exists in
  Chromium (Chrome/Edge/Opera), not Firefox/Safari. It is picker-and-handle
  based and user-gesture gated — it never resolves an arbitrary path string, so
  it does **not** make `path:` work. What it enables is an "Open local file…"
  button: the user picks a file, we load it as a molecule, and the scene then
  references that molecule by name like any other loaded structure.
- **Electron:** the renderer can read a path directly over IPC (`fs.readFile` in
  the main process). This is the clean route on desktop and could back a real
  `path:` implementation there if we decide we want one.

Either way the conclusion for *scene* refs is the same: keep the grammar's
portable forms (`pdb`/`url`/`bundle`, plus the ccp4i2-internal ones) and treat
local-file loading as a picker-driven action the user performs, after which the
scene refers to the result by name.

## 2. The missing ingredient: a contents summary

A model can write the structure of a scene from the grammar alone. It cannot
guess that the ligand is ATP on chain A at residue 600, or which chains form the
dimer. Those facts live in the file, and gemmi reads them in milliseconds. So an
authoring prompt needs three inputs, not two:

1. the grammar (`types/moorhen-scene.md`),
2. a **contents summary** of the loaded structure — chains, ligand three-letter
   codes *with positions*, residue ranges,
3. the natural-language request.

We already have what we need to produce (2) on the client: the same gemmi WASM
module the slab walk uses (`window.CCP4Module`). It can introspect the
*actually-loaded* molecule regardless of where it came from, so the summary is
always ground truth, never a guess from the request.

## 3. The delivery gradient: don't gate anything on API keys

There are several ways the model call can happen. They share one pluggable
draft-service interface; only the credential/compute source differs. Build the
floor first — it needs no key and serves everyone.

| Tier | Requires | Experience |
|------|----------|------------|
| Managed deployment (e.g. Azure OpenAI) | admin config | one-click in-app drafting |
| Platform on-device (Apple / Windows) | supported hardware + OS | one-click in-app drafting, no key, no cost |
| BYOK desktop | user's own API key | one-click in-app drafting |
| **Manual (floor)** | **nothing** | "Copy prompt" → paste into any chatbot → paste YAML back |

Notes on each:

- **Manual floor.** A "Copy prompt" button assembles grammar + contents summary +
  request onto the clipboard. The user pastes it into whatever model they already
  have and pastes the YAML back. This needs no key, no backend, and no provider
  choice — and it stays relevant even for capable users, because someone with a
  flat chatbot subscription would rather not pay metered API rates for a quick
  view. For many people this is the cheaper primary, not a fallback.
- **BYOK desktop.** A motivated subset of users will paste an API key. Electron's
  `safeStorage` (OS keychain) encrypts it at rest; the main process holds it and
  makes the call. Reasonable to offer, wrong to depend on — capability is not the
  barrier, willingness and metered cost are.
- **Managed.** Where an LLM endpoint is configured (Azure OpenAI is the concrete
  case in our deployment), the key lives server-side and the user does nothing.
  Gate the feature on config presence — the same runtime-config pattern we use
  for auth — so un-configured instances degrade to the manual floor.
- **Platform on-device.** Apple's Foundation Models framework (on-device model
  behind Apple Intelligence, Apple silicon + 2025 OS) and Windows' AI APIs /
  Phi Silica (Copilot+ PCs) expose an on-device LLM to developers — key-free,
  private, no metered cost. They attack the BYOK barrier directly by removing
  both the key and the bill. Three caveats keep them an accelerant for a subset,
  not a floor: (a) **access needs a small native shim, not a generic bridge** —
  see below; (b) the on-device models are small (~3B class), so reliability for
  this task is empirical — though their constrained / guided-generation features
  are a good fit for forcing schema-valid output and lighten the repair loop;
  (c) they are hardware/OS gated (Intel Macs, older Windows, all Linux get
  nothing). Over time these may make BYOK largely unnecessary on Mac/Windows. A
  local-runtime model (Ollama / llama.cpp) fills the same key-free slot
  cross-platform for unsupported hardware and Linux.

  **On the bridge, and where it lives.** PyObjC does *not* reach Foundation
  Models: it reflects the Objective-C runtime, and Foundation Models is
  Swift-only — its `@Generable`/`@Guide` macros, guided generation and
  `LanguageModelSession` have no Objective-C surface. So you write a thin **Swift
  shim** that calls the framework and exposes a C ABI (`@_cdecl`) or JSON over
  stdio, then drive it from Python (`ctypes`/`subprocess`) or Node alike. Windows
  is the parallel case: Phi Silica is WinRT, reached via `pywinrt` or a small
  C#/WinRT shim. Either way it is a per-platform shim, and it only runs where the
  framework does — **on the user's Mac/Windows machine, never on the Linux
  deployment**. The good news is the desktop app already runs a local
  `ccp4-python` Django server (`main/ccp4i2-django-server.ts`), which is the
  natural host: on Mac desktop the same `/api/scenes/draft/` endpoint is backed
  by the Swift shim (a `subprocess`/`ctypes` call) instead of a cloud call —
  same contract, platform-specific implementation, and no Electron-renderer
  native-addon work. So "where the model lives" — Azure OpenAI on the Linux
  deployment, a Swift shim on Mac desktop, a WinRT shim on Windows desktop — is a
  per-deployment backend behind one endpoint; each platform needs its own small
  shim, none of them is PyObjC.

### v2 architecture (managed / BYOK)

- A **thin server endpoint** (`/api/scenes/draft/`) holds the endpoint+key
  server-side, builds the system prompt, calls the model, returns YAML. The key
  never reaches the browser.
- The **contents summary is computed client-side** (gemmi WASM) and posted with
  the request, so the server stays a pure keyed proxy.
- The **validate→repair loop closes in the client**, which already owns
  `parseScene` — it validates the returned YAML and, on error, posts the error
  back for a one-shot repair turn. Keeping the loop client-side avoids porting
  the TypeScript validator to Python.
- Keep the service **provider-agnostic** behind a small interface keyed on
  config — Azure OpenAI where that is what an instance has, a platform on-device
  model (Apple Foundation Models, Windows Phi Silica) or a local runtime (Ollama)
  on desktop, another endpoint elsewhere. Each is one implementation of the same
  interface (NL request + contents summary → YAML); the cloud ones are a web
  call, the platform ones are a per-platform Swift/WinRT shim hosted by the local
  desktop `ccp4-python` server (not PyObjC — see the delivery-gradient section),
  but the prompt scaffold, contents summary, and validate→repair loop are
  identical across all of them. (For new cloud build-out one would normally reach for the latest
  Claude models; the point is that the architecture does not care, so don't
  hardcode a provider.)

A useful property throughout: because the resolver is deterministic and the
scene is validated before apply, a bad generation is caught rather than executed.
The blast radius is a failed parse or a clamped selection, not an arbitrary API
call.

## 4. The ccp4i2 dialect: reference files by job and role

Outside ccp4i2, `pdb`/`url` are the natural forms for a model to emit. Inside a
ccp4i2 project they are the *wrong* ones. The user does not think "the file at
`/api/proxy/ccp4i2/files/1207/download/`" — they think "the output coordinates
of job 32." The grammar already has the form for that thought: `job` + `param`
(and `fileId` + `projectId`). And the server already holds the mapping — `File`
carries `job` (FK) and `job_param_name` (`"XYZOUT"`), plus `type`/`sub_type`/
`annotation` — so resolving "job 32 / XYZOUT" is a one-line ORM query.

This is the contents-summary idea one level up. Instead of "what is inside this
structure," the app injects "what files exist here, by job and role":

| Summary the app injects | Lets the model reference… | It emits |
|-------------------------|---------------------------|----------|
| gemmi structure contents | chains / ligands / ranges | `selection:` CIDs |
| **project / job manifest** | **files by job and role** | **`job: N, param: XYZOUT`** |

References by job and role are better than `url`/`path` here: they are stable
across re-runs, human-meaningful, and they match the mental model the user
already has. A wrong ref (`job: 99, param: NOPE`) simply fails the lookup and
feeds the repair loop.

Two surfaces follow:

- **Job view (current job implicit).** The page already knows its job, so the
  manifest is just this job's outputs and the model only names the role — "the
  output coordinates" → `param: XYZOUT`, job supplied by context.
- **Project-scoped Moorhen page (cross-job).** The powerful case: pull job 32's
  refined model and job 35's difference map into one view. There is a precedent —
  the campaign-scoped Moorhen page at
  `app/ccp4i2/(authed)/moorhen-page/campaign/[id]/` — so a general project-scoped
  surface extends an existing pattern rather than inventing one.

The dialect is therefore **context-aware**: the app advertises the right
reference forms and the right "here's what you can reference" summary for the
surface — one structure's internals, one job's outputs, or a whole project's
jobs.

## 5. Decisions to settle

1. **Authoring dialect per context.** Confirm the constraint sets: `pdb`/`url`
   for generic Moorhen; `job`+`param`/`fileId` first inside a ccp4i2 project, with
   `pdb`/`url` for genuinely external structures. Forbid `path:` in either.
2. **Where the contents summary is computed.** Client-side gemmi WASM (proposed)
   versus a server endpoint. Client-side reuses code we have and always reflects
   the loaded molecule.
3. **Manual floor first?** Whether "Copy prompt" + the summary ships before any
   in-app model call. Proposed: yes — it is universal and provider-free.
4. **Validator location for the repair loop.** Keep it client-side (proposed) or
   port a minimal validator to Python.
5. **Provider abstraction.** Shape of the pluggable draft-service interface and
   how config (and platform detection) selects a provider — cloud endpoint,
   platform on-device model, local runtime, or none. Whether to invest in the
   per-platform native bridges (Apple Foundation Models / Windows Phi Silica)
   now or treat them as later implementations of the same interface.
6. **Project-scoped page scope.** Whether to generalise the campaign page into a
   project-scoped Moorhen surface now, or defer until the job/param refs land.

## 6. Suggested order

The work that pays off in every tier is the summary plumbing, so do it first.

1. **Contents summary (client, gemmi WASM)** + a **"Copy prompt" button** in the
   Scenes tab. This is the whole manual floor and it is provider-free.
2. **Wire `job` + `param`** in the coordinate and dictionary fetchers (lookup →
   download URL), plus a **project/job manifest** the Scenes tab can fold into the
   prompt. This unlocks the ccp4i2 dialect.
3. **The draft endpoint and in-app box** (managed + BYOK), reusing the same
   prompt scaffold and the same client-side validate→repair.
4. **Project-scoped Moorhen page**, if cross-job authoring proves worth a
   dedicated surface.

## 7. Concrete touch points and known gaps

- `components/moorhen/moorhen-wrapper.tsx` — `handleFetchSceneFile` /
  `handleFetchSceneDictionary`: add the `job`+`param` branch (gap today).
- `lib/moorhen-scene.ts` — validator: the context-aware dialect constraints
  (which ref forms are permitted) if we choose to enforce them at parse time.
- A server endpoint to resolve `job`+`param` → a `File` (or a manifest the client
  resolves against), backed by `File.job` / `File.job_param_name`.
- A server **manifest** endpoint: jobs + their output params + file types for a
  project (and the single-job slice for the job view).
- A `/api/scenes/draft/` endpoint for v2, holding the model credential
  server-side; provider-agnostic behind config.
- Scenes tab UI: "Copy prompt" (floor) and, later, a "Draft scene" box.
- Out of scope here: a real `path:` implementation (Electron IPC) and an "Open
  local file…" picker (File System Access) — both useful, both separate from
  scene refs.
