# Export MTZ — design & wiring plan

Status: **Phase 1, 2a, and 2b (a1/a2) IMPLEMENTED** (2026-07-10, `django` branch).
COMPLETE_MTZ is a tracked output on refmac, servalcat, prosmart_refmac,
servalcat_pipe, aimless_pipe and i2Dimple; the toolbar Export MTZ button
surfaces and downloads it. See the Phase 2b section for the per-pipeline table.

The job-panel toolbar (`client/renderer/components/tool-bar.tsx`) has an **Export
MTZ** button that is currently a no-op (`onClick: () => {}`), gated by a
misleading blanket `job.status === 6` check. This doc captures the truth we
discovered about how MTZ export *should* work, and the plan to wire it.

## The pre-existing per-plugin contract

Some plugin scripts already declare how to export their MTZ, via two
**module-level functions** (a legacy Qt-era contract, never re-wired after the
Django move):

```python
def exportJobFileMenu(jobId=None):
    # -> [[mode, label, mimeType], ...]
    return [['complete_mtz', 'MTZ file', 'application/CCP4-mtz']]

def exportJobFile(jobId=None, mode=None, fileInfo={}):
    # -> absolute path to the file to serve
    ...
```

**6 tasks implement it** (`grep -rl "def exportJobFileMenu"`):
`aimless_pipe`, `import_merged`, `servalcat_pipe`, `prosmart_refmac`, `parrot`,
`AlternativeImportXIA2`.

They fall into two behaviours:

| Kind | Tasks | What `exportJobFile` returns |
|------|-------|------------------------------|
| **Locator** | `servalcat_pipe`, `prosmart_refmac` | The **natural, never-split** output — burrows into child jobs and returns `refined.mtz` / `hklout.mtz` as-is. |
| **Reconstructor** | `aimless_pipe`, `import_merged` | A **catenation** — there is no single natural file, so it recombines subjob outputs (ctruncate intensities/amplitudes + FreeR column) via `runCad` into a fresh `exportMtz.mtz`, with intensity-vs-amplitude logic and FreeR-as-base-dataset handling. |

`parrot` / `AlternativeImportXIA2` are further locator-style variants.

### Design principle: natural unsplit file first

The priority within any job is:

1. **Serve the natural unsplit output** where one genuinely exists (locator
   tasks, and the generic fallback below). Lossless, no recombination.
2. **Reconstruct** only where the data was genuinely produced across subjobs and
   there is no single file (import/scaling pipelines). Reconstruction is the
   *correct* answer there, not a compromise.
3. **Generic fallback** for any task that declares nothing: offer the job's
   tracked output-MTZ `File`(s) directly (still natural, just not task-curated).
   Button **hidden** only when there is genuinely no MTZ to export.

## What already exists (reuse, don't rebuild)

- **API routes** in `server/ccp4i2/api/JobViewSet.py`:
  - `export_job_file_menu` (`GET /api/jobs/{id}/export_job_file_menu/`) — **stub**:
    hardcoded to log "not available" and return `[]`.  *(the reason the button
    is dead)*
  - `export_job_file` (`GET /api/jobs/{id}/export_job_file/?mode=…`) — imports
    `export_job_file` from `lib/utils/jobs/export.py`, **which does not define
    it** → would `ImportError` at call time.
- **Django compat layer** `server/ccp4i2/db/ccp4i2_django_dbapi.py`
  (`class CCP4i2DjangoDbApi`) already exposes the exact methods the legacy
  functions call, backed by the ORM:
  - `getChildJobs(jobId, details=True) -> [(number, uuid, task_name), ...]`
    (same tuple shape the plugins unpack)
  - `jobDirectory(jobId=…) -> path`
  - `getJobFilesInfo(jobId=…, jobParamName=…)`
- **`PROJECTSMANAGER().db()`** (`core/CCP4ProjectsManager.py:120`) returns a
  `CCP4i2DjangoDbApi()` instance. So the legacy `exportJobFile` functions —
  which call `CCP4Modules.PROJECTSMANAGER().db().getChildJobs(...)` — resolve
  against Django **essentially unmodified**. The "rewrite" collapses to "verify
  they run + tidy"; do **not** touch the crystallographic reconstruction logic.
- Task registry (`core/tasks.py: get_plugin_class`, plugin module import) to
  reach the module-level functions.
- Streaming-download pattern to mirror: the existing `export_job` ZIP action.

## The one real gap: `runCad` is gone

The reconstructor tasks call `m.runCad(...)` on a `CMtzDataFile`. **`runCad` has
no definition anywhere in the Django branch** — so `aimless_pipe` /
`import_merged` export would crash today. Nearest replacements:

| Utility | Direction | gemmi? | Notes |
|---------|-----------|--------|-------|
| `lib/utils/formats/mtz_utils.extract_columns` | split | ✅ gemmi | replacement for `splitHklout`/`sftools` |
| `CCP4Utils.merge_mtz_files_cad` | join/merge | ❌ | **still shells out to the `cad` binary** (`shutil.which('cad')`) |
| `CCP4PluginScript.makeHklin` | join | ✅ gemmi | targets mini-MTZ *program input*, not general column-preserving CAD catenation |

So there is **no gemmi CAD-equivalent for the join direction yet**. This matters
because a CCP4-binary dependency conflicts with the slim-server architecture
(server is CCP4-free; jobs run in ccp4-python — see
`memory: project_slim_env_separation_guard`).

**Consequence for phasing:** locators + generic fallback serve a real file
directly and have **no binary dependency** → ship first. Reconstructors need
either (a) a new `merge_mtz_files_gemmi` (preferred), or (b) to run export in the
ccp4-python job environment, not the CCP4-free server.

## Plan

### Phase 1 — locators + generic fallback (no binary dependency)

1. **`lib/utils/jobs/export.py`** — add the missing glue:
   - `export_job_file_menu(job)`: import the plugin *module* via the registry;
     if it defines `exportJobFileMenu`, call it. **Else** synthesise a menu from
     the job's tracked output-MTZ `File`s (generic fallback).
   - `export_job_file(pk, mode)`: for a plugin mode, call the module's
     `exportJobFile(jobId, mode)` (runs against `CCP4i2DjangoDbApi`); for a
     fallback mode, resolve `File.path` directly. Stream via `FileResponse`.
     Guard against reconstructor tasks until Phase 2 (detect the dead `runCad`
     path / catch and report a clear "reconstruction not yet available" error).
2. **`JobViewSet.export_job_file_menu`** — un-stub: call the util instead of
   returning `[]`.
3. **`tool-bar.tsx`** — fetch `export_job_file_menu` for the job; **hide the
   button when the menu is empty** (replaces the blanket `status===6` gate — the
   plugin menu already checks file existence, which implies finished). Single
   mode → direct download; multiple → small pick menu. Wire `onClick` to
   `export_job_file?mode=…`.

Net after Phase 1: Export MTZ correct for `servalcat_pipe` / `prosmart_refmac`
(the screenshot cases) + any task with a tracked output MTZ; truthfully
**absent** where there is nothing to export.

### Phase 1 — IMPLEMENTED (2026-07-10)

Shipped:

- `core/tasks.py`: new `get_plugin_module(task_name)` accessor (imports the
  module holding the plugin class, to reach module-level export functions).
- `lib/utils/jobs/export.py`: `export_job_file_menu(job)` and
  `export_job_file(pk, mode)` + helpers. Menu prefers the plugin contract,
  **validates each declared mode by resolving it to an on-disk file** (the
  plugin `exportJobFileMenu` is optimistic and advertises a mode even when the
  file is absent), and falls back to the job's tracked output MTZ `File`s
  (`_FALLBACK_PREFIX = "file:"`). Reconstructor tasks return HTTP 501.
- `api/JobViewSet.py`: un-stubbed `export_job_file_menu`; cleaned debug prints
  from `export_job_file`.
- `servalcat_pipe.exportJobFile`: made robust — searches **all** `servalcat`
  subjobs (not the fragile `childJobs[-1]`/`[-2]` positional assumption) and
  covers **SPA mode** (`refined_diffmap.mtz`) as well as xtal (`refined.mtz`).
- `tool-bar.tsx`: fetches the menu; **hides the button when empty** (replaces
  the misleading blanket `status===6` gate); single mode → download, multiple →
  pick menu.

Verified end-to-end against a real Django DB (in-memory test runner): servalcat
xtal + SPA locators, refmac generic fallback (FileResponse 200), parrot
(menu-but-no-`exportJobFile`) → empty, aimless reconstructor → 501.

### Phase 2 — reconstructors + tracked unsplit MTZ

**2a. Reconstructors — IMPLEMENTED (2026-07-10).** Rather than a new util, the
branch already had a gemmi join: `CCP4Utils.merge_mtz_files` (pure gemmi despite
its "using CAD" docstring; distinct from the binary-shelling
`merge_mtz_files_cad`). Added a thin wrapper `combine_mtz_files(sources,
output_path)` in `lib/utils/jobs/export.py` that copies all data columns from
each source (first wins on clash). Repointed `aimless_pipe` /`import_merged`
`exportJobFile` from the dead `m.runCad(...)` to `combine_mtz_files(...)`
(observed data first, FreeR second) and removed the `_RECONSTRUCTOR_TASKS` 501
guard. No `cad` binary — satisfies the slim-server guard. Verified: a real
gemmi merge of F/SIGF + FreeR mini-MTZs produces the combined file.

**2b. Track a `COMPLETE_MTZ` output on the pipelines — the finalized 3-layer
design.** *(makes Export MTZ predictable regardless of the cleanup flag)*

The naive "keep the child subjob's `hklout.mtz`" approach is a **dead end** —
proven by a real prosmart_refmac run (`REFMAC_CLEANUP=True`, gamma data): the
refmac `hklout.mtz` was still deleted. Why:

- **Untracked.** `refmac.def.xml` declares **no** HKLOUT output; `refmac.py:238`
  uses `hklout.mtz` only as the split *source*. `servalcat.py` likewise. So the
  file is never gleaned → cannot be protected as a modelled output.
- **Child-dir layout defeats the keep pattern.** The pipeline calls
  `purgeJob(firstRefmac.jobId, context="extended_intermediate")` on the **child**
  refmac job, whose dir holds `hklout.mtz` at top level. prosmart's keep pattern
  is `refmac%*/hklout.mtz` — written as if purging the *parent* — so it never
  matches; the default `['hklout.mtz', 1]` rule deletes it.
- Purge context `extended_intermediate = [1,2,5,7]` includes the categories the
  file falls under.

So the child file is fundamentally fragile. The fix is to **model it as a real
output — at BOTH the wrapper and the pipeline** (each serves a different role):

**Layer a1/a2 — IMPLEMENTED (2026-07-10).**

| Pipeline | PR | How COMPLETE_MTZ is produced |
|----------|-----|------------------------------|
| refmac + prosmart_refmac | #251 (merged) | a1 wrapper: setFullPath to `hklout.mtz`. a2 pipeline: inherits def.xml via `<file>` include; existing `dataOrder()` copy-up in `finishUp` copies it up + annotate. |
| servalcat + servalcat_pipe | #252 (merged) | a1 wrapper: setFullPath to `refined.mtz`/`refined_diffmap.mtz`. a2 pipeline: inherits def.xml; `dataOrder()` copy-up + annotate in `_applyAnnotations`. |
| aimless_pipe | #254 | Reconstructor (pipeline-only): `process_finish` recombines `HKLOUT[0]` + FreeR via `combine_mtz_files`; tracked outputData `CMtzDataFile`. |
| i2Dimple | #254 | **Locator, single wrapper** (plan correction — dimple is NOT a reconstructor): it splits mini-MTZs from `final.mtz` in its own work dir, so `final.mtz` IS the unsplit file. `processOutputFiles` setFullPath to `final.mtz` (no copy). Tracked under its native name in the job's own dir. |

Also fixed the frontend bug that hid the button (#253): `tool-bar.tsx` fetched
`export_job_file_menu` via `get_endpoint`, whose fetcher does not unwrap the
`{success, data}` envelope `api_success()` adds — so `exportMenuData.result` was
`undefined` and the button was filtered out. Now reads `.data?.result ?? .result`.

Verification: TDD i2run tests for all four pipelines pass end-to-end; refmac &
servalcat assert survival past `REFMAC_CLEANUP=True`. The full refmac+servalcat
i2run files are green (8 passed, 2 pre-existing clipper skips).

Design notes below retained for context.

**Layer a1 — model COMPLETE_MTZ on the WRAPPER (refmac, servalcat).** This is the
primary, correct fix: declare the unsplit file as an output where it is
*produced*.
  1. **wrapper def.xml** (`refmac.def.xml`, `servalcat.def.xml`) — add
     `outputData CMtzDataFile COMPLETE_MTZ`.
  2. **wrapper `processOutputFiles`** — right where it already loads
     `hklout.mtz` (refmac.py:238) / `refined[_diffmap].mtz` (servalcat),
     `setFullPath` + annotate `COMPLETE_MTZ` to that path (no copy; the file is
     already in the wrapper's own work dir). The gleaner then tracks it as a
     `File` of the wrapper (child) job.
  Why the wrapper: (i) honest data model — declared at point of production;
  (ii) **#247 then protects it during the CHILD purge automatically** (it is now
  one of the child job's modelled outputs), closing the child-purge gap at its
  root with no keep pattern; (iii) every pipeline that uses refmac/servalcat
  (prosmart_refmac, servalcat_pipe, dr_mr_modelbuild, SubstituteLigand, lorestr,
  …) inherits it.

**Layer a2 — expose COMPLETE_MTZ on the PIPELINE (prosmart_refmac,
servalcat_pipe).** The Export MTZ button acts on the *top-level pipeline* job and
the generic fallback queries *that* job's output Files, so the pipeline must also
carry it:
  1. **pipeline def.xml** — add `outputData CMtzDataFile COMPLETE_MTZ`.
  2. **pipeline finish handler** — `shutil.copyfile` the terminal refmac/servalcat
     child's COMPLETE_MTZ up to the pipeline job dir (the existing idiom, cf.
     LIBOUT/PSOUT copies at prosmart_refmac.py:300-303), set + annotate,
     **before** the child `purgeJob`. #247 then protects it in the pipeline dir
     too. Copy (not move) — cleanup-independent.

For **aimless_pipe / dimple** there is no single wrapper producing an unsplit
file (it is reconstructed), so those are pipeline-only (a2), building the file
via `combine_mtz_files` (2a). Layer a1 applies specifically to refmac & servalcat.

**Layer b — tracked outputs are unpurgeable. DONE (PR #247).** `purgeJob` now
blacklists the real paths of a job's modelled OUT `File`s; `_deleteFile` skips
them whatever their category/name. So once COMPLETE_MTZ is a modelled parent
output, the parent purge cannot touch it — no per-task keep rule needed.

**Layer c — category-0 keep now suppresses same-task rules. DONE (PR #248).**
`_buildSearchList` let a task's `['pat',0]` keep coexist with its own
`['pat',7]` delete (the prosmart bug). Fixed so cat-0 means keep unconditionally.
Belt-and-braces; layers a1/a2 already make export safe.

Result: COMPLETE_MTZ is a tracked parent output → generic fallback finds it →
survives cleanup and project export. The plugin **locator becomes a fallback for
legacy (pre-tracking) jobs**. For aimless/dimple the copied file is produced via
`combine_mtz_files` (2a).

**Verification harness:** the i2run tests (`test_refmac`, `test_servalcat`,
`test_aimless`, `test_dimple`) run the real pipelines and the required binaries
are present. TDD loop per pipeline: add COMPLETE_MTZ → run its i2run test with
`REFMAC_CLEANUP`/cleanup on → assert `(job / "COMPLETE_MTZ.mtz").is_file()`
(proves both production AND survival past cleanup). Add the assertion as a
regression guard.

## Open questions

- Multiple output MTZs in the generic fallback: pick-list (current) vs. zip-all.
- Naming of the tracked output (`COMPLETE_MTZ.mtz`) + whether to backfill legacy
  jobs (probably not — locator fallback covers them).
