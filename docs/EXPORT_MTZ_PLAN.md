# Export MTZ — design & wiring plan

Status: **design, not yet implemented** (2026-07-10, `django` branch).

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

**2b. Track the unsplit MTZ as a `CMtzDataFile` output.** *(recommended — makes
export robust and purge-safe)*

Today the natural unsplit files our locators serve are **untracked** and, worse,
**on the default purge kill-list**:

- `refmac.py:238` reads `hklout.mtz` only as the *source* to split; the def.xml
  `outputData` declares only the split mini-MTZs (FPHIOUT/DIFFPHIOUT/…). So
  `hklout.mtz` is never gleaned.
- `servalcat.py` likewise splits `refined.mtz` / `refined_diffmap.mtz` and
  declares only the mini-MTZs.
- `core/CPurgeProject.py SEARCHLIST` puts **`hklin.mtz` and `hklout.mtz` in
  category 1** ("scratch, always safe to delete") — purged in **every** context
  including the gentlest `script_finish` `[1,2]`. So the instant a job/project is
  purged, the locator's target vanishes and Export MTZ silently breaks.

Fix per task (**refmac, servalcat, aimless, dimple**):

1. **def.xml** — add an `outputData` `CMtzDataFile` (e.g. `COMPLETE_MTZ` /
   `HKLOUT_UNSPLIT`).
2. **wrapper** — in `processOutputFiles`, set + annotate it to the unsplit path
   *before* splitting, so the gleaner tracks it as a `File` row.
3. **purge** — ensure its filename is **not** matched by the category-1
   SEARCHLIST globs. `hklout.mtz` currently IS — so either write the tracked copy
   under a non-scratch name (e.g. `complete.mtz`) or add an explicit keep rule.
   Reconcile the `hklout.mtz`/`hklin.mtz` category-1 entries so we never declare
   an output and then purge it.

Once tracked, the **generic fallback finds it** and the plugin **locator becomes
a fallback for legacy (pre-tracking) jobs**. For aimless/dimple this depends on
2a (there is no single unsplit file until the gemmi join produces one).

## Open questions

- Multiple output MTZs in the generic fallback: pick-list (current) vs. zip-all.
- Naming of the tracked unsplit output + whether to backfill legacy jobs.
