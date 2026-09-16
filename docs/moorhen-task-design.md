# A recorded Moorhen task

A design note for a `moorhen` task: a job that opens Moorhen on chosen data,
records what was loaded, and harvests what was saved, so that Moorhen model
building leaves the same provenance trail in the database as the Coot 1, Coot
0.9 and CCP4mg tasks do today.

**Baseline:** commit `ff746109c` (branch `django`), surveyed 2026-09-16.
**Status:** design note, revised the same day for the in-browser (Azure)
deployment, then to defer scene capture, then to record the decisions in
[Decisions](#decisions). **PR 1 (server) is implemented on branch
`moorhen-task`**: `wrappers/moorhen/`, `lib/utils/jobs/interactive.py`, the
`JobInteractiveSession` model (migration 0022), the four `interactive_*`
endpoints on `JobViewSet`, the `force_dispatch` hook in
`run_job_context_aware`, the `Task.interactive` flag, and the shared
`cootbridge/harvest.py` filing helpers. PR 2 (client) and PR 3 (Electron)
are not started.

Companion documents:

- [docs/authoring-a-task.md](authoring-a-task.md) — the task-authoring spine.
- [docs/interrupt-and-resume.md](interrupt-and-resume.md) — glean is gated on
  success; this note works inside that rule rather than changing it.
- `server/ccp4i2/cootbridge/` — the Coot bridge whose contracts this note
  reuses (load plan, drop-directory harvest, handshake).

## The problem in one paragraph

Coot 1, Coot 0.9 and CCP4mg are classical tasks: the job's process blocks on
the GUI, the populated input lists decide what is loaded, and models saved
into a drop directory are harvested into `XYZOUT` and gleaned when the GUI
exits. Moorhen has none of that. It is opened from a file or job context
menu into a fresh Electron window, infers what to load from the job's output
files, and the only way back is "push to CCP4i2", which creates a *new*
`coordinate_selector` job. Nothing records that a Moorhen session happened,
what it was opened on, or which saved model came from which session. The
missing verb is *save into the job this window represents*, and for that
there has to be a job.

## What exists today

| | Coot 1 task (`coot1`) | Moorhen preview |
|---|---|---|
| A job exists | yes | no |
| Inputs chosen by | `XYZIN_LIST`, `FPHIIN_LIST`, `DELFPHIIN_LIST`, `DELFPHIINANOM_LIST`, `DICT` | inferred from a file id, a job's output files, or a project |
| Load order | `api_client.INPUT_PARAM_KINDS`: dictionaries, coordinates, maps | `fetchJobFiles`: dictionaries, one coordinate file (mmCIF preferred), maps |
| Runs as | grandchild process of uvicorn; `subprocess.run` is the wait | a `BrowserWindow` in the Electron app, opened by `window.open` |
| Save contract | `COOT_FILE_DROP/output<N>.pdb|cif`, harvested on exit | push to a new `coordinate_selector` job via `upload_file_param` |
| Provenance | `XYZOUT[n]` gleaned, `params.xml` records inputs | the pushed file is an *import* of a new job |
| Launch from CLI | `i2run coot1 --XYZIN_LIST ...` | none |

The server side already mirrors Moorhen's loader: `cootbridge/api_client.py`
declares `LOADABLE_TYPES`, `_MAP_SUBTYPE_KINDS` and `_LOAD_PRIORITY` with a
comment citing `fetchJobFiles`. Moorhen's own reverse direction is also nearly
there: `push-to-ccp4i2-panel.tsx` serialises a molecule, sniffs PDB against
mmCIF and uploads it. What is missing is a job for that upload to land in,
and a way for anything outside the renderer to open the window.

## Design

### The task

`server/ccp4i2/wrappers/moorhen/script/moorhen.def.xml` carries the same
input container as `coot1`, with two differences decided below:
`XYZIN_LIST`, `FPHIIN_LIST`, `DELFPHIIN_LIST`, `DELFPHIINANOM_LIST`, and

- `DICT_LIST` in place of `coot1`'s single `DICT`: a `CList` of
  `CDictDataFile`. Moorhen re-reads every dictionary against every molecule
  (the "everything is called LIG" fix in `fetchJobFiles`), and fragment
  campaigns routinely need several. A `coot1` job still clones into a
  `moorhen` one: the clone path's `adopt_legacy_container` matches typed
  inputs by name, so it needs one rule, `DICT` into `DICT_LIST[0]`. The
  reverse direction (a `moorhen` job cloned to `coot1`) takes the first
  dictionary and warns if there were more. `coot1` should grow the same
  list when it is next touched.
- `MAPIN_LIST`, a `CList` of `CMapDataFile`, for real-space maps and masks
  (cryo-EM half maps, the `molrep_map` outputs). Moorhen already loads
  these through `fetchMapFile` with mask handling. Optional;
  `allowUndefined`.

`outputData` matches `coot1`: `XYZOUT` and `DICTOUT`, both
`COutputFileList`.

**Scene capture is deliberately out of scope.** The liftable scene
description (`client/renderer/MOORHEN_SCENES_SCHEMA_V1_DESIGN.md`) is
expected to be ratified upstream, but that will take time, and a `SCENEOUT`
recorded against a format that may still move would have to be migrated.
Nothing in this design depends on the scene format, so a scene output can be
added later as one more `COutputFileList`, one more drop kind and one more
mime type, without touching the session mechanics.

Registry entry in `core/tasks.py`: title "Moorhen", description "Interactive
model building with Moorhen", same chooser category as `coot1`
("Model building and Graphics" in `task-chooser.tsx`). Because it needs no
CCP4 binary at all, it is `ccp4_free=True`, which is also what lets it run
on the slim server.

### Execution model: a session row, and a plugin that waits only if it must

The Coot tasks work because the plugin process blocks on a child. Moorhen
has no child process, and the job may not even run on the machine showing
the window (see *The in-browser deployment* below). So the session is
separated from the process:

- **The session lives in the database**, not in a process. A small
  `JobInteractiveSession` row (or fields on `Job`) holds `requested_at`,
  `last_heartbeat`, `finished` and `finished_at`. Both the desktop job
  process and an Azure worker have Django and the database; neither is
  guaranteed a fast, shared filesystem for polling marker files.
- **Run does not dispatch the job.** For a task registered as interactive,
  `POST jobs/{id}/run/` validates, saves `params.xml`, computes the load
  plan, sets the status to `RUNNING`, creates the session row, and returns.
  No process exists yet. The client opens the session route.
- **Finish dispatches it.** `interactive_finish` marks the session finished
  and then dispatches the job through the ordinary `run_job_context_aware`,
  so it runs locally on the desktop and through Service Bus on Azure. The
  runner does what it always does: `track_job`, `plugin.process()`, glean on
  success, `FINISHED`.
- **The plugin waits only if the session is not finished.**
  `moorhen.startProcess()` (38 wrappers already override this for
  pure-Python work) checks the session row; if finished, it returns at once
  and `processOutputFiles()` harvests. If not finished, it polls the row and
  waits. That second branch is what `i2run moorhen ...` uses: i2run
  dispatches immediately, so the plugin is the thing that blocks, exactly
  like `i2run coot1`.

**Why this is the right shape for the desktop, judged on the desktop
alone.** The web deployment is secondary and must not add complexity here;
the session row has to earn its place on desktop cleanliness and
robustness, and it does:

- *Saved work survives an app quit.* Quitting the app kills the Django
  process tree, and `run_job_safe.sh` marks any job whose process died as
  `FAILED`. With a plugin that blocks for the session, a user who has saved
  two models and quits loses them from the database (the files stay on disk,
  unharvested). With the session row there is no process to kill; the job is
  still `RUNNING` on the next launch, the window reopens from the job menu,
  and Finish harvests. This is the one property the Coot tasks lack that
  users will notice.
- *No idle interpreter per window.* A blocking plugin is a full
  `ccp4-python` with Django loaded, doing nothing, for as long as the window
  is open. Coot pays this because Coot must be a child process; Moorhen
  has no such need.
- *One glean path, unchanged.* The runner still does validity, `track_job`,
  `process()`, glean, `FINISHED`. Nothing about how outputs become `File`
  rows is new.
- *Cancel needs nothing new.* The existing endpoint already handles a job
  with no `process_id`.

The complexity it adds is two small hooks and one branch: Run skips
dispatch for an interactive task, Finish dispatches, and `startProcess()`
waits if the session is not yet finished. The simpler "plugin always waits"
shape, identical to `coot1`, remains a legitimate choice if those three are
judged too much; its desktop costs are exactly the first two bullets above.

This gets the good property of both shapes considered earlier. There is no
idle process per open window, a `RUNNING` session survives an app quit or a
container restart (it is just a row and a directory, and the window can be
reopened from the job menu), the CLI still blocks, and there is still **one
glean path**: the runner's. The only new runner-side behaviour is "do not
dispatch on Run for interactive tasks", a switch in `run_job_context_aware`
keyed on the task's registration.

**Ending the session: only the user ends it.** No timer decides anything
on the desktop. A laptop lid closed on an open window stops the heartbeat
for hours; if a timeout had harvested and finished the job in the meantime,
the window on waking would find its saves refused, and work done after the
lid closed would be lost. The rules are therefore:

1. *Finish* — from the window's "Finish session" button, or from the job
   menu's "Finish session" entry on any `RUNNING` interactive job. Marks the
   row and dispatches the job, which harvests and succeeds. If nothing was
   ever saved, Finish does not dispatch: it sets `MARK_TO_DELETE`, so an
   "opened and closed" session does not litter the job list.
2. *Window closed* — `beforeunload` sends a beacon saying the window
   detached. If nothing was saved, the session is finished as in rule 1
   (a status flip, no process, safe even during an app quit). If something
   was saved, the session stays open: the job remains `RUNNING`, the job
   menu offers "Open session window" to reconnect and "Finish session" to
   harvest. Closing a window is not finishing.
   `navigator.sendBeacon` cannot set headers, but the Next proxy already
   accepts the token as an `access_token` query parameter (for anchor-click
   downloads), and on Azure the Easy Auth header is added by the platform,
   so the beacon authenticates either way. Best-effort; if it is lost, the
   session simply stays open, which is the safe outcome.
3. *Cancel* — from the job menu, as today. The endpoint already handles a
   job with no process; it needs only to mark the session finished. Nothing
   is gleaned, which is the rule `coot1` lives under.

The heartbeat (every 10 s from an attached window) therefore decides
nothing. It is information: the job menu uses "no heartbeat in the last
30 s" to offer "Open session window" rather than "Bring window to front",
and the main window can show a toast for a `RUNNING` interactive job that
no window is attached to, which is also how a CLI-started session gets
noticed. This keeps the desktop rule simple: a session ends when the user
says so, and everything the user saved is always still there to harvest.

**The stale-job sweep must know about sessions.** `cleanup_stale_jobs`
fails any job `RUNNING` for more than two hours. It is invoked only from
`server/worker.py`, the Azure worker, so the desktop never runs it; but an
open modelling session is routinely longer than two hours, and so is a Coot
1 session, so the sweep must skip interactive jobs whose session is open.
On Azure, where nobody owns the machine, the sweep may additionally apply
rule 1 to sessions with no heartbeat for 24 hours; that is a tidy-up on a
shared server, not a desktop behaviour.

### Choosing the data to load at startup

The server, not the window, decides what to load. At Run, the API turns the
job's `params.xml` into a **load plan** and stores it with the session:
an ordered list of `{kind, file_id, label, sub_type}` where `kind` is one of
`dictionary`, `coordinates`, `map_2fofc`, `map_fofc`, `map_anom`, `map`.
`GET jobs/{id}/interactive_session/` hands it to the window.

This is not new code. `cootbridge/api_client.py` already does the
`params.xml` walk with the right ordering (`INPUT_PARAM_KINDS`: dictionaries,
coordinates, 2Fo-Fc, Fo-Fc, anomalous) for the Coot bridge, and it is
stdlib-only, so the CCP4-free API process can call it. The producer is that
function with two entries added to `INPUT_PARAM_KINDS`: `DICT_LIST` as
`dictionary` and `MAPIN_LIST` as `map`, the latter with its sub-type carried
through for mask handling.

The window applies the plan with the loaders it already has: `fetchDict`,
`fetchMolecule`, `fetchMap` (whose `mapSubType` argument is exactly the
2Fo-Fc / Fo-Fc / anomalous distinction) and `fetchMapFile`. That is the same
dictionaries-first order `fetchJobFiles` uses, expressed by the server so
both viewers agree with it.

A load plan was chosen over a startup scene on purpose. A scene would let
the server dictate representations and view, which is what `ccp4mg_general`
does by hand with `script.mgpic.xml`, but it would tie the task's startup
contract to the scene format before that format is settled. When it is, the
plan can be replaced by a scene with the same file references and the
window's existing scene resolver takes over; nothing else changes.

The window's existing project browser stays. A file loaded from the browser
mid-session is not an input of the job. A later refinement can append it to
the matching input list through `set_parameter` while the job is `RUNNING`,
which is provenance the Coot tasks do not record either.

### Saving into the job: the drop contract over REST

This is the missing verb. The window must not use `upload_file_param` for it.
That endpoint is the *import* path: it copies into `CCP4_IMPORTED_FILES`,
writes a `FileImport` row, deduplicates against earlier imports and replaces
whatever the parameter held. An output of the session is not an import; it
must be set on `outputData` by the plugin and gleaned once, by the one glean
path, when the job succeeds.

So the window writes into the drop directory and the plugin harvests, which
is the `output<N>` contract every GUI task already uses:

| Endpoint | Method | Does |
|---|---|---|
| `jobs/{id}/interactive_session/` | GET | `{state, load_plan, outputs_so_far}` from the session row and the drop directory |
| `jobs/{id}/interactive_drop/` | POST, multipart | writes `MOORHEN_FILE_DROP/output<N>.<pdb|cif>` via `api_client.next_output_number`; body carries `kind` (`model`, `dictionary`) and an optional `annotation`; returns `{number, path}` |
| `jobs/{id}/interactive_heartbeat/` | POST | updates `last_heartbeat` on the session row |
| `jobs/{id}/interactive_finish/` | POST | `{"finished": bool}`; marks the row and dispatches the job (or `MARK_TO_DELETE` if nothing was saved) |

All four refuse unless the job's task is registered as interactive and its
status is `RUNNING`. They touch the database and the project store only,
need no CCP4, and are fully exercisable from a test without a browser,
which is the point of the i2run-tier test below. The drop directory is
written by the API process into the job directory, which is the same place
`upload_file_param` already writes, so it needs no new storage assumption.

Harvest reuses `coot1.processOutputFiles` wholesale: `harvestable_outputs`,
the model-versus-dictionary content sniff in `cootbridge/harvest.py`, the
`makeItem()`/`pop()` list handling (not `set(slice)`, which drops annotation
and sub-type), and `mergeDictToProjectLib`. There is no new kind to harvest. That harvest
should be lifted out of `coot1.py` into `cootbridge/harvest.py` as
`harvest_session_outputs(...)` and called from `coot1`, `coot_rebuild` and
`moorhen` alike.

Client side, the "Save to this job" button reuses the serialisation in
`push-to-ccp4i2-panel.tsx` (`getAtoms` plus `detectCoordinateFormat`),
which should be extracted into a shared `serialiseMolecule()` so push and
save cannot drift. The annotation defaults to the molecule name.

### Launching the window: the route is the job

The moorhen-page router already resolves what to show from the URL:
`file-by-id/[id]`, `job-by-id/[id]`, `project/[id]`, `campaign/[id]`. A
session is one more segment, `/ccp4i2/moorhen-page/session/[jobId]`, and
**the URL is the whole job-specifying part of the launch**. Every way of
starting a session reduces to "make the app open that route", so nothing
downstream needs to know how it was launched.

The segment is a sibling of `job-by-id`, not a flag on it: `job-by-id` loads
a *finished* job's *output* files (the menu gates it on status 6), whereas a
session loads a *running* job's *inputs* from the load plan and adds the
save, finish and heartbeat behaviour. The client component reads `jobId` from
the route, fetches `interactive_session/`, and hands `moorhen-wrapper.tsx` a
`session` prop alongside the existing `jobId` one.

Three ways to reach the route:

**From the app.** Run is posted from four places (`job-context-menu.tsx`,
`job-card.tsx`, `tool-bar.tsx`, `new-project-content.tsx`). They should share
one `runJob()` helper, and that helper opens the session route when the
task is interactive, using the same `window.open` the preview uses. The
job context menu gains "Open session window" and "Finish session" for a
`RUNNING` interactive job; the first is the reconnect path and the path for
CLI-started jobs, the second is how a session whose window was closed with
saved work gets harvested.

**From the command line.** `i2run moorhen ...` dispatches at once and the
plugin blocks in `startProcess()`, but it cannot open a renderer window
itself. The app can
be told to open a route, with one small Electron change and no protocol
registration:

- Electron main takes `app.requestSingleInstanceLock()` and, on
  `second-instance`, reads `--open-route <path>` from the forwarded argv and
  calls `createWindow` on it. This is the standard single-instance idiom,
  works on all three platforms, and is not Moorhen-specific: any route the
  app serves can be opened this way.
- Electron main exports `CCP4I2_DESKTOP_EXECUTABLE` (`process.execPath`) into
  the Django child's environment, next to `CCP4I2_LOCAL_SESSION_TOKEN`. The
  handshake already proves this reaches the job process.
- The plugin, if that variable is set, spawns
  `[$CCP4I2_DESKTOP_EXECUTABLE, "--open-route",
  "/ccp4i2/moorhen-page/session/<id>"]` once, then waits. If it is not set (a
  bare `runserver`, a web deployment), it logs the route and waits for a
  window to attach by any other means.

In a web deployment the same route works in an ordinary browser, because
auth there is the web session, not the local token. On the desktop a plain
browser cannot attach, because the token lives in the preload; that is a
feature.

A cheaper fallback that needs no Electron change: the main window already
polls jobs, so it can show a toast "Moorhen session requested for job N" for
any `RUNNING` interactive job with no heartbeat. Worth having anyway for the
case where the user closed the window and forgot the job is open.

### Reports

`coot1` writes no program XML and its report says "Happily finished".
`coot_rebuild` writes `number_output_files` and has a real report. The
`moorhen` plugin should do the latter, and can offer
a running report (`runningReport=True`, `watchedFile` on the drop directory)
that shows "session open, N models saved" while the job is `RUNNING`.

### The existing "Open in Moorhen" preview

Keep it. It is the analogue of `files/{id}/preview/` with `viewer: coot`:
cheap, creates no job, right for looking. The recorded task is for building.
The job context menu therefore gains a second entry, "Model build in
Moorhen", which creates a `moorhen` task with `context_job_uuid` set so the
inputs auto-populate from the job, runs it, and opens the window. That is
one `create_task` plus the shared `runJob()`; no new API.

## The in-browser (Azure) deployment

The web deployment is secondary to the desktop. This section records that
the desktop design, as chosen above on desktop grounds, also serves the
deployment where the API, the worker and the browser are three different
machines, and it lists what would have gone wrong if the simpler shape had
been chosen. Nothing here is permitted to add a code path the desktop does
not need.

**What the session row buys on Azure.** In Azure mode
`run_job_context_aware` queues the job to a worker over Service Bus. A
plugin that blocked for the length of a modelling session would occupy a
worker slot per open browser tab, on a machine the browser cannot talk to.
Deferring dispatch to Finish means the worker sees a job that harvests in
seconds, and nothing waits anywhere. The choice was made for the desktop
(saved work surviving a quit); Azure gets it for free.

**Storage.** The drop directory is written by the API process into the job
directory. That is the mount `upload_file_param`
already writes and the worker already reads, so the design assumes nothing
the import path does not. Worth confirming on each instance that the
server and worker containers mount the same project store, because a
session saved on one and harvested on another is the whole mechanism.

**Auth.** The heartbeat and the close beacon go through the Next proxy,
which takes the token from the Authorization header, the Azure Easy Auth
header, or an `access_token` query parameter. On Azure the platform adds
the Easy Auth header to every request including beacons, so the query
parameter fallback is not needed there and should not be used, to keep
tokens out of access logs.

**Launch.** Nothing from PR 3 applies. The browser opens the session route
itself when the user presses Run; there is no command line and no argv.
The `--open-route` argument and `CCP4I2_DESKTOP_EXECUTABLE` are desktop
conveniences the web build never sees.

**Ownership and concurrency.** Several tabs, or several users with access
to the project, can attach to one session. The session row records who
requested it; saves record nothing further, because the job's owner is the
job's owner. If per-save attribution ever matters (a campaign where several
people model the same series), the drop endpoint can stamp the annotation
with the requesting user, which the harvest already carries through.

**Campaigns.** The campaign Moorhen wrapper already creates and runs
`servalcat_pipe` jobs from inside the viewer. "Model build in Moorhen" from
the campaign tables should create the recorded task the same way, with
`context_job_uuid` pointing at the member project's latest refinement, so
that a campaign's modelling history is a chain of `moorhen` and
`servalcat_pipe` jobs rather than refinements with untraceable inputs.

**Stale sweep.** `cleanup_stale_jobs` runs from `server/worker.py` at
startup with a two-hour threshold. Without the exemption above, every
worker restart would fail every open session on the instance. The exemption,
and the 24-hour tidy-up of heartbeat-less sessions, are filters on the
sweep, not a change to the desktop.

## Shared code to extract first

- `cootbridge/harvest.py: harvest_session_outputs(work_dir, drop_dir,
  xyz_list, dict_list, scene_list, annotate)` from `coot1.processOutputFiles`.
- `cootbridge/api_client.py: load_plan(...)` extended for `MAPIN_LIST`;
  it is already CCP4-free and unit tested against `params.xml` fixtures.
- `client/renderer/lib/moorhen-serialise.ts: serialiseMolecule(mol)` from
  `push-to-ccp4i2-panel.tsx`.
- `client/renderer/lib/run-job.ts: runJob(job)` from the four call sites.

`ccp4mg_general` is the natural follow-up: it predates the bridge, builds its
scene by hand and re-implements `output<N>` in `ccp4i2CCP4MGInterface.py`. It
is not in scope here but should move onto the same harvest helper.

## Work packages

Three pull requests, each shippable alone.

**PR 1, server (no client needed).** Task, def.xml, registry entry, the
session row and the deferred-dispatch switch in `run_job_context_aware`,
the `startProcess()` wait, the four endpoints, the stale-sweep exemption,
the load-plan producer and harvest extraction. No new mime type and no
schema change. Tests: unit tests for the load-plan producer and
harvest helper (CCP4-free), API unit tests for the endpoint guards, and two
i2run-tier tests: one that runs the job through the API (Run, drop a
demo-data PDB, Finish) and asserts `XYZOUT[0]` is gleaned with its
annotation, and one that dispatches first, as i2run does, and finishes the
session from a thread while the plugin is waiting. Together they prove both
branches of `startProcess()` without a browser.

**PR 2, client.** Session route and session mode in the wrapper, save and
finish controls, heartbeat and beacon, `runJob()` extraction with the
interactive hook, "Open session window" and "Model build in Moorhen" menu
entries, task interface (a copy of `coot1.tsx` with `DICT_LIST` and `MAPIN_LIST`), chooser
entry. Manual test: run from the task panel, save twice, finish, see two
`XYZOUT` files in the job's file list; close a window without
finishing and see the job stay `RUNNING` with "Finish session" offered;
open and close a session without saving and see `MARK_TO_DELETE`.

**PR 3, Electron.** Single-instance lock, generic `--open-route` argv handling,
`CCP4I2_DESKTOP_EXECUTABLE` export, plugin-side spawn. Manual test: with the
app open, `i2run moorhen --XYZIN_LIST ...` from a terminal opens the window
and blocks until Finish.

## Decisions

Taken 2026-09-16, after the survey and the two revisions above.

1. **Session row plus deferred dispatch**, not an always-waiting plugin.
   Decided on desktop grounds: saved work survives an app quit, no idle
   interpreter per window, one glean path, Cancel unchanged.
2. **A list of dictionaries** (`DICT_LIST`), with the one-line clone rule
   from `coot1`'s `DICT`.
3. **No scene capture for now.** Startup is a load plan, not a scene; the
   scene output waits for the format to be ratified upstream.
4. **No abandonment timeout on the desktop.** Sessions end by Finish,
   by closing an empty window, or by Cancel; a closed window with saved work
   leaves the job `RUNNING` for reconnect or Finish from the job menu. The
   heartbeat only informs the job menu and the toast. Azure may tidy
   heartbeat-less sessions after 24 hours.

## Side findings from the survey

Worth fixing regardless of this design:

- `COOT_EXECUTABLE` does not reach the `coot1` task: its `TASKCOMMAND` is
  `coot-1`, and `program_discovery._EXECUTABLE_PREF` only maps `coot`.
- `ASYNCHRONOUS = True` on the GUI wrappers is dead for top-level jobs:
  `async_run_job.create_plugin_for_job` forces `doAsync = False`, and the
  blocking `subprocess.run` is what waits for the GUI. If the async path is
  ever revived, `process()` calls `postProcess()` unconditionally and would
  double-harvest.
- The Moorhen `job-by-id` route loads through `fetchFile`, which has no
  mmCIF branch and no `directory == 1` filter, so it behaves differently from
  the in-viewer "load all job outputs" path that uses `fetchJobFiles`.
