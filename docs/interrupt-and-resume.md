# Interrupting and Resuming Jobs

A design note on three job-lifecycle verbs CCP4i2 conflates — **Cancel**,
**Interrupt** and **Resume** — arising from an investigation into the DUI2 task.

**Baseline:** commit `b1cf282f3` (branch `django`), investigated 2026-08-26
against `ccp4-20260702` (DUI2 `2026.6.23`, modelcraft `6.1.1`).
**Status:** design note. Nothing here is implemented. The live defects in
[Defects found along the way](#defects-found-along-the-way) are real today and
are worth fixing regardless of what is decided about the verbs.

Companion documents:

- [docs/authoring-a-task.md](authoring-a-task.md) — the task-authoring spine.
- [docs/def-xml-reference.md](def-xml-reference.md) — output container shapes,
  which turn out to decide Resume-ability.
- [docs/error-handling-remediation.md](error-handling-remediation.md) — the
  neighbouring tracker; the glean-gating finding below sits next to its concerns.

## The finding in one paragraph

CCP4i2 has one verb, *Run*, and one terminal status, *Interrupted*, that is only
ever reached by *Cancel*. But three different things want to happen to a job:
abandon it (**Cancel**), stop it and keep what it has built (**Interrupt**), and
re-enter it to carry on (**Resume**). Legacy CCP4i2 implemented Interrupt — and
gleaned output files on *every* finish status, so an interrupted job published
its artefacts. The Django branch gleans only on success, so an interrupted job
publishes nothing, which removes the reason to interrupt at all. Meanwhile
`POST /jobs/{id}/run/` has no status precondition, so re-running a *finished* job
in place is already possible and would silently rewrite the `File` rows that
downstream jobs point at. The verbs need separating, and two mechanical rules
decide which tasks can support which.

## Three verbs, not two

| Verb | Publishes artefacts? | Job afterwards | Status here |
|---|---|---|---|
| **Cancel** — abandon it, I don't want it | no | dead | exists, misnamed as `INTERRUPTED` |
| **Interrupt** — stop, I want what you've built so far | **yes** | resumable or terminal | legacy had it; lost |
| **Resume** — re-enter and carry on | — | continues | never existed |

The hazards live in different places for each, which is why conflating them is
expensive:

- **Cancel** is safe by construction — it publishes nothing.
- **Interrupt** is where the overwrite hazard lives, because it is the verb that
  publishes *and* leaves the job resumable.
- **Resume** is safe if and only if the artefacts it will overwrite have not been
  consumed.

## What the code does today

**Gleaning is gated on success.** All three call sites require it:

- [`db/async_db_handler.py:331`](../server/ccp4i2/db/async_db_handler.py#L331) —
  `if status == models.Job.Status.FINISHED and container is not None:`
- [`db/async_db_handler.py:971`](../server/ccp4i2/db/async_db_handler.py#L971) —
  `if plugin_status == CPluginScript.SUCCEEDED:`
- [`db/async_db_handler.py:1070`](../server/ccp4i2/db/async_db_handler.py#L1070) —
  same, for subjobs

The `status_map` at
[`:1061-1065`](../server/ccp4i2/db/async_db_handler.py#L1061-L1065) does not even
map `INTERRUPTED` to a status. **An interrupted job publishes nothing.**

**`cancel` is the only producer of `INTERRUPTED`.**
[`api/JobViewSet.py:2029-2088`](../server/ccp4i2/api/JobViewSet.py#L2029-L2088)
kills the process tree with `proc.kill()` (SIGKILL, no grace period) and sets the
status. It is Cancel wearing Interrupt's name.

**Re-gleaning is idempotent, which is a hazard as well as a convenience.**
`register_output_file` uses `get_or_create` keyed on `(job, job_param_name)` and
**updates the existing row in place**
([`:378-408`](../server/ccp4i2/db/async_db_handler.py#L378-L408)). So re-running
a job does not duplicate File rows — it silently redirects them. Any downstream
job holding a `FileUse` on that row has its input changed underneath it, with no
record that it happened.

**Nothing gates `run` on status.** Neither
[`api/JobViewSet.py:873`](../server/ccp4i2/api/JobViewSet.py#L873), nor
`run_job_context_aware`, nor the local runner checks the job's current status. So
re-running a FINISHED job in place is possible today via the API. The UI does not
expose it — `client/renderer/components/tool-bar.tsx:212-213` gates the Run
control on `job?.status === 1` (PENDING) — so it is reachable from i2run, scripts
and the campaigns code rather than by clicking.

**Files live in their producing job's directory.**
[`db/models.py:392-397`](../server/ccp4i2/db/models.py#L392-L397):
`path` is `self.job.directory / self.name`, and consuming jobs reference the same
row via `FileUse(role=IN)`. i2 never copies an antecedent's outputs. Note also
that `name` is a bare basename
([`:545`](../server/ccp4i2/db/async_db_handler.py#L545)) — **a File row cannot
point into a subdirectory**, which constrains how session-based tasks can harvest.

## The legacy precedent, and its defect

`main:pipelines/buccaneer_build_refine_mr` declared `INTERRUPTABLE = True` and
`RESTARTABLE = True`, kept its resume state in a `container.interruptStatus`
subcontainer, and restarted from `LASTCYCLE + 1` via `handleRestart()`, with
`parseXmlFromBeforeInterrupt()` so the report continued rather than restarting.

Legacy published on interrupt, deliberately:

- `main:core/CCP4PluginScript.py:1031-1060` — `reportStatus` calls
  `setOutputFileContentFlags()` (walking `outputData`, `saveToDb()` +
  `setContentFlag()` on every set `CDataFile`) for *every* finish status.
- `main:core/CCP4PluginScript.py:2085-2093` — the glean inside `updateJobStatus`
  is gated on `if container is not None:` and nothing else. It runs before the
  status is written and irrespective of what it is.
- `main:core/CCP4PluginScript.py:1266-1273` — on `INTERRUPTED`, `saveParams`
  additionally writes the `interruptStatus` subcontainer to an INTERRUPT params
  file, read back by `loadInterruptStatus()` when `RESTARTABLE`.

And buccaneer's own interrupt branch calls `self.alignment()` first — the same
call the success path makes immediately before `reportStatus(SUCCEEDED)`, and one
that reads `self.container.outputData.XYZOUT`. The pipeline completed its outputs
and its report *before* declaring INTERRUPTED. Interrupt was designed to make
artefacts available.

**Its defect:** buccaneer's outputs are all scalar — `XYZOUT` (`CPdbDataFile`),
`FPHIOUT` and `DIFFPHIOUT` (`CMapCoeffsDataFile`), `ABCDOUT` (`CPhsDataFile`).
Interrupt published them; a downstream job could consume them; restart resumed at
`LASTCYCLE + 1` and overwrote all four in place. Declaring a scalar-output task
both interruptible and restartable is the bug.

One loose end, unverified: buccaneer's `reportStatus(INTERRUPTED)` has no `break`
or `return` after it, and `reportStatus` returns normally rather than terminating.
Whether the cycle loop actually stopped depended on the parent tearing the child
process down. Moot here — the pipeline is not in this branch — but it smells.

## Why nothing consumes the mechanism now

The mechanism has **zero live consumers** on this branch.

- **modelcraft pushed the iteration into the program.**
  [`wrappers/modelcraft/script/modelcraft.py`](../server/ccp4i2/wrappers/modelcraft/script/modelcraft.py)
  is 112 lines with five methods and no cycle loop; it passes `--cycles` and
  `--auto-stop-cycles` and lets the binary iterate. There are no cycle boundaries
  left for i2 to hook.
- **The one apparent survivor is dead code.**
  [`pipelines/dr_mr_modelbuild_pipeline/script/dr_mr_modelbuild_pipeline.py:729-790`](../server/ccp4i2/pipelines/dr_mr_modelbuild_pipeline/script/dr_mr_modelbuild_pipeline.py#L729-L790)
  — `startModelCraftProcess` returns unconditionally from both branches of an
  if/else, and the ~50 lines after it (buccaneer's cycle loop, including
  `testForInterrupt()` / `interruptStatus.LASTCYCLE` / `reportStatus(INTERRUPTED)`)
  are unreachable. They also call `plugin.handleRestart()`,
  `plugin.buildRefineCycle()`, `plugin.alignment()` and friends on a *modelcraft*
  plugin, which has none of those methods; it would `AttributeError` immediately
  if it ever ran. It is buccaneer's `startProcess` pasted in and then bypassed.
- **crank2 uses `testForInterrupt()` as a stop signal only**, alongside its own
  `stop_file`. It declares neither `INTERRUPTABLE` nor `RESTARTABLE`.
- **No def.xml anywhere declares an `interruptStatus` container**, so even the
  state-persistence half has no schema behind it.

What remains is `testForInterrupt()` in core (a cancel signal), `Status.INTERRUPTED`
(set only by `cancel`), and an unreachable paste. That is good news for design
freedom: there is no installed base to preserve.

## Two rules

**Rule 1 — Resume requires no dependents.** Because re-gleaning updates File rows
in place, resuming a job whose outputs have been consumed silently redirects a
downstream job's inputs. `find_dependent_jobs`
([`lib/utils/navigation/dependencies.py:17-42`](../server/ccp4i2/lib/utils/navigation/dependencies.py#L17-L42))
already walks exactly the `FileUse` graph needed; enforcement is a guard, not new
machinery.

**Rule 2 — Resume requires append-only outputs.** A list output can append; a
scalar output can only overwrite. No stable-keying trick rescues a scalar, because
there is nowhere to put the new artefact except on top of the old one. So a task
is Resume-able only if every output is list-shaped **and** keyed on something
intrinsic (a node number, a filename index) rather than on list position.

Both rules are mechanically checkable — Rule 1 from the database, Rule 2 from the
def.xml plus the harvest code.

## Per-task assessment

| Task | Launches | Session state | Output shape | Verdict |
|---|---|---|---|---|
| `dials_image` | `dials.image_viewer` | none (`MARK_TO_DELETE`) | none | neither — the job deletes itself |
| `dials_rlattice` | `dials.reciprocal_lattice_viewer` | none (`MARK_TO_DELETE`) | none | neither |
| `qtpisa` | `qtpisa` | none; deterministic from `XYZIN` | list | neither — re-running *is* resuming |
| `ccp4mg_edit_model` | `ccp4mg` | `-norestore`; empty `.mgpic.xml`, never harvested | list, keyed by filename | Resume conceivable, nothing built |
| `coot1` | `coot-1` | none captured | `COutputFileList`, **keyed by position** | Resume conceivable; needs Rule 2 fix |
| `coot_rebuild` | `coot` | **writes `state<N>.py` per save; has a consuming input** | `COutputFileList`, keyed by filename | **Resume: strong, mostly built** |
| `dui` | `dui2` | `run_data` DAG — the whole point | `CList`, **keyed by position** | **Resume: strong**; needs Rule 2 fix |
| `modelcraft` | `modelcraft` | none needed | scalars | **Interrupt: nearly free**; Resume not wanted |
| `buccaneer_build_refine_mr` | *(not in this branch)* | `interruptStatus.LASTCYCLE` | **scalars** | the Rule 2 counter-example |

## Modelcraft: Interrupt is nearly free

modelcraft 6.1.1 has no stop file, no `SIGINT`/`SIGTERM` handler, no `atexit`, no
`KeyboardInterrupt` catch. Its `terminate(reason)` (`pipeline.py:61`) fires only
on internal conditions. **But it does not need one to be interrupt-safe.**

`modelcraftxray.py:315-336` — every time the model improves, `process_cycle_output`
writes `modelcraft.cif` and `modelcraft.mtz`, sets `report["final"]`, and flushes
`modelcraft.json`. After any improving cycle a complete, mutually consistent
output set is on disk. Killing the process between cycles leaves exactly what a
shorter successful run would have left.

And the wrapper's `processOutputFiles`
([`modelcraft.py:96-112`](../server/ccp4i2/wrappers/modelcraft/script/modelcraft.py#L96-L112))
harvests precisely those three files and nothing else. **A killed modelcraft could
be harvested by the existing code with no change to the harvest logic.**

What would need doing:

1. **Ungate the glean** for INTERRUPTED. This is the blocker; the harvest is ready.
2. **Guard the harvest.** Killed before the first improving cycle, none of the
   three files exist and `shutil.copy` raises `FileNotFoundError`; `result["final"]`
   is likewise absent. Needs an existence check and a graceful "nothing built yet".
3. **Mind a narrow race.** `cancel` uses SIGKILL. `process_cycle_output` writes the
   cif, then the mtz, then the JSON; a kill inside that window leaves a new cif
   beside an old mtz. An upstream stop-file check at the top of `run()`'s cycle
   loop would avoid it entirely, and is a small, natural ask.

**Continuation is already free and needs no new verb.** `makeCommandAndScript`
passes `--model` when `inputData.XYZIN` is set, so "carry on from where modelcraft
got to" is: clone the job, feed the previous `XYZOUT` as `XYZIN`. Ordinary file
passing, a real `FileUse` edge, full provenance — and it sidesteps Rule 2 entirely,
because a new job never overwrites the old one's scalars.

## DUI: what it needs

Current state (see also the defects section):

- Input `DUI2_RUN_DATA` is a `CString` with `onlyEnumerators=True` and **no
  enumerators**, because the population logic lived in the Qt widget
  (`main:wrappers/dui/script/CTaskdui.py`, from PR #142) and did not survive.
- The React interface renders `itemName="DUI_DIR"` — the pre-#142 name — so the
  task panel shows the literal text "No item".
- `run_data` is `<jobdir>/run_dui2_nodes/run_data`, a JSON manifest whose
  `step_list` nodes carry **absolute** paths (`_base_dir`, `_run_dir`,
  `_lst_expt_in/out`, `_lst_refl_in/out`, `log_file_path`, `_html_rep`,
  `_predic_refl`) into sibling `runN/` directories. The wrapper copies the manifest
  alone, which works only while the previous job directory survives at the same
  absolute path.
- Outputs are `CList`s harvested by `sorted(nodes.glob("**/*.mtz"))` with names
  `"_".join(relative parts)`, so `run10_*.mtz` sorts before `run2_*.mtz` and the
  indices reshuffle. **Violates Rule 2.**

Design options, in increasing order of ambition:

1. **Restore parity.** Fix the itemName and repopulate the enumerators
   server-side. Makes the task work; changes nothing structural.
2. **Patch the paths on copy-in.** Rewrite only the known path-bearing keys from
   harvested typed data, leaving the rest of DUI2's private format untouched, and
   assert every rewritten path resolves before starting. Small contract surface.
3. **Harvest node artefacts as typed data and regenerate `run_data`.** Turns the
   manifest from authoritative into derived state, which is the right direction —
   but note that `File.name` is a bare basename, so `.expt`/`.refl` would have to
   be flattened and copied to the job-dir root, and the DAG topology
   (`cmd_dict_ini`, `full_cmd_lst`, `number`, parent/child links, `bigger_lin`)
   is not file-shaped and would need persisting separately. It is also an
   unversioned contract with DUI2's `_recover_state`.

Option 2 is the recommended middle path.

## Coot: 90% built and unwired

The Resume loop for `coot_rebuild` is nearly complete and entirely disconnected:

- **Producer exists.**
  [`ccp4i2CootInterface.py:228`](../server/ccp4i2/wrappers/coot_rebuild/script/ccp4i2CootInterface.py#L228)
  — `saveToDatabase()` calls `save_state_file_py(dropDir/state<N>.py)` on every
  save-to-i2, paired 1:1 with `output<N>.pdb`.
- **Nothing harvests it.** `processOutputFiles` globs only `output*.pdb` and
  `output*.cif`. The state files die with the job directory.
- **Consumer exists.** Input `COOTSTATEFILE` (`CCootHistoryDataFile`) plus
  `copyStateFile()`, which rewrites the molecule-loading line to repoint at the
  current `XYZIN` — the same path-patching pattern DUI needs.
- **Nothing can feed it.** `COOTSTATEFILE` is `saveToDb=False`, and no task
  anywhere outputs one. Every other `CCootHistoryDataFile` in the tree (refmac,
  servalcat, validate_protein, privateer, edstats, lorestr, i2Dimple, morda) is
  declared as `COOTSCRIPTOUT` in `outputData` — a *coot script* ("here is how to
  view my results"), consumed as `COOTSCRIPTFILE`. A different thing.

Two defects in the dormant path, to fix before trusting it:

- `copyStateFile()` accumulates the patched text into `newText` and then **writes
  `text`**, the unpatched original
  ([`coot_rebuild.py:314-322`](../server/ccp4i2/wrappers/coot_rebuild/script/coot_rebuild.py#L314-L322)).
- It writes `0-coot-history.scm` (Scheme) while what Coot saves is `state<N>.py`
  (Python), and the command line passes `--no-state-script` regardless.

On the plus side, `coot_rebuild` derives its output index from the filename
(`iFile = int(fname[6:-4])`), so it already satisfies Rule 2.

## Defects found along the way

Independent of any decision about the verbs:

1. **`dui.tsx` renders a nonexistent parameter.** `itemName="DUI_DIR"` versus the
   def.xml's `DUI2_RUN_DATA`;
   [`task-element.tsx:240-245`](../client/renderer/components/task/task-elements/task-element.tsx#L240-L245)
   falls through to `<Typography>No item</Typography>`. The DUI task panel is
   visibly broken.
2. **DUI's "continue from a previous session" is dead** — enumerators are never
   populated on this branch.
3. **`POST /jobs/{id}/run/` has no status precondition**, so a finished job can be
   re-run in place via the API and its File rows silently redirected. Worth
   guarding regardless.
4. **`dui` and `coot1` key list outputs by position**, so any re-harvest can
   redirect existing File rows at different files. Latent today (each run is a
   fresh job); a prerequisite for Resume.
5. **`coot_rebuild.copyStateFile()` discards its own patched output** (writes
   `text`, not `newText`).
6. **`dr_mr_modelbuild_pipeline.startModelCraftProcess` carries ~50 lines of
   unreachable code** that would `AttributeError` if reached. Dead-code removal.
7. **`dui_report` is empty** (`addResults()` and nothing else) while DUI2 writes a
   `dials.report.html` per node and a `run_data` DAG with per-node command,
   status, log path and report path — good report material going unused.

## Suggested sequencing

1. Fix defects 1 and 2 — make the DUI task work again. Self-contained.
2. Fix defect 3 — the `run` status guard. Cheap, and it is Rule 1's enforcement
   point whatever else happens.
3. Fix defect 4 — stable output keying for `dui` and `coot1`. A live latent bug
   *and* the precondition for Resume.
4. Decide **Interrupt**: ungate the glean for INTERRUPTED, guard modelcraft's
   harvest, and consider asking upstream for a cycle-boundary stop file. Benefits
   a task people run daily and does not commit to anything on the Resume side.
5. Decide **Resume**, if at all — DUI and `coot_rebuild` are the only genuine
   customers. Wiring Coot's existing producer/consumer is close to free; DUI needs
   option 2 above.

## Open questions

- Should Resume re-enter the *same* job, or is "clone and continue" (as modelcraft
  already does via `--model`) good enough everywhere? Cloning gives real provenance
  edges and sidesteps both rules; the cost is job-list noise and, for DUI, a
  cross-job reference into the previous job's `run_dui2_nodes`.
- Should Interrupt leave a job resumable, or is "publish and close" sufficient?
  Only the former needs Rule 2.
- If Interrupt publishes, what does the UI say about a job that is Finished-ish
  but was stopped early? `Status.UNSATISFACTORY` already exists and may be the
  honest label.
- Does `ccp4mg_edit_model` want any of this, or is its `.mgpic.xml` session
  support not worth wiring?
