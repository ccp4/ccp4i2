# Authoring a Task

This is the **end-to-end guide** to adding a new task (a wrapper around a
crystallographic program) to CCP4i2. It ties together the individual references
so you can follow one path from nothing to a working, UI-rendered task.

A task is four things:

1. A **`def.xml`** — declares the data model (inputs, outputs, parameters).
2. A **wrapper** — a `CPluginScript` subclass that builds the command line and
   harvests outputs.
3. A **registration** — one entry in the `TASKS` dict.
4. (Optional) A **custom UI** — only if the auto-generated one isn't enough.

> **You get a working interface for free.** If you stop after steps 1–3, the
> frontend auto-renders the task from its `def.xml` via `GenericInterface`. Write
> a bespoke UI (step 4) only when you need conditional visibility, custom layout,
> or reactive behaviour.

---

## Step 0 — Choose how to declare parameters

There are **two first-class ways** to declare a task's parameters:

| | Use when | Guide |
|---|---|---|
| **`def.xml`** | A conventional program with a manageable set of options | [def.xml Reference](def-xml-reference.md) |
| **PHIL** | The tool already exposes a `master_phil` (Phenix, PhaserTNG, DIALS…) | [PHIL_TASK_GUIDE.md](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md) |

With PHIL, `PhilPluginScript` ingests the tool's own parameter definitions at
run time and builds the **same** parameter model — rendered by the **same** UI,
with expert-level filtering for free — so you don't hand-transcribe hundreds of
options, and the wrapper auto-syncs when the tool updates. Even a PHIL task keeps
a small `def.xml` for its `inputData`/`outputData`.

The rest of this guide follows the `def.xml` path; the PHIL guide slots in at
steps 1–2 where noted.

---

## Step 1 — Write the `def.xml`

Create `wrappers/<Name>/script/<Name>.def.xml` with three containers —
`inputData`, `outputData`, `controlParameters`. Each parameter is a
`<content id="…">` with a `<className>` (`CData` type) and `<qualifiers>`.

See the **[def.xml Reference](def-xml-reference.md)** for the full skeleton, the
common `CData` types, the qualifiers (`mustExist`, `allowUndefined`, `default`,
`enumerators`, `min`/`max`, `sameCrystalAs`, …), and a complete worked example
(`freerflag`).

*(PHIL path: put only `inputData`/`outputData` here; declare control parameters
via `get_master_phil()` in the wrapper.)*

## Step 2 — Write the wrapper

Create `wrappers/<Name>/script/<Name>.py` — a `CPluginScript` subclass. The
essentials:

- Class constants (`TASKNAME`, `TASKCOMMAND`, whether it runs async) — plus
  `AUXILIARY_PROGRAMS` and `LOG_FAILURES` if they apply; see
  [Declaring the programs your task needs](#declaring-the-programs-your-task-needs).
- `process()` — build the command line from the container and run the program.
- `processOutputFiles()` — register/convert outputs so they can be gleaned.
- Optionally override `validity()` (fast, polled during editing) and
  `runTimeValidity()` (expensive, at submission) — see the validation section of
  the project `CLAUDE.md` and [Validity Patterns](../mddocs/pipeline/VALIDITY_PATTERNS.md).

Also add empty `__init__.py` files in the task dir and its `script/`.

References:
- Patterns & lifecycle: [QUICK_REFERENCE.md](../mddocs/QUICK_REFERENCE.md),
  [pipeline_best_practices.md](pipeline_best_practices.md),
  [ERROR_HANDLING_PATTERNS.md](../mddocs/pipeline/ERROR_HANDLING_PATTERNS.md).
- Simple examples to copy: `wrappers/freerflag/`, `wrappers/nucleofind/`.
- **PHIL path:** [PHIL_TASK_GUIDE.md](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md)
  — includes [porting a task that already has a generated `.def.xml`](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md#porting-a-task-that-already-has-a-generated-defxml) and the interface half of a PHIL task.
  (subclass `PhilPluginScript`; implement `get_master_phil()`,
  `get_shim_definitions()`, `get_command_target()`, `processOutputFiles()`).
  Reference wrappers: `wrappers/phasertng_picard/`, `wrappers/phasertng_riker/`.
- **Web-service path:** [WEB_API_TASK_GUIDE.md](../server/ccp4i2/wrappers/WEB_API_TASK_GUIDE.md)
  — tasks that submit to a remote API needing a token or password. Reference
  wrapper: `wrappers/pdb_redo_api/`.

> **A credential is never a task parameter.** If your task needs a token,
> password or passphrase, it does **not** go in the `def.xml`. Task parameters
> are written to `input_params.xml`, stored in the project database, and
> included in every exported project zip. Register the credential in
> [`config/credentials.py`](../server/ccp4i2/config/credentials.py) and read it
> with `credentials.get_credential(...)`; see the guide above.

> **Outputs are persisted automatically.** The gleaner saves any `outputData`
> `CDataFile` that is set and exists on disk after the job — there is no
> `saveToDb` gate.

### Declaring the programs your task needs

`TASKCOMMAND` names the one program the base class launches. Two class
attributes cover what it does not.

**`AUXILIARY_PROGRAMS`** — binaries your task drives from *inside*
`TASKCOMMAND`, or (for a pipeline that overrides `process()`) binaries it runs
itself. Without this they are invisible to CCP4i2 until the job fails:

```python
class arcimboldo(CPluginScript):
    TASKCOMMAND = 'ARCIMBOLDO_LITE'
    # ARCIMBOLDO_LITE is a driver: phaser and shelxe do the real work, and
    # shelxe is licence-restricted so it is NOT shipped with CCP4.
    AUXILIARY_PROGRAMS = ('phaser', 'shelxe')
```

Declaring them buys two things:

- The pre-run check (`runTimeValidity`) reports a missing one *before* the job
  starts, naming the program and pointing at Preferences → Program locations.
  Where CCP4i2 will run the job itself — desktop/Electron, i2run — this
  **blocks** submission, because a missing binary there is a certain failure.
  Where the job is queued for a remote worker, or this is the CCP4-free API
  server, it stays advisory: "not found here" says nothing about the run host.
- They appear on the Preferences → Program locations page, so a user can point
  CCP4i2 at their own copy.

**Resolve paths, never hardcode them.** If you write a program path into a
config file, get it from `resolve_program()` so the user's `SHELXDIR` /
`exePaths` preferences are honoured — a hardcoded `$CCP4/bin/<name>` silently
ignores a correctly-configured install:

```python
from ccp4i2.config.program_discovery import resolve_program

path = resolve_program('shelxe') or os.path.join(ccp4_home, 'bin', 'shelxe')
```

**`LOG_FAILURES`** — patterns that mean the run failed even though it exited 0.
A great many CCP4 programs report a fatal error to their log and exit cleanly
anyway, so the exit code is not evidence of success. Declare a table rather
than writing a log parser:

```python
LOG_FAILURES = (
    # (pattern, error code, details)  -- details=None means "use what the
    # pattern captured": group(1) if it has one, else the whole line.
    (r'^\s*FATAL(?:\s+ERROR)?\b[:\s]*(.*)$', 301, None),
    (r'minimum or no description', 302,
     'No geometry dictionary for a ligand -- run Make Ligand first'),
)
```

The base `postProcessCheck()` applies these after a successful exit, matching
per line with ANSI colour codes stripped. Any match fails the job and files the
message in the error report. Override `reportLogFailures(failures)` if your
report class needs the text in `PROGRAMXML` to render something better than an
empty panel. Tasks that declare no patterns are unaffected — the scan is
skipped entirely.

## Step 3 — Register the task

Add **one entry** to the `TASKS` dict in
[`server/ccp4i2/core/tasks.py`](../server/ccp4i2/core/tasks.py). That is the
entire registration — there is **no registry-regeneration step** (the old
`plugin_registry.py` / `plugin_lookup.json` scan was removed).

```python
"mytask": Task(
    title="Human-readable title shown in the task chooser",
    description="One-line description",
    shortTitle="Short label",
    pluginPath="ccp4i2.wrappers.mytask.script.mytask:mytask",  # module:class
    defXmlPath="wrappers/mytask/script/mytask.def.xml",
    reportPath="ccp4i2.wrappers.mytask.script.mytask_report:mytask_report",  # optional
    # runningReport=True,   # if the report updates live while the job runs
    # watchedFile="...",    # file to watch for a running report
    # ccp4_free=True,       # ONLY if execution needs no CCP4 binary/$CCP4 (pure gemmi/python)
),
```

Fields (from the `Task` dataclass): `title`, `description`, `shortTitle`,
`pluginPath`, `defXmlPath` are the usual ones; `reportPath` only if you wrote a
report class; `runningReport`/`watchedFile` for live reports; `ccp4_free` only
for tasks that provably need no CCP4 install (verified by the CCP4-free guard).

At this point the task runs (via i2run and the API) and **renders a default UI
automatically**.

## Step 3.5 — Cite the program your task runs

If your task wraps a published program, its report should say so. Reports build
their bibliography from the MedLine files in `server/ccp4i2/references/`, and a
task that cites nothing renders an empty references section — which used to
happen silently and now logs a warning naming your task.

There is one thing to get right, and it is not the obvious one: **the lookup is
by citation key, not by task name.** Register the mapping in
`server/ccp4i2/core/citations.py`:

```python
TASK_CITES = {
    ...
    "my_new_task": ["theprogram"],   # -> references/theprogram.medline.txt
}
```

- Reuse an existing key wherever the program is already cited — don't copy a
  reference file per task. The map is one-to-many, so a task that drives several
  programs cites all of them.
- Only add a new `references/{key}.medline.txt` when the program is not yet
  cited anywhere, and name it for the *program*, not your task.
- If your task is pure i2 plumbing with nothing to cite, add it to `NON_CITABLE`
  in the same file. That is the only correct way to be quiet.

Full rules, and the three bugs that motivated them, are in
[Bibliography](BIBLIOGRAPHY_PLAN.md#how-a-citation-resolves).

## Step 4 — (Optional) Custom frontend interface

Only if the auto-generated interface isn't sufficient. Register a React
interface in
[`task-container.tsx`](../client/renderer/components/task/task-interfaces/task-container.tsx)
and add the task to a category in
[`task-chooser.tsx`](../client/renderer/components/task/task-chooser.tsx).

The definitive, worked guide (layout, conditional visibility, `onChange`,
auto-sync, pitfalls, two full examples) is the **[Task Interface Implementation
Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)** —
read it before building an interface.

---

## Step 5 — Test it

Add an end-to-end CLI test under `server/ccp4i2/tests/i2run/` (one i2run job per
test — running two in one test can hit a read-only DB). See
[writing-i2run-tests.md](writing-i2run-tests.md). Run:

```bash
cd server
ccp4-python -m pytest ccp4i2/tests/i2run/test_mytask.py -v
```

---

## Checklist

- [ ] `wrappers/<Name>/script/<Name>.def.xml` (three containers)
- [ ] `wrappers/<Name>/script/<Name>.py` (`CPluginScript` / `PhilPluginScript`)
- [ ] `__init__.py` in the task dir and its `script/`
- [ ] One entry in `core/tasks.py` `TASKS`
- [ ] (optional) report class + `reportPath`
- [ ] (optional) React interface + task-chooser category
- [ ] (if it calls an authenticated service) credential registered in
      `config/credentials.py` — **not** in the `def.xml`
- [ ] (if it needs binaries beyond `TASKCOMMAND`) `AUXILIARY_PROGRAMS`
      declared, and any program paths written into config files resolved with
      `resolve_program()` rather than hardcoded
- [ ] (if the program can exit 0 after a fatal log message) `LOG_FAILURES`
      declared
- [ ] (if it wraps a published program) citation registered in
      `core/citations.py` — `TASK_CITES`, or `NON_CITABLE` if there is nothing
      to cite
- [ ] An i2run test

## See also

- [def.xml Reference](def-xml-reference.md)
- [PHIL Task Guide](../server/ccp4i2/wrappers/PHIL_TASK_GUIDE.md)
- [Web-API Task Guide](../server/ccp4i2/wrappers/WEB_API_TASK_GUIDE.md) ·
  [Handling Secrets](CREDENTIALS_DESIGN.md)
- [Bibliography — how a citation resolves](BIBLIOGRAPHY_PLAN.md#how-a-citation-resolves)
- [Task Interface Implementation Guide](../client/renderer/components/task/task-elements/TASK_INTERFACE_IMPLEMENTATION_GUIDE.md)
- [Pipeline Best Practices](pipeline_best_practices.md)
- [Error Handling Patterns](../mddocs/pipeline/ERROR_HANDLING_PATTERNS.md) ·
  [Validity Patterns](../mddocs/pipeline/VALIDITY_PATTERNS.md)
