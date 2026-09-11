# Driving CCP4i2 as an agent

A briefing for an AI agent — or a developer pairing with one — that needs to
operate CCP4i2 through its API or CLI, e.g. to run a molecular-replacement and
refinement from merged data. Read this once and you should be productive
without rediscovering the mechanics.

> **Living document.** It describes what the CCP4i2 *substrate* lets an agent do
> **today**. The faculties CCP4i2 does **not** yet provide — holding a goal,
> planning across tasks, judging whether a result is right — you must supply
> yourself for now (§6). As we bridge those gaps, the matching sections here
> gain "now built" notes (§7). Keep it honest and current.

---

## 1. What CCP4i2 gives you, and what it doesn't

CCP4i2 is a strong **substrate**, not (yet) an agent. Concretely:

| You get, ready to use | You must bring |
|---|---|
| A deep, validated library of crystallographic operations (tasks + pipelines) | A **goal** and a **success test** ("solved" = ?) |
| Every job's inputs, parameters, outputs and provenance, recorded & reproducible | **Planning** across tasks, and choosing at branch points |
| Strong input/pre-flight **validation** (`validity()` / `runTimeValidity()`) that stops nonsense before it runs | **Judgement** of results — CCP4i2 shows quality metrics; deciding "good / failed / retry" is on you |
| Rich metrics, reports and KPIs on every job | Interpreting those metrics |

Pipelines are your friend: `aimless_pipe`, `phaser_pipeline`, the model-building
pipelines etc. are *pre-baked expert workflows* that make many low-level choices
internally. Prefer a pipeline over hand-chaining wrappers where one exists.

---

## 2. Get connected (desktop app — jobs appear live in the UI)

The desktop app runs a local Django server, **token-authenticated**, on a port
chosen at launch. Both the port and the token are surfaced for you:

> **Help → About CCP4i2** shows copyable fields: **Backend (Django) port**,
> **Front-end (Next.js) port**, **Session user**, **Session token**, and an
> **Example request**. (Requires v3.1.0a51 or later.)

Talk **directly to the Django port** (not the Next proxy):

```bash
export CCP4I2_PORT=<Backend (Django) port from About>
export CCP4I2_TOKEN=<Session token from About>
BASE="http://localhost:${CCP4I2_PORT}/api/ccp4i2"

curl -s -H "Authorization: Bearer ${CCP4I2_TOKEN}" "$BASE/projects/"
```

Two rules that will otherwise cost you an hour:

- **Auth on every request**: `Authorization: Bearer <token>`. Without it you get 401.
- **Trailing slashes**, direct to Django (DRF): `…/projects/`, `…/jobs/47/run/`.
  (The *Next proxy* path `/api/proxy/ccp4i2/…` is the opposite — it rejects a
  trailing slash with a 308. Going direct to the Django port, keep the slash.)

The token grants full access to the user's projects and dies when the app
closes. Never paste it into a chat, a bug report, or a shared document.

---

## 3. The core REST flow

All paths below are under `BASE = http://localhost:<port>/api/ccp4i2`, all
`POST`s carry the Bearer header. Endpoints verified against `api/urls.py`,
`ProjectViewSet`, `JobViewSet`.

| Step | Call | Notes |
|---|---|---|
| Create a project | `POST /projects/` `{ "name": "...", "directory": "<optional path>" }` | `ProjectViewSet` is a ModelViewSet; returns the project incl. `id`. |
| Create a task/job | `POST /projects/{id}/create_task/` `{ "task_name": "phaser_pipeline" }` | Returns the new job (`data.new_job.id`). Task names are the keys of `core/tasks.py`. |
| Set a parameter | `POST /jobs/{id}/set_parameter/` `{ "object_path": "<task>.container.inputData.XYZIN", "value": ... }` | Dot-path into the task's def.xml container. Discover paths from the job's container/params or the def.xml. |
| Import a file into a param | `POST /jobs/{id}/upload_file_param/` (multipart: `object_path`, `file`) | The import path; splits/validates as needed. |
| Pre-flight validation | `GET /jobs/{id}/validation/` and `GET /jobs/{id}/run_time_validation/` | **Check these before running.** Errors here are the guardrail; heed them. |
| Run | `POST /jobs/{id}/run/` (or `/run_local/`) | Submits the job. |
| Poll status + KPIs | `GET /projects/{id}/job_tree/` | Returns the job tree with statuses and KPIs (`float_values` / `char_values`) — your window on progress and quality. |
| Clone a job (as template) | `POST /jobs/{id}/clone/` | Re-run with tweaks: clone, then `set_parameter`, then `run`. |

Job status codes you'll poll for: `1` pending, `2` queued, `3` running,
`6`/finished vs `5` unsatisfactory vs failed — read them off `job_tree`.

---

## 4. Alternative: the `i2run` CLI (when you have the source + a CCP4 env)

If a CCP4i2 checkout and a sourced CCP4 environment are present, `i2run` is the
most robust driver (it *is* how the repo is developed):

```bash
source <ccp4>/bin/ccp4.setup-sh
cd server
env CCP4I2_BACKEND=django DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
    ccp4-python manage.py i2run <task> --project_name <proj> --PARAM value …
```

- **Set `DJANGO_SETTINGS_MODULE` per command — never `export` it.** A shell with
  it exported will make **pytest** flush the live database. Running `i2run` is
  safe; running the test suite against real settings is not.
- Jobs land in the project's database, visible in the desktop app if it points
  at the same home (`~/.ccp4i2x`).

---

## 5. A worked shape: merged data → model

A reliable demo is MR + refinement on a **known-good** dataset (use one from
`server/ccp4i2/demo_data`, e.g. `beta_blip` or `gamma`, that the `i2run` tests
already solve — do **not** improvise a hard case live).

1. `POST /projects/` → project P.
2. `create_task` `phaser_pipeline` (or `phaser_simple`) → job J1.
3. Import the merged MTZ and the search model, and set composition/sequence,
   via `upload_file_param` / `set_parameter` (confirm the object paths from J1's
   container).
4. `GET /jobs/J1/run_time_validation/` — fix anything it flags.
5. `POST /jobs/J1/run/`; poll `job_tree` until finished; read LLG/TFZ from KPIs.
6. `create_task` `refmac` or `servalcat`; bind J1's output coordinates + the
   reflections (see the `sameCrystalAs` note below); run; watch **R / R-free**
   fall in the KPIs.

**Note (a real guardrail you'll meet):** refmac's `XYZIN` and servalcat's `XYZIN`
declare `sameCrystalAs` against the reflection data, so `run_time_validation`
will *warn* (overridable) if the model and data cells don't match and *block* on
an incompatible point group. That is the substrate protecting the run — surface
it, don't paper over it.

Tip: **provenance is the demo.** A handful of finished projects with visible
job trees and dropping R-free is often more convincing than one live run —
and it's what `job_tree` gives you for free.

---

## 6. What you (the agent) must supply today

These are the audited gaps — the cognitive faculties CCP4i2 does not yet own.
Until §7 says otherwise, they are your job:

- **Goal & success test.** CCP4i2 has no "solve this structure" objective. Hold
  the goal and the acceptance criteria yourself (e.g. "R-free < 0.30 and an
  interpretable map").
- **Planning across tasks.** Pipelines decide *within* themselves; nothing plans
  the *whole* path. You sequence the tasks and handle branch points (space-group
  choice, MR solution selection, MR vs experimental phasing).
- **Judging results.** CCP4i2 *reports* quality; it rarely *judges* it. Read the
  KPIs and reports and decide pass / fail / backtrack. Treat a finished job as
  "ran", not "succeeded".

Etiquette: creating projects and jobs is safe and normal. Confirm before
deleting or overwriting. Respect `run_time_validation` — it is the safety net
that makes an external agent trustworthy.

---

## 7. How this evolves (mapped to the anatomy audit)

Each gap gets a "now built" note here as it closes:

| Anatomy faculty | Today | Will be closed by |
|---|---|---|
| Tools / actions | ✅ REST API + `i2run` | (add a curated **MCP facade** so any agent connects with no bespoke glue) |
| Memory & provenance | ✅ recorded & reproducible | (expose it as agent-readable *working* memory, not just a record) |
| Safety & validation | ✅ strong input checks | — |
| Expert knowledge | ◐ locked in pipelines | a **task/recipe catalogue** the agent can consult |
| Perceive / Evaluate & critique | ○ metrics shown, not judged | result-critique helpers that turn KPIs/reports into pass-fail signals |
| Goal & planning | ○ absent | a `JobPlan` resolver (see `docs/NLP_JOB_CONSTRUCTION_DISCUSSION.md`) |

When any row moves from ○/◐ to ✅, update §1, §6 and the relevant flow above so a
fresh agent inherits the new capability instead of re-deriving it.

---

## 8. Related docs

- Strategy & anatomy: `docs/NLP_JOB_CONSTRUCTION_DISCUSSION.md` (§0 substrate stance).
- Task registry & names: `server/ccp4i2/core/tasks.py`.
- Validation discipline: `CLAUDE.md` → "Task Validation".
- Running & tests: `CLAUDE.md` → "Running", "Tests".
