# Remote Job Execution (ssh / qsub / SLURM) — Implementation Plan

## Goal

Restore the legacy CCP4i2 "Run on server" UX (see
`server/ccp4i2/docs/general/servers.html`) as a **superset** in the new
Django/React stack: a user with a plain desktop install, running a **local**
project, can dispatch an individual job to a named server over ssh/scp — with no
central deployment, the server chosen and credentials entered at run time.

The same machinery also covers HPC/cluster submission (qsub/SLURM), which the
new stack does *not* currently provide at all.

## Why this is a clean fit (not a port of the legacy design)

The legacy design was **client-side paramiko** inside the Qt GUI. We are *not*
porting that. The new stack already separates client → Django → execution, and
already has a pluggable dispatch (`local` vs `azure`). Three things make the
remote case fall out naturally:

1. **The detached `run_job` process is exactly the legacy "controlling ccp4i2
   process."** The legacy doc describes one controller process per server job
   holding the connection. In the new stack `run_job_local` already spawns a
   detached controller ([context_run.py:373](../server/ccp4i2/lib/utils/jobs/context_run.py#L373)).
   We repurpose it as the **ssh shepherd**: it builds the payload, runs the job
   remotely, holds the channel / polls the scheduler, pulls results back, and
   gleans — all while updating `Job.status` in the **local** DB.

2. **Frontend polling and the gleaner need zero changes.** The client polls
   `GET /api/jobs/{id}/report_xml/` etc. on a 5 s interval
   ([api.ts](../client/renderer/api.ts)). Because the shepherd updates the local
   DB just like a local job, the UI Just Works. Results are imported into
   `job.directory` and run through the existing `glean_job_files`
   ([async_db_handler.py](../server/ccp4i2/db/async_db_handler.py)).

3. **The data model already exists** — `ServerJob` and
   `Job.Status.RUNNING_REMOTELY` are dormant in the schema (see
   [models.py:154](../server/ccp4i2/db/models.py#L154),
   [models.py:198](../server/ccp4i2/db/models.py#L198)).

4. **CCP4 10 bundles ccp4i2** (cf. `docs/CCP4_10_CCP4I2_DISTRIBUTION_SWAP.md`),
   so a remote machine with a CCP4 install already has the ccp4i2 code needed to
   execute a task — the historical requirement that "the server has CCP4" now
   also satisfies "the server has the runner."

## The one hard constraint: credentials & where orchestration runs

Legacy held the per-session password in Qt-client memory, never persisted. In
the new stack the React client must POST credentials to Django, which hands them
to the shepherd. This is **safe only when Django is the user's own localhost**
(the desktop case). It is **not** safe against a shared/central Django — that
would route one user's credentials for a third machine through a multi-tenant
host.

**Gate:** runtime-credential ("enter machine + password now") dispatch is
allowed only when the request is local/loopback (or an explicit
`CCP4I2_ALLOW_INTERACTIVE_SSH=1` opt-in for trusted single-user deployments).
Central deployments use **key-based auth** via `ServerJob.key_file_name` only.
This boundary happens to line up exactly with who wants which feature
(laptop-offload vs institutional HPC).

---

## Architecture

```
React "Run on server…" dialog
        │  POST /api/jobs/{id}/run_remote/ { server_group, machine, username,
        │                                    secret?, ccp4Dir, mechanism }
        ▼
JobViewSet.run_remote   ──►  creates ServerJob row, validates cred gate
        │
        ▼
run_job_context_aware(job)  ── sees ServerJob with remote mechanism
        │                       ⇒ always spawn DETACHED local shepherd
        ▼
manage.py run_job -ju UUID  (the shepherd / "controlling process")
        │
        ▼
async_run_job.run_job_async ── branch: ServerJob present?
        │                          ├─ no  → plugin.process()  (today's path)
        │                          └─ yes → run_job_remote_async(job, server_job)
        ▼
RemoteTransport (ssh | ssh_shared | slurm | qsub)
   build payload ─► transfer ─► launch ─► poll/wait ─► retrieve ─► import+glean
```

### Components

#### 1. `RemoteTransport` abstraction (new)
`server/ccp4i2/lib/remote/transport.py`

```python
class RemoteTransport(ABC):
    def connect(self, server_job, secret): ...
    def push(self, local_path, remote_path): ...     # scp/sftp
    def run(self, argv, cwd, env): ...               # exec, returns handle
    def poll(self, handle): ...                       # status of remote job
    def pull(self, remote_path, local_path): ...
    def cleanup(self, remote_path): ...
```
Implementations: `SshTransport` (paramiko — already shipped with CCP4 — or
`asyncssh`), `SshSharedTransport` (no push/pull), `SlurmTransport` /
`QsubTransport` (run = submit, poll = `squeue`/`qstat`).

Selected by `ServerJob.mechanism` ∈ the existing legacy enum
(`ssh`, `ssh_shared`, `qsub_local`, `qsub_remote`, `slurm_remote`, `custom`).

#### 2. Payload (the "jobball")
`server/ccp4i2/lib/remote/payload.py`

Two strategies — start with (A):

- **(A) Self-contained MVP.** Build a scratch slice of the project containing
  only the job dir + its input files (walk the plugin's `inputData` container —
  mirror of how `glean_job_files` walks `outputData`), plus a **freshly migrated
  SQLite** seeded with the `Project`, `Job`, and input `File` rows. Zip it. On
  the remote, set `DATABASE` to that SQLite and run the **unmodified**
  `manage.py run_job -ju UUID`. This reuses the entire execution stack remotely
  with no special "DB-less" mode. Pull back the job dir + the now-populated
  SQLite; import outputs into the local DB and glean. This is the proven legacy
  shape (jobball + DB fragment) minus the Qt coupling.

- **(B) Shared-filesystem.** `mechanism=ssh_shared`: skip the file zip; still
  ship/point at a scratch SQLite (legacy did this too, to avoid NFS lock issues),
  `cd` to the shared project path, run in place, pull back only the DB.

Reuse `export_project_to_zip` ([JobViewSet.py:49](../server/ccp4i2/api/JobViewSet.py#L49))
as the starting point for the zip builder; trim it to a single-job slice.

#### 3. Shepherd branch
`server/ccp4i2/lib/async_run_job.py` — at the top of `run_job_async`, look up
`ServerJob` for the job. If present with a remote mechanism, delegate to
`run_job_remote_async`; otherwise the existing local path. The shepherd sets
`Job.status = RUNNING_REMOTELY`, records the remote PID / scheduler id in
`ServerJob.server_process_id`, and on completion imports + gleans, then sets
`FINISHED`/`FAILED`. Crash-safety: the existing `run_job_safe.sh` wrapper already
marks the job FAILED if the shepherd dies.

#### 4. Server-group configuration (new) — replaces `serverSetup.params.xml`
`ServerGroup` model: `name`, `mechanism`, `machines` (JSON list), `ccp4_dir`,
`temp_dir` (with `$USER` substitution), `allow_user_machines`, `scope`
(user|installation). Admin CRUD via a small viewset. The "Run on server" dialog
reads these to populate dropdowns. (Legacy stored this as XML in
`.CCP4I2/configs/`; we use the DB + an optional import of the legacy XML.)

#### 5. API surface
- `POST /api/jobs/{id}/run_remote/` — body: server params + secret; enforces the
  cred gate; creates `ServerJob`; calls `run_job_context_aware`.
- `GET/POST/PUT /api/server-groups/` — admin config CRUD.
- `POST /api/server-groups/{id}/test/` — legacy "Test server group": ssh in,
  check `ccp4Dir` and temp dir exist. Returns per-machine results.

#### 6. Frontend
- Split/secondary action next to **Run**: "Run on server…" (visible only when
  ≥1 server group is configured — mirrors legacy "option appears when
  configured").
- Dialog: server group → machine dropdown (+ free entry if allowed) → ccp4Dir
  (prefilled) → username → password **or** key file → mechanism.
- Status badge renders `RUNNING_REMOTELY` ("Running remotely"); the existing job
  list already needs only a label for status 7.
- Secret is sent once over the POST and **never** stored client- or server-side.

---

## Phasing

| Phase | Scope | Outcome |
|------|-------|---------|
| **0** | Confirm `ServerJob` is migrated; add `ServerGroup` model + migration + admin API; add status-7 label to the frontend. | Config plumbing, no execution yet. |
| **1 (MVP)** | `mechanism=ssh` (no shared FS): `SshTransport`, payload strategy (A), shepherd branch, results import + glean, `RUNNING_REMOTELY`, localhost cred gate, "Run on server…" dialog. **Final results only** (no live running report). | A laptop user dispatches a job to a Linux server over ssh and gets results back. |
| **2** | `mechanism=ssh_shared` (skip file copy). | Cluster-with-shared-FS users avoid large transfers. |
| **3** | Live running report: shepherd periodically rsyncs `diagnostic.xml` / `report_xml.xml` from the remote job dir so the polling UI shows progress. | Parity with legacy "running report" polling. |
| **4** | `slurm_remote` / `qsub_*`: submit instead of exec, poll `squeue`/`qstat`. Optional qsub option file (`$CCP4/share/ccp4i2/local_setup/qsub_options`). | HPC scheduler submission. |
| **5** | "Test server group" endpoint, key management UI, per-user vs installation config scope, import of legacy `serverSetup.params.xml`. | Admin ergonomics + migration path. |

## Open questions / risks to resolve in Phase 1

1. **Isolated-SQLite run.** Confirm `manage.py run_job` can run against a
   `DATABASE` override pointing at the shipped scratch SQLite, and that
   `migrate` on an empty file produces a runnable schema (or ship a pre-migrated
   template). This is the load-bearing assumption of payload (A).
2. **Code-version skew.** Remote ccp4i2 must be new-stack and close enough to the
   client to run the task. CCP4 10 bundling helps; add a version handshake on
   connect and refuse/ warn on mismatch.
3. **paramiko vs asyncssh.** paramiko ships with CCP4 (zero install, matches
   legacy) but is sync — fine inside the shepherd. asyncssh fits `async_run_job`
   better but adds a dependency. Lean paramiko in a thread executor for MVP.
4. **Secret handling.** Pass via the transport in-memory only; never log; scrub
   from any serialized job representation.

## Testing

- **Unit:** payload builder (correct input-file slice + seeded DB rows);
  transport interface against a fake; cred-gate logic (reject non-loopback +
  runtime password).
- **e2e (opt-in, needs sshd):** ssh to `localhost` as the current user with a
  small fast task (e.g. a freerflag-class job), assert outputs gleaned into the
  local project identically to a local run. Skippable in CI without sshd.
