# Run targets: where a job runs, and how a deployment adds a place

CCP4i2 ships one place for a job to run, a local subprocess, and a hook for a
deployment to register another. Nothing in CCP4i2 names a platform. This is
decision 18 of [the PanDDA campaign design note](pandda-campaign-design.md),
under one rule: **CCP4i2 builds generic capabilities with hooks; each
environment registers what it has.** And one invariant every change is
judged by: the desktop developer setup and the built desktop apps must keep
working with no configuration at all.

## The hook

```python
# ccp4i2/config/settings.py — what CCP4i2 ships
CCP4I2_RUN_TARGETS = {"local": "ccp4i2.lib.dispatch.local.LocalTarget"}
CCP4I2_JOB_TARGET = os.environ.get("CCP4I2_JOB_TARGET", "local")

# a deployment's settings overlay — a target CCP4i2 has never heard of
from ccp4i2.config.settings import *
CCP4I2_RUN_TARGETS = {**CCP4I2_RUN_TARGETS,
                      "azure": "azure_extensions.dispatch.ServiceBusTarget"}
CCP4I2_JOB_TARGET = "azure"
```

A target is resolved by name through `ccp4i2.lib.dispatch`: `get_target(name)`
imports the class the first time it is asked for, so a platform SDK is loaded
only where that platform is configured. `available_targets()` reports every
registered name with what it can do, and reports a target that fails to load
with its error rather than hiding it. `CCP4I2_JOB_TARGET` may also come from
the environment; the setting wins when both are present. With neither
present, everything is local.

## Two interfaces, in `ccp4i2/lib/dispatch/base.py`

| Axis | What moves | Methods | Job status while it runs |
|---|---|---|---|
| A, a **job target** | the whole CCP4i2 job (worker needs CCP4i2 and CCP4) | `run_job(job, *, synchronous=False) -> dict` | `QUEUED`, then whatever the worker sets |
| B, a **program target** | one heavy external program (worker needs only that program) | `submit(tree, argv, out_dir, sizing_hint) -> handle`, `poll(handle) -> state`, `cancel(handle)`, `logs(handle) -> path` | `RUNNING_REMOTELY`, holding the handle |

They are `typing.Protocol`s: a class implements either or both by having the
methods, subclassing nothing. `run_job` returns the dict every caller of
`run_job_context_aware` already handles, `{"success": True, "data": job,
"status": 200}` or `{"success": False, "error": ..., "status": ...}`, sets the
job `QUEUED` once handed over, and never raises. A program target never
classifies a failure: `logs` returns the path of the program's stderr under
the job directory, present by the time `poll` is terminal, and CCP4i2's own
failure catalogue reads it. So a new target inherits the whole taxonomy.

## What resolves through it today

Axis A: `lib/utils/jobs/context_run.run_job_context_aware` picks the target
(`local` when `force_local`, else `CCP4I2_JOB_TARGET`) and calls `run_job`.
Every job the API, the interactive-session finish, and the PanDDA fan-out
start goes this way. `program_checks_are_authoritative()` is true only on the
local target with CCP4 present, as before.

Axis B is defined here and consumed by the PanDDA orchestrator's dispatch
mode (design note §14.3, item 16).

## What a deployment does

1. Put the target class in the deployment's own package (for Materia,
   `azure_extensions`).
2. Register it in that package's settings overlay, as above.
3. Nothing else. `EXECUTION_MODE` and `SERVICE_BUS_*` are not read by CCP4i2;
   a target that wants them reads them itself.

## What CCP4i2 refuses

`tests/unit/slim/test_no_platform_imports.py` fails CI if any module under
`ccp4i2/` imports a cloud SDK (`azure`, `boto3`, `google.cloud`,
`kubernetes`). Mounting a deployment app's URLs when that app is installed is
a hook, not an SDK, and stays.
