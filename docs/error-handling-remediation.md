# Error Handling Remediation

Tracking document for the work arising from the error-handling audit of
`server/ccp4i2/wrappers/` and `server/ccp4i2/pipelines/`.

**Baseline:** commit `b3a83af8a` (branch `django`), scanned 2026-08-24;
re-measured at `7c54f4a5b` and again at `59579c932` (3.1.0a28, 2026-08-25).
**Status:** not started, except the [Bibliography](#bibliography) case, which is
closed and kept as a worked example.

Companion documents:

- [docs/pipeline_best_practices.md](pipeline_best_practices.md) — how pipelines
  *should* be written. Several sections need amending once the core fixes below
  land; see [Best-practices amendments](#best-practices-amendments).
- [docs/authoring-a-task.md](authoring-a-task.md) — the task-authoring spine.
- [docs/test-coverage-tracker.md](test-coverage-tracker.md) — per-task test
  coverage, all of it happy-path today.

## The finding in one paragraph

The wrappers do report their failures. The core machinery discards most of those
reports before anyone sees them. A wrapper author who does everything correctly —
detects the failure, calls `appendErrorReport`, returns `FAILED` — still ends up
with a job marked *Finished* and a blank Diagnostics panel. So the ordering of
this work is: fix the core first, because it is currently punishing correct
wrapper code; then the wrapper cleanup has somewhere to land.

**Do not begin the per-pipeline work before C1 and C2 land.** The i2run suite is
green partly *because* failures are being downgraded to warnings and discarded.
Fixing the core will reveal which pipelines are genuinely broken, and that list is
a far better work queue than the table in
[Per-pipeline triage](#part-3--per-pipeline-triage).

Phase 1 will turn some currently-passing i2run tests red. That is the fix
working, not a regression — the affected tests are predictable in advance and
the handling is set out in [The Phase 1 transition](#the-phase-1-transition).

## A second axis — nothing to discard in the first place

The paragraph above describes reports that are produced and then thrown away.
There is a second shape, found on 2026-08-25 and not covered by any metric in
this document, where **no report is produced because nothing goes wrong at the
point where it goes wrong**.

`CDataFile.setFullPath()` decided whether a path was project-relative by
looking for `CCP4_JOBS` or `CCP4_IMPORTED_FILES` anywhere in the string. A file
chosen from *another* CCP4i2 project matches, so its real location was
discarded and `getFullPath()` rebuilt it under the *current* project, where
nothing exists. No exception, no error report, and the corrupted path was
written into `input_params.xml` at parameter-set time. The failure surfaced
much later as xia2 saying `Could not find … in …` — a true statement, about a
path the core had invented, attributed to the wrong layer entirely. It was not
xia2-specific: it cost inputs in any task fed a file from another project.

The same commit turned up the purest instance of the shape. The database-aware
branch it lives in runs only when the plugin carries a `_dbProjectId`, and the
parent chain holds the plugin *weakly* — so if nothing else holds a reference,
`_find_plugin_parent()` returns `None` and the branch **silently does not run**.
No test reached it for exactly that reason.

Three consequences for this plan:

1. **The scan cannot find this class.** A silently-skipped conditional has no
   `except:`, no `[0]`, no `return FAILED`. Every count in the burn-down is a
   count of *visible* mishandling. Treat the scan as a reading order for Part 3,
   never as a measure of how much is left.
2. **It belongs to the guardrails, not the core-defect list.** The defence is
   G2-style negative-path tests over the paths that rewrite user-supplied data
   — and, specifically, a test that a conditional fast-path is actually taken,
   not merely that its outcome looks plausible when it is skipped.
3. **It raises the priority of C7's provenance work.** A message naming a file
   the user never chose should say who chose it. When the core rewrites an
   input, the rewrite is the thing worth reporting.

## Measuring progress

All counts in this document are reproducible:

```bash
python3 scripts/scan_error_handling.py            # summary + per-pipeline table
python3 scripts/scan_error_handling.py --files    # list the offending files
python3 scripts/scan_error_handling.py --json     # machine-readable

python3 scripts/scan_error_handling.py --predict-red-list   # which i2run tests C1 may flip
```

Stdlib only — no CCP4 environment needed. Re-run it when updating this document
and paste the new figures into [Burn-down](#burn-down). The counts are proxies,
not defects: a bare `except:` around a cleanup path is harmless; one wrapping a
parse of the program's primary output is not. Use them to pick a reading order,
not to declare a file guilty.

## Burn-down

Measured at `59579c932` (2026-08-25) unless noted.

| Metric | Baseline (2026-08-24) | Current | Target |
|---|---:|---:|---|
| `bare except:` | 381 (92 files) | 381 (92) | 0 outside a reviewed allowlist |
| `except …: pass` / `continue` | 83 (49 files) | 84 (50) | 0 |
| `.findall(…)[0]` / `.xpath(…)[0]` | 924 (85 files) | 924 (85) | 0 (via M4) |
| `print()` | 917 (146 files) | 913 (142) | 0 (via M6) |
| `logger.*()` | 56 (5 files) | 54 (5) | rising |
| `return FAILED` with no error report nearby | 128 (56 files) | 128 (56) | 0 |
| `processOutputFiles` returning `FAILED`/`UNSATISFACTORY` | 33 + 3 | 35 + 3 | n/a once C1 lands |
| `appendErrorReport` calls passing explicit severity | 7 of 320 | 7 of 326 | n/a once C2 lands |
| Tasks overriding `runTimeValidity` | 6 of 173 | 6 | rising |
| Tasks overriding `postProcessCheck` | 3 of 173 | 3 | n/a once M1/M2 land |
| Tasks declaring `LOG_FAILURES` | 0 of 173 | 1 (`arcimboldo`) | rising |
| Tasks declaring `AUXILIARY_PROGRAMS` | 0 of 173 | 1 (`arcimboldo`) | rising |
| i2run negative-path tests | 0 | 0 | ≥ 12 |
| unit negative-path tests (hand-counted) | 0 | 28 | rising |
| `C1:`/`C2:` xfail markers outstanding | 0 | 0 | 0 (rises during Phase 1, then falls) |
| Pre-existing `@pytest.mark.skip` in i2run | 18 | 18 | falling — see [transition](#the-phase-1-transition) |

The two counts that moved on their own are worth reading correctly. `print()`
fell by four because the PHIL and xia2 refactors of 2026-08-25 deleted
boilerplate, not because anything was remediated. `processOutputFiles` returning
`FAILED` **rose by two**: authors keep writing the correct check into the one
hook that still discards it. Nothing in this document has yet made a wrapper
author's correct code count for anything.

---

# Part 1 — Core machinery

Seven defects, ordered by how much correct wrapper code each one currently
defeats. Each discards information a wrapper has already gone to the trouble of
producing.

## C1 — `processOutputFiles()` returning `FAILED` does not fail the job

- [ ] Honour the int return in `process()`
- [ ] Delete the duplicate post-process path
- [ ] Triage the i2run failures this exposes

**Where:** `server/ccp4i2/core/CCP4PluginScript.py:853`,
`server/ccp4i2/core/base_object/error_reporting.py:71`

`process()` handles the legacy int-vs-`CErrorReport` dual return convention in
three of its four hooks. It forgets in the fourth:

```python
# process(), the synchronous path that actually runs
error = self.processOutputFiles()
if error:                          # FAILED == 1, so this is True
    self.errorReport.extend(error) # extend(1) -- silent no-op, see below
    print("Warning: processOutputFiles() returned errors "
          "but process completed successfully")
    # status = self.FAILED  <- commented out in the source

# CErrorReport.extend
def extend(self, other):
    if isinstance(other, CErrorReport):   # an int isn't; nothing happens
        self._errors.extend(other._errors)
```

**Blast radius:** 33 of the 91 `processOutputFiles` implementations return
`FAILED` this way — including `refmac`, `servalcat`, `aimless`, `modelcraft`,
`buster`, `acedrg`, `xia2_dials`, `morda_i2`. A further three return
`UNSATISFACTORY`, likewise discarded, so that status is effectively unreachable
from a wrapper. Get the current list with
`python3 scripts/scan_error_handling.py --files`.

The checks being thrown away are exactly the ones the framework cannot do
generically: *"refmac exited 0 but wrote no XMLOUT"*, *"the ligand had no
dictionary"*, *"phaser found no solutions"*. The exit-code case is caught
independently by `postProcessCheck`, which is why this has stayed hidden.

**Related:** `postProcess()` at `CCP4PluginScript.py:2138` contains a *correct*
version of this logic (`status = self.FAILED`), but it is only reachable via the
async path, and `lib/async_run_job.py:281` forces `plugin.doAsync = False`. Two
post-processing paths have diverged and only the broken one runs. Fix the live
one and delete the other rather than repairing both.

**Expect this to turn currently-green i2run tests red.** That is the point.
Capture the list — it is the work queue for Part 3.

## C2 — Error severity is decided by searching the message for the word "exception"

**Split into two steps — see [Why C2 must be staged](#why-c2-must-be-staged).
C2a is safe to land whole; C2b is a per-wrapper migration, not a core commit.**

**C2a — honour what is already declared (safe, land with Phase 1)**

- [ ] Use `ERROR_CODES[code]['description']` as the message when no details are given
- [ ] Honour an explicitly declared `severity` where one exists (3 files today)
- [ ] Leave the default for undeclared codes at `SEVERITY_WARNING` for now
- [ ] Add an import-time check that every code used is declared (20 are not)

**C2b — change the default (per-wrapper, Phase 3)**

- [ ] Add an explicit `severity` to every entry in each wrapper's `ERROR_CODES`
- [ ] Once a wrapper is fully annotated, flip its undeclared-code default to `SEVERITY_ERROR`
- [ ] Retire the substring heuristic when the last wrapper is annotated

**Where:** `server/ccp4i2/core/CCP4PluginScript.py:4462`

```python
def appendErrorReport(self, code=0, details='', ..., severity=None):
    if severity is None:
        if 'exception' in details.lower():
            severity = SEVERITY_ERROR      # 4
        else:
            severity = SEVERITY_WARNING    # 2 -- everything else
```

313 of 320 `appendErrorReport` calls pass no severity, so their level is set by
this substring test. In practice `self.appendErrorReport(201, 'Exit status: 3')`
in the refmac wrapper is filed as a *warning*.

The same method ignores the `ERROR_CODES` dictionaries that **73 wrapper files**
define. Those carry both the human-readable description and, in some cases, a
deliberate per-code severity:

```python
# wrappers/refmac/script/refmac.py -- declared, never read
ERROR_CODES = { 201: {'description': 'Refmac returned with non zero status'},
                203: {'description': 'New library created',
                      'severity': CCP4ErrorHandling.SEVERITY_WARNING},
                204: {'description': 'Program completed without generating XMLOUT.'} }
```

38 calls pass a code and nothing else, on the assumption the lookup happens.
Those produce an error report with an empty `details` field — the user sees a
bare code number. 20 calls use a code their own file does not declare.

> **This defect also disables the project's own safety net.**
> `server/ccp4i2/tests/i2run/utils.py:226` asserts that `diagnostic.xml` contains
> no entry with severity ≥ 4, and explicitly filters warnings out first. Because
> nearly every wrapper-reported failure lands at severity 2, that assertion almost
> never fires. The suite prints `Note: 1 warning(s) in diagnostic.xml (non-fatal)`
> and passes.

### Why C2 must be staged

Two measurements taken while planning Phase 1 changed the shape of this fix.

**`appendErrorReport` severity does not affect job status anywhere.** The only
place a severity gates the outcome is `process()` line 688, for
`runTimeValidity` — and no `validity()` or `runTimeValidity()` override in the
tree uses `appendErrorReport`; they all build their report with
`error.append()` directly. The `maxSeverity()` check in
`lib/async_run_job.py:81` only logs. So C2 governs *display* in the Diagnostics
panel and *sensitivity* of the i2run harness assertion — nothing else. C1 is
what governs job status. The two are independent, and C1 is the correctness fix.

**Defaulting undeclared codes to ERROR would flip 268 of 320 calls at once:**

| Transition | Calls | Note |
|---|---:|---|
| WARNING → ERROR (code declared in `ERROR_CODES`) | 230 | declared, but with a *description* only — no severity key |
| WARNING → ERROR (code not declared at all) | 38 | would be promoted by the proposed default |
| ERROR → ERROR (details contain "exception") | 33 | unchanged |
| non-literal code | 12 | needs manual review |
| explicit `severity=` passed | 7 | unaffected |

Only 3 files declare a per-code severity at all, so "declared" almost always
means "declared with a description and nothing else". Promoting all 268 in one
commit would fail nearly every i2run test that emits any error report, for
reasons unrelated to whether the job worked — and would bury the C1 signal,
which is the one worth reading.

The safe order is therefore: land C2a (pure gain — 38 blank `details` fields get
their text back, and the 3 deliberate severities start being honoured, with no
severity promotion at all), then annotate `ERROR_CODES` wrapper by wrapper in
Phase 3, where the person doing the work knows whether code 201 in *that*
wrapper is fatal. The test-harness guard becomes meaningful progressively rather
than all at once.

Regenerate the table above with
`python3 scripts/scan_error_handling.py --json`.

## C3 — Exceptions in pipeline continuation callbacks are swallowed into the server log

- [ ] Record slot exceptions on the emitting plugin's error report
- [ ] Mark the emitting plugin failed before continuing to other slots
- [ ] Keep the inter-slot isolation

**Where:** `server/ccp4i2/core/base_object/signal_system.py:355`

Pipelines chain their steps by connecting handlers to a sub-plugin's `finished`
signal. `Signal.emit` wraps every slot invocation:

```python
try:
    result = slot(*args, **kwargs)
    results.append(result)
except Exception as e:
    logger.error(f"Slot execution error in {self._name}: {e}")
    logger.error(f"Traceback:\n{traceback.format_exc()}")
    # Continue with other slots even if one fails
```

Because the top-level plugin runs with `doAsync = False`, the whole pipeline body
executes inside these slots. Any unhandled Python error in a continuation — the
`IndexError` from a `.findall(…)[0]` against XML a failed program never wrote, a
missing attribute, a bad path — is caught here. It never reaches
`plugin.errorReport`, never reaches `diagnostic.xml`, and never changes the job
status. The pipeline simply stops partway through and reports whatever status was
last emitted.

This is the mechanism behind the class of report that arrives as *"the job says
it finished but there's nothing there"*. The traceback exists, in the server log,
disconnected from the job.

The generic-exception capture already written in `lib/async_run_job.py:155` is the
right shape; it needs to also apply one level down, inside `emit`.

## C4 — `reportStatus()` does not stop anything, but authors use it as if it did

- [ ] Add `self.fail(code, details='')` raising `CPluginFailure`
- [ ] Catch it at the top of `process()` → clean FAILED with the report intact
- [ ] Migrate the `reportStatus`-as-abort sites
- [ ] Guard against `reportStatus` being called twice on one plugin

**Where:** `server/ccp4i2/core/CCP4PluginScript.py:2254`;
specimens at `pipelines/phaser_pipeline/script/phaser_pipeline.py:292,302` and
`pipelines/prosmart_refmac/script/prosmart_refmac.py:419`

`reportStatus()` saves params, updates the database row, emits `finished`, and
returns. It is not terminal. But its name reads like one, and it is called in
positions where the author clearly expected execution to end. The clearest
specimen is 8 lines long and contains four instances of the pattern:

```python
def checkSolutionsFound(self, finishStatus, failedErrCode):
    if finishStatus == CPluginScript.FAILED:
        self.appendErrorReport(failedErrCode)
        self.reportStatus(finishStatus)      # does not return -- falls through
    self.appendXML(self.phaserPlugin.makeFileName('PROGRAMXML'),
                   'PhaserMrResults')        # runs even after "failure"
    if self.xmlroot.xpath('//solutionsFound')[0].text == 'False':
        #                 ^ IndexError when phaser wrote no XML
        self.reportStatus(CPluginScript.UNSATISFACTORY)   # also falls through
```

A failed phaser run proceeds to parse XML that doesn't exist, raises
`IndexError`, and that exception is then eaten by **C3**. The user gets a job
stuck in whatever state the first `reportStatus` left it. The same shape appears
in `prosmart_refmac.checkFinishStatus()`, whose callers continue straight into a
`try` block after the "abort".

This is the clearest single example of the discoverability problem: there was no
obvious way to say *"stop, this failed"*, so people invented one that doesn't
work. `self.fail()` is the fix, and is the addition most likely to change how new
wrappers get written.

## C5 — The severity threshold that fails a job is inverted between hooks

- [ ] One rule: `maxSeverity() >= SEVERITY_ERROR` fails; below that is collected and shown
- [ ] Document it in one sentence at the top of `process()`
- [ ] Amend `pipeline_best_practices.md` §5.2

**Where:** `server/ccp4i2/core/CCP4PluginScript.py:667-845`

`CErrorReport.__bool__` is `maxSeverity() >= WARNING`. Within one method,
`process()` applies two different rules to that fact:

| Hook | Test applied | A returned WARNING… | A returned ERROR… |
|---|---|---|---|
| `runTimeValidity` | `if error.maxSeverity() >= 4` | passes | fails the job |
| `processInputFiles` | `elif result:` | **fails the job** | fails the job |
| `makeCommandAndScript` | `if error:` | **fails the job** | fails the job |
| `startProcess` | `elif result:` | **fails the job** | fails the job |
| `processOutputFiles` | — | ignored | ignored (C1) |

Net effect: a wrapper that appends a benign advisory while building its command
line kills the job; a wrapper that detects a genuine failure while reading its
outputs does not. Authors have no way to reason about which is which.

## C6 — The Diagnostics panel reads three field names the backend never writes

- [ ] Settle the `diagnostic.xml` element list in one place
- [ ] Emit `className`, `description`, `stack`, and the `ccp4i2_header` block
- [ ] Surface `<name>` in the UI
- [ ] Contract test at both ends

**Where:** `server/ccp4i2/core/base_object/error_reporting.py:224`,
`client/renderer/components/diagnostic.tsx:86`

`CErrorReport.getEtree()` emits `<class> <code> <details> <name> <severity>
<severityName>`. The React component asks for:

```ts
className:   getTextContent("className")    // never written -> ""
description: getTextContent("description")  // never written -> ""
stack:       getTextContent("stack")        // never written -> ""
severity:    getTextContent("severityName") || getTextContent("severity")  // ok
details:     getTextContent("details")                                     // ok
```

Every accordion header renders as ` - 201`, and the field labelled "Description"
is blank. The `<name>` element — which carries *which parameter* the error
concerns, and is what the validation system depends on for field-level
highlighting — is written by the backend and read by nothing. The component also
renders a header card from a `<ccp4i2_header>` element that
`write_diagnostic_xml()` does not produce, so the provenance block never appears.

Drift bugs recur, so the fix includes a test at each end: a Python test asserting
`getEtree()` emits the fields, and a TypeScript test asserting the component
queries only those fields.

## C7 — stderr goes to the server console; no execution has a timeout

- [ ] Attach the stderr tail (and log tail) to the failure report
- [ ] Wire a real timeout to the handler that already exists for it

**Where:** `server/ccp4i2/core/CCP4PluginScript.py:1839,1878`

When a program exits non-zero, `_startProcessSync` reads its stderr and then
`print()`s it — to the server's stdout, not into the error report. The user is
shown `Process refmac5 exited with code 1` with no further detail, while the
actual message sits in a terminal they cannot see.

**Appending the last ~40 lines of stderr into `details` is a few lines of code
and is probably the highest ratio of user-visible improvement to effort in this
entire document.**

Separately, `subprocess.run()` is called with no `timeout=`, which makes the
`except subprocess.TimeoutExpired` branch immediately below it unreachable. A
hung program hangs the worker indefinitely with no diagnostic. A generous default
(configurable, per-task overridable) raising a clear "timed out after N minutes"
would close that.

---

# Part 2 — Mechanisms

The audience is a distributed set of scientist-authors who will not read a style
guide and will copy whichever nearby wrapper looks closest to their problem.
Anything depending on authors *knowing* a convention will decay. What survives is
defaults in the base class, and CI that refuses the bad shape.

## M1 — A declared output contract, checked by the base class

- [ ] Base `postProcessCheck` walks `outputData` declarations from the def.xml
- [ ] Fail with a specific message for any non-optional output missing or zero-length
- [ ] Mark genuinely-optional outputs in the def.xml, not in code
- [ ] Delete the hand-rolled equivalents from wrappers

Verifying that a program actually produced its outputs is left to each wrapper's
`processOutputFiles`, done inconsistently, and then discarded (C1). The def.xml
already declares every `outputData` file. Roughly 90 wrappers get a real
completeness check for free — *"refmac completed but produced no FPHIOUT; see log
line 4213"*.

## M2 — Declarative log scanning

- [x] `LOG_FAILURES` class attribute honoured by base `postProcessCheck`
- [ ] Harvest the strings from existing inline `if … in logFileText` checks

**Landed** (arcimboldo case, 2026-08-25). `CPluginScript.LOG_FAILURES`,
`scanLogForFailures()` and the `reportLogFailures()` surfacing hook are in
`core/CCP4PluginScript.py`; `tests/unit/plugins/test_log_failures.py` covers
them. Deliberately placed in `postProcessCheck`, not `processOutputFiles`, so it
works on the synchronous and asynchronous paths alike and does not depend on C1.
Tasks declaring no patterns are untouched — the scan is skipped entirely.

Only 3 tasks override `postProcessCheck`, so the very common CCP4 case of a
program exiting 0 while writing a fatal message to its log is essentially
unhandled. Rather than asking 173 task authors to write log parsers, let a
wrapper declare patterns as class data:

```python
LOG_FAILURES = [
    (r'Your coordinate file has a ligand which has either minimum or no description',
     301, 'No geometry dictionary for a ligand -- run Make Ligand first'),
    (r'\*\*\* Fatal error', 302, None),
]
```

The base `postProcessCheck` applies them after a successful exit. This is a shape
scientists reliably get right — it is a table, not control flow.

## M3 — `self.fail()`

Covered by C4. One import, one call, unambiguous semantics.

## M4 — Safe accessors for XML scraping

- [ ] `self.xmlValue(root, xpath, default=None)`
- [ ] `self.xmlRequired(root, xpath, code, details)` → raises `CPluginFailure`
- [ ] Worked example in `pipeline_best_practices.md`
- [ ] Burn down the 924 unguarded indexes

924 instances of `.findall(…)[0]` / `.xpath(…)[0]` is not carelessness so much as
the absence of an alternative. With these two helpers, the distinction between
"this is optional" and "this being absent means the run failed" becomes something
the author states, instead of something decided by an `IndexError`.

`pipelines/tableone` — 25 unguarded indexes in 304 lines, pure XML aggregation —
is the natural first conversion and worked example.

## M5 — Report classes that degrade section by section

- [ ] Guard each top-level report section in the report base class
- [ ] Emit a visible "this section could not be generated: <reason>" block
- [ ] Land the failure in the error report at WARNING
- [ ] Say so in the UI when a served report is a preserved rendering rather than
      a current one

141 of the 173 tasks have a report class, and reports are where the unguarded XML
indexing concentrates. A report that raises produces no report at all — the user
loses the 90% that would have rendered because of one absent table.

**Amended 2026-08-25.** `get_job_report_xml()` no longer deletes a cached
`report_xml.xml` when regeneration produces an error report; it keeps and serves
the cached copy. That is right — a report class is written against the
`program.xml` its wrapper produced at the time, so a newer report class can fail
on an older job's output through no fault of the job, and only renderings that
worked are ever cached. It is also why the cache is deliberately not
version-stamped: invalidating on a report-format version would regenerate every
finished job against report classes that may no longer understand their
`program.xml`.

But as it stands the regeneration failure is recorded in a `logger.warning` that
no user will ever see, and the viewer is shown a rendering from an earlier
report class with nothing saying so. That is a loud, destructive failure traded
for a quiet one — the thing this document is against. M5 is therefore not
"guard the sections" alone: **a preserved rendering must announce itself**, with
the reason regeneration currently fails. The per-section guards then reduce how
often the whole-report fallback is needed at all.

## M6 — Logging that reaches the job, not the console

- [ ] `self.log(msg, level)` writing to the job's log *and* the Python logger
- [ ] Lint `print()` in `wrappers/` and `pipelines/` to an error

917 `print()` calls against 56 `logger` calls in 5 files. `print()` is not merely
untidy here: per `CLAUDE.md`, a non-ASCII `print()` raises `UnicodeEncodeError`
on Windows, and inside one of the 381 bare `except:` blocks that becomes a silent
job failure — which is exactly how the aimless_pipe Windows bug behaved.

## M7 — Make `runTimeValidity` discoverable

> Partly addressed 2026-08-25: the base `runTimeValidity` program check now
> blocks rather than merely warns wherever it is authoritative, so authors get
> a visible reason to reach for the hook. See the changelog.


- [ ] Ship a stubbed `runTimeValidity` in the task-authoring template
- [ ] Name the three questions worth asking in its comment

Six tasks use it. The mechanism is good — `_checkProgramAvailable` and
`_checkSameCrystalAs` in the base class are exactly the right idea, well
executed, and `_checkProgramAvailable` now covers `AUXILIARY_PROGRAMS` too (M8). The gap is that nothing prompts an author to consider it. Put the hook
in front of them at the moment they are writing, with the three questions: *are
the inputs mutually consistent? is everything the program needs present? is this
combination one the program supports?*

## M8 — Declare the tool; let the base class find it, and say so when it is missing

- [x] `AUXILIARY_PROGRAMS` — binaries a task drives from *inside* `TASKCOMMAND`
- [x] `PHIL_SCOPE` / `PHIL_PARAMS_FILE` / `PHIL_PROGRAM` — where a PHIL task's
      parameters come from
- [ ] State the missing-dependency policy once, in `authoring-a-task.md`
- [ ] Apply the same shape to monomer dictionaries and reference data

Two changes on 2026-08-25 landed the same idea from different directions, and
together they are a mechanism in their own right rather than incidents of M2.

`AUXILIARY_PROGRAMS` exists because `TASKCOMMAND` names only the leaf binary:
arcimboldo declares `ARCIMBOLDO_LITE` but drives phaser and shelxe from within
it, so the pre-run availability check saw nothing to check. Declaring the extra
binaries as class data brings them under the same check and onto Preferences →
Program locations, and made a correctly-configured `SHELXDIR` effective for that
task for the first time.

`PhilPluginScript` did the same for parameters. Seven wrappers each implemented
`get_master_phil()`, all seven the same import-and-return in one of three
shapes, the bodies carrying no information beyond the name of the thing to
import. They now declare the shape and the base class resolves it — and in doing
so fixed a disagreement of exactly the audited kind: the phasertng pair caught
`ImportError` and returned `None`, while the xia2 tasks let it escape, so a
missing tool stopped the task **loading at all** rather than opening so it could
say what it needed.

That disagreement is the general point, and the policy is worth stating once for
every mechanism in Part 2 to follow:

> A missing external dependency resolves to a defined "absent" value, logs one
> warning naming the task and the declaration, and lets the task open. The task
> tells the user what it needs. It does not fail to load, and it does not
> pretend the dependency is there.

Declaring nothing, or declaring two sources, is an error at the point of
declaration rather than a silent preference for one — the same reasoning as
M1's output contract: an author who states an intention gets a specific message,
and an author who states none gets told to.

---

---

# Part 3 — Per-pipeline triage

Measured per pipeline directory by `scripts/scan_error_handling.py`. The ranking
is a reading order, not a severity list — and it is superseded by the red list
from C1/C2 once those land.

Legend: **swall** = `except …:` immediately followed by `pass`/`continue`;
**[0]idx** = `.findall(…)[0]` or `.xpath(…)[0]` without a length guard;
**retFAIL** = `return self.FAILED` / `return CPluginScript.FAILED`. Counts include
each pipeline's nested wrappers and report modules.

| Pipeline | LOC | bare | swall | print | [0]idx | retFAIL | Done | First move |
|---|---:|---:|---:|---:|---:|---:|:--:|---|
| `phaser_pipeline` | 3126 | **33** | **9** | 23 | **48** | **62** | [ ] | **Read first.** `checkFinishStatus` and `checkSolutionsFound` both fall through after "aborting" (C4), straight into unguarded xpath indexing (C3). Highest concentration of every systemic defect. |
| `prosmart_refmac` | 1404 | **20** | 6 | 5 | **45** | 4 | [ ] | Same `checkFinishStatus` fall-through. Several `except Exception: pass` around dictionary merge and molprobity hide genuine failures behind "optional" comments — make them WARNING-level reports instead of silence. |
| `crank2` | 16159 | 12 | 3 | 8 | 0 | 3 | [ ] | Largest by far, but sparse in the risky idioms — mostly vendored upstream code. Scope a boundary: harden the CCP4i2-side adapter, leave the interior alone. |
| `import_merged` | 3676 | 4 | 1 | **81** | **51** | 0 | [ ] | Import paths see the widest variety of malformed user files, and this one reports nothing back (`retFAIL` = 0). Convert the `print()` diagnostics into error-report entries so a bad mmCIF produces an actionable message. |
| `aimless_pipe` | 3123 | 9 | 3 | **51** | **57** | 1 | [ ] | Known Windows-fragile (the `print()`/`UnicodeEncodeError` case in `CLAUDE.md`). Nine bare `except:` around scaling-statistics parsing; audit which can hide a genuinely failed scaling run. |
| `servalcat_pipe` | 1707 | 10 | 8 | 29 | 20 | 0 | [ ] | Good `appendErrorReport` density (33) but zero `FAILED` returns — it reports and continues. Once C2 lands those become real errors; check each is at the severity intended. |
| `dr_mr_modelbuild_pipeline` | 1507 | 11 | 1 | **81** | 29 | 6 | [ ] | One of four pipelines using `doAsync`/`connectSignal` chaining, so fully exposed to C3. Re-check step transitions once slot exceptions surface. |
| `MakeLink` | 564 | 6 | 0 | 46 | 3 | **21** | [ ] | 21 `return FAILED` and no `appendErrorReport` anywhere — every failure path is silent by construction. Small file, high yield: give each return a code and a message. |
| `SubstituteLigand` | 945 | 0 | 3 | 23 | 0 | 0 | [ ] | Cleanest of the large pipelines — 33 error reports, no bare excepts. Already the reference implementation in `pipeline_best_practices.md`; keep it that way. |
| `MakeProjectsAndDoLigandPipeline` | 318 | **10** | 3 | 16 | 18 | 1 | [ ] | Very high defect density for its size — 10 bare excepts in 318 lines. Cheap to rewrite outright rather than patch. |
| `tableone` | 304 | 0 | 0 | 1 | **25** | 0 | [ ] | Pure XML aggregation; 25 unguarded indexes in 304 lines. Prime candidate for M4 as a worked example. |
| `import_xia2` | 574 | 6 | 0 | 2 | 6 | 0 | [ ] | Import path; same argument as `import_merged` but smaller. Reports 11 errors, returns FAILED never. |
| `phaser_ep` | 318 | 0 | 0 | 0 | 5 | 7 | [ ] | Handle with `phaser_pipeline` — shared idioms, shared fixes. |
| `phaser_rnp_pipeline` | 318 | 2 | 1 | 2 | 3 | 3 | [ ] | As above. |
| `phaser_simple` | 125 | 0 | 0 | 0 | 0 | 0 | [ ] | As above. |
| `molrep_pipe` | 513 | 5 | 0 | 2 | 3 | 0 | [ ] | Low density. Sweep after the core fixes land. |
| `PrepareDeposit` | 318 | 1 | 0 | 6 | 3 | 1 | [ ] | Low density. Sweep after the core fixes land. |
| `pisapipe` | 305 | 2 | 0 | 1 | 1 | 4 | [ ] | Low density. Sweep after the core fixes land. |
| `LidiaAcedrgNew` | 277 | 1 | 0 | 6 | 0 | 1 | [ ] | Low density. Sweep after the core fixes land. |
| `import_serial_pipe` | 183 | 0 | 0 | 4 | 1 | 0 | [ ] | Low density. Sweep after the core fixes land. |
| `shelx` | 40 | 0 | 0 | 0 | 0 | 0 | [ ] | Nothing to do. |

## Wrappers

The same defects apply to the 127 wrapper directories, but there is no useful
per-wrapper ranking until C1 lands and the i2run suite says which ones actually
break. The known list to start from is the 33 wrappers whose
`processOutputFiles` returns `FAILED` — `python3 scripts/scan_error_handling.py
--files` prints it.

---

# Part 4 — Guardrails

## G1 — A conformance test over the task registry

- [ ] Parametrised test over `TASKS` in `server/ccp4i2/core/tasks.py`

`core/tasks.py` gives a single enumerable list of 173 tasks. A parametrised test
can assert, without running any program, that each task: loads its def.xml;
instantiates; declares `ERROR_CODES` covering every code its module passes to
`appendErrorReport`; has a `TASKCOMMAND` that resolves or is explicitly marked as
a pipeline; and (once M1 lands) declares which outputs are optional. Fast, needs
no CCP4 binaries, runs on every commit alongside the existing unit tier.

## G2 — Negative-path tests

- [ ] Shared fixture set: truncated MTZ, PDB with undescribed ligand,
      non-matching sequence, deliberately absent binary
- [ ] Run against the dozen most-used tasks
- [ ] Assert `diagnostic.xml` contains an ERROR whose `details` is non-empty and
      names the actual problem
- [ ] Cover the paths that **rewrite** user-supplied data — a file from another
      project, a file in a sub-directory of `CCP4_IMPORTED_FILES` — asserting
      the recorded path, not merely that the job ran

All 96 i2run tests are happy-path; the failure modes users report are the unhappy
ones. That last-but-one assertion is the one that prevents regression to blank
error cards.

The last bullet comes from [A second
axis](#a-second-axis--nothing-to-discard-in-the-first-place). Where a fast path
is conditional, assert that it was *taken*: the `setFullPath()` bug hid behind a
branch that no test reached because the plugin was held weakly and the branch
therefore silently did not run. A test whose fixture drops the reference passes
while testing nothing.

## G3 — Lint rules scoped to `wrappers/` and `pipelines/`

- [ ] `E722` (bare `except:`) → error, with a short reviewed allowlist
- [ ] `except …: pass` → error; require at least a logged reason
- [ ] `print()` → error; direct to `self.log()`
- [ ] Non-ASCII literals inside `print()` → error (the Windows
      `UnicodeEncodeError` class documented in `CLAUDE.md`)
- [ ] Custom check for `.findall(…)[0]` / `.xpath(…)[0]` → warning, pointing at M4

Introduce these with a baseline file so existing code doesn't block anyone, and
burn the baseline down as the Part 3 work proceeds. New code is held to the rule
from day one; that alone stops the numbers growing.

## G4 — A schema contract for `diagnostic.xml`

Covered by C6.

---

# Sequence

### Phase 0 — Capture the baseline

Before touching anything, record what the suite does today, so the effect of
Phase 1 is a diff rather than an impression:

```bash
CCP4_SETUP=/path/to/ccp4-XXXX/bin/ccp4.setup-sh \
CCP4_LABEL=pre-C1 \
bash server/run_i2run_baseline.sh
# -> server/.test-baselines/pre-C1/{results.xml,pytest.log}
```

That runner already exists (it was written for the CCP4 version migration) and
writes JUnit XML, which diffs test-by-test. Re-run it with `CCP4_LABEL=post-C1`
after Phase 1 and compare the two `results.xml`.

### Phase 1 — Stop discarding what already exists

C1, **C2a** (not C2b), C5, C7, C6. Small, contained, entirely within the core.

See [The Phase 1 transition](#the-phase-1-transition) for what "going red" means
here and how to handle it — it is the part of this plan most likely to be
misread as a regression.

### Phase 2 — Give authors somewhere to land

C3, C4/M3, M1, M4, M5, G1, G3. Everything here is additive — no wrapper has to
change for it to start paying off, which is what makes it stick.

### Phase 3 — Work the pipelines

Start from the Phase 1 red list, not from the Part 3 table. **C2b** rides along
here: as each wrapper is worked, annotate its `ERROR_CODES` with explicit
severities and flip that wrapper's default. `phaser_pipeline` and
`prosmart_refmac` first (the fall-through pattern), then `MakeLink` (21 silent
failure returns, small file), then the import paths (widest input variety). M2 log
tables harvested from existing inline checks. G2 negative-path fixtures. Burn the
G3 lint baseline down as each pipeline is touched, rather than as a separate
campaign.

---

# The Phase 1 transition

Fixing C1 will turn some currently-passing i2run tests red. This section exists
because that will look like a regression and is not one, and because "we'll see
what breaks" is not a plan.

## Nothing breaks; C1 reveals

A test that flips red under C1 was already running a wrapper whose own
`processOutputFiles()` had concluded the run failed. The conclusion was being
discarded, so the job was recorded as *Finished* and the test passed. After C1
the same run produces the same conclusion and the job is recorded as *Failed*.

The software's behaviour has not got worse. The **reporting** has got truthful,
and the test is now asserting the thing it always meant to assert. Each red test
is one of:

1. **A latent bug**, now visible — the run really was failing. Fix the wrapper
   or the pipeline. This is the valuable outcome and the reason to do the work.
2. **An over-strict wrapper check** — the wrapper declares failure on a
   condition that is actually acceptable (an optional output absent, a log
   phrase that is not fatal). Fix the wrapper's check, or downgrade it to a
   warning. Also valuable: that check was wrong and nobody could tell.
3. **A test asserting too little** — the test only checked the job ran, and the
   job never produced a usable result. Strengthen the assertion.

None of the three is "C1 broke it".

## The list is predictable in advance

You do not have to discover it empirically:

```bash
python3 scripts/scan_error_handling.py --predict-red-list
```

This walks the `makePluginObject()` subjob graph, so it catches wrappers reached
through a pipeline as well as those invoked directly. As of the baseline:

- **27 active i2run test files** run a task that reaches a C1-affected wrapper
- **4 already-skipped files** would also be affected if re-enabled
- **52 test files** run only unaffected tasks

**At risk is not the same as will fail.** A wrapper's failure branch fires only
on a specific condition — no XMLOUT written, an undescribed ligand, no solution
found — that a happy-path test usually does not hit. The genuinely red set will
be a subset, probably a small one. What matters is that the boundary is known
before the work starts, so an unexpected red *outside* these 27 files is
immediately identifiable as a real regression in the C1 change itself.

## Mechanism: strict xfail, not skip

When a red test is triaged as category 1 or 2 above and cannot be fixed in the
same commit, mark it:

```python
@pytest.mark.xfail(
    strict=True,
    reason="C1: aimless reports no XMLOUT on this dataset; "
           "see docs/error-handling-remediation.md#part-3--per-pipeline-triage",
)
```

`strict=True` is the whole point:

- an expected failure keeps CI green, so the suite stays usable during the migration;
- an *unexpected pass* is reported as a failure, so when someone fixes the
  underlying wrapper the marker cannot be silently left behind — the suite forces
  it to be removed;
- the set of markers is a live, greppable backlog:
  `grep -rn "C1:" server/ccp4i2/tests/i2run/ | wc -l` is the burn-down.

**Do not reach for `@pytest.mark.skip`.** The suite already carries 18 skips,
several of which record real unfixed bugs (`test_ctruncate.py`,
`test_merge_mtz.py`, `test_editbfac.py`). A skip hides the test forever and
reports nothing when the problem is fixed; that is exactly the failure mode this
whole document is about, reproduced in the test suite. Every marker added during
this migration should carry a `C1:`/`C2:` prefix and a link to the relevant
section here.

## Reading the diff

```bash
CCP4_LABEL=post-C1 bash server/run_i2run_baseline.sh
# then compare server/.test-baselines/{pre-C1,post-C1}/results.xml
```

Three questions for each newly-failing test, in order:

1. Is it in the 27-file predicted list? If not, suspect the C1 change itself
   before suspecting the wrapper.
2. Does `diagnostic.xml` now carry an error whose `details` names a real
   problem? If the details are empty, that is a C2a bug, not a wrapper bug.
3. Did the job produce usable output despite the reported failure? If yes, the
   wrapper's check is wrong (category 2); if no, it was right all along
   (category 1).

## What to tell people

The honest headline is *"the test suite got more sensitive, and here is what it
found"* — not *"we broke 27 tests"*. The number that should go in the release
notes is the count of category-1 findings: latent failures that were shipping to
users as successful jobs.

---

# Bibliography

<a id="bibliography"></a>

A worked instance of this document's whole thesis, found the day the CI job in
PR #276 first ran and closed two PRs later. Kept here because it is the shortest
complete illustration of the argument: a silent failure, a fix that went the
wrong way, and a guard that caught it.

## What happened

Two callers read `references/{X}.medline.txt` and build `X` differently:

| Caller | Builds the name from | For AceDRG wants |
|---|---|---|
| `report/metadata.py` `loadFromMedLine` | report class `TASKNAME` | `Acedrg.medline.txt` |
| `lib/utils/jobs/bibliography.py` | citation key / alias | `acedrg.medline.txt` |

On a case-insensitive filesystem both resolve, so the two conventions were never
forced to agree. On Linux only one can, whichever way the file is spelled.
Renaming the file fixed one caller and broke the other; the landed fix
(PR #276) resolves case-insensitively via `CCP4Utils.findReferenceFile` and
both callers use it.

**The reason it went undetected for as long as the file existed** is the part
that belongs in this document: `loadFromMedLine` appended to `self.errReport`
and returned, and *nothing reads that errReport*. The report rendered normally,
just without the citation. That is C3/M5 exactly — a failure recorded into a
channel with no consumer.

## Status — closed

- [x] Filename resolution, case-insensitively rather than by renaming one file
      — PR #276 (`c2b485386`)
- [x] A missing reference file warns, naming the task; the non-citable allowlist
      moved to import-free `core/citations.py` so `report/metadata.py` stays
      Django-free
- [x] Report side routed through one citation resolver — `ReferenceGroup.loadFromTask()`
      — PR #277 (`7c54f4a5b`), which also moved the alias map to
      `core.citations.TASK_CITES` and refreshed 12 expected-report fixtures

## What the fix was worth

Of the 144 report classes that declare a `TASKNAME`, before #277:

| | Before | After |
|---|---:|---:|
| resolve a reference file directly | 66 | 66 |
| declared in `core.citations.NON_CITABLE` | 26 | 26 |
| **render an empty bibliography** | **52** | **0** |
| …of those, resolvable through the alias map | 52 | — |
| …of those, genuinely lacking a citation | 0 | — |

Every one of the 52 was empty for a single reason: the report side had no alias
map while the Bibliography button did. `phaser_MR_FRF` → `phaser`,
`xia2_dials` → `xia2`, all ten `crank2_*` → `crank2`, `imosflm` → `mosflm`,
`shelxeMR` → `shelxe`. The citations existed and were already correctly mapped;
the per-report bibliography just never consulted the mapping.

So it was not 52 missing reference files. It was one missing call — which is
why the warning count was worth publishing as a burn-down rather than
suppressing. The count went to zero in a single change, and the codebase now has
one key→file resolver instead of two.

Reproduce by walking `*_report.py` for `TASKNAME` and resolving each through
`CCP4Utils.findReferenceFile` and `core.citations.TASK_CITES`.

## What it says about the wider work

Three things this document argues, all visible in one afternoon:

1. **A check that exists is worth nothing until something enforces it.** The
   reference files had been wrong on Linux since they were added; a CI job on a
   case-sensitive runner found it on its first run.
2. **The first fix went the wrong way**, and could not have been known to be
   wrong on the machine it was written on. Sequencing the guard *before* the fix
   is what caught it.
3. **The silence was the real defect**, not the filename. An `errReport` nobody
   reads is the same shape as `processOutputFiles()` returning `FAILED` into a
   `process()` that discards it (C1), and as a slot exception logged and
   swallowed (C3).

---

# Best-practices amendments

`docs/pipeline_best_practices.md` needs updating once the core fixes land. Its
current advice is correct *given* the defects, which is itself a sign of the
problem:

| Section | Change | Blocked on |
|---|---|---|
| §5.1 "Dual Error Reporting" | The two-channel dance exists because `appendErrorReport` loses severity and `process()` ignores `processOutputFiles`. Replace with `self.fail()` and a single channel. | C1, C2, C4 |
| §5.2 "Severity Levels" | The table says WARNING → "Continue". False for `processInputFiles` / `makeCommandAndScript` / `startProcess`, where a WARNING fails the job. Correct once C5 lands. | C5 |
| §5.3 "Error Checking Pattern" | `>= 4` becomes the documented rule everywhere rather than one call site's convention. | C5 |
| §5.4 "Exception Wrapping" | Keep, but note that unwrapped exceptions in signal slots are now captured rather than lost. | C3 |
| §6.1 `processOutputFiles()` pattern | Show M1's declarative output contract instead of hand-rolled existence checks. | M1 |
| §7 "XML Management" | Add the M4 safe accessors. | M4 |
| §9.3 "Bare Except Clauses" | Already correct; add that it is now lint-enforced. | G3 |
| §9.4 "Debug Prints in Production" | Already correct; point at `self.log()`. | M6, G3 |
| new §13 | Log-failure tables (M2) and `runTimeValidity` (M7). | M2, M7 |
| new §14 | Declaring a task's external dependencies — `AUXILIARY_PROGRAMS`, the PHIL source declarations, and the missing-dependency policy. | M8 (landed) |

`mddocs/pipeline/ERROR_HANDLING_PATTERNS.md` is largely sound — its worked
example does pair `reportStatus(FAILED)` with `return CPluginScript.FAILED`, and
step 5 of its checklist says "Return early", so it is not the source of the C4
pattern. Two smaller amendments:

| Section | Change | Blocked on |
|---|---|---|
| "Key Elements" step 1 | Teaches `print()` for phase logging; becomes a lint error. Point at `self.log()`. | M6, G3 |
| "Key Elements" step 4 | `reportStatus(FAILED)` followed by `return FAILED` emits `finished` twice — once directly, once from `process()`. Replace both lines with `self.fail(102, …)`. | C4/M3 |

---

## Changelog

| Date | Change |
|---|---|
| 2026-08-24 | Document created from the audit of `b3a83af8a`; baseline counts recorded; nothing started. |
| 2026-08-24 | Added [Bibliography](#bibliography): a live instance of the silent-failure class, found by PR #276's new CI job and closed by #277. 52 report classes rendered an empty bibliography; all 52 had one cause and went to zero in a single change. Kept as a worked example. |
| 2026-08-25 | Pre-run program-availability check made **blocking where it is authoritative**. Severity is now decided by deployment: local execution against a mounted CCP4 (desktop/Electron, i2run) blocks submission, because a missing binary there is a certain failure; Azure mode (job queued for a worker whose filesystem we cannot see) and the slim CCP4-free API server stay advisory. Orthogonally, `TASKCOMMAND` only blocks for plugins that leave `process()` to the base class — `crank2` declares `crank2.py` but runs in-process, and `buster` declares `refine` but sources `$BUSTERDIR/setup.sh` to find it, so blocking either would break a working task. `AUXILIARY_PROGRAMS` is opt-in and always blocks when authoritative. Helper: `context_run.program_checks_are_authoritative()`. |
| 2026-08-25 | **M2 landed** for its motivating case. ARCIMBOLDO exits 0 after printing `FATAL` (`ARCIMBOLDO_LITE.main()` wraps the run in `except SystemExit: pass`), so a missing shelxe produced a job marked *Finished* with no outputs and no message — verified by rerunning the failing job's `setup.bor`: exit code 0, empty stderr, `FATAL` in the log only. Base `postProcessCheck` now applies `LOG_FAILURES`. Also added `AUXILIARY_PROGRAMS`, so binaries a task drives from *inside* `TASKCOMMAND` are covered by the pre-run availability check and listed on Preferences → Program locations; arcimboldo now resolves phaser/shelxe through `resolve_program()` instead of hardcoding `$CCP4/bin`, which is what made a correctly-configured `SHELXDIR` ineffective for this task. |
| 2026-08-25 | Re-measured at `59579c932` (3.1.0a28). Counts essentially flat; `processOutputFiles` returning `FAILED` **rose** 33 → 35, which is the argument for doing C1 next in one line. |
| 2026-08-25 | Added [A second axis](#a-second-axis--nothing-to-discard-in-the-first-place): the core silently rewriting what a wrapper was told, from the `setFullPath()` cross-project path bug (`043774c4c`) — no exception, no report, corruption written into `input_params.xml`, surfacing later as a true message about an invented path, attributed to the wrong layer. Includes the weakly-held-parent branch that silently does not run. No scan metric can see this class; it is a guardrails (G2) problem, and it raises C7's provenance work. |
| 2026-08-25 | **M5 amended** after `363d60e1c`. Keeping a cached rendering when regeneration fails is right, and the report cache is deliberately not version-stamped for the same reason. But the failure now reaches only a `logger.warning`, and the viewer sees an older rendering with nothing saying so — one silent failure traded for another. M5 gains: a preserved rendering must announce itself and say why regeneration fails. |
| 2026-08-25 | **M8 added** — declare the tool, let the base class find it, and say so when it is missing. Covers `AUXILIARY_PROGRAMS` and the `PHIL_SCOPE`/`PHIL_PARAMS_FILE`/`PHIL_PROGRAM` declarations (`3fd357f40`), which removed seven identical `get_master_phil()` bodies and settled a real disagreement: a missing tool used to stop the xia2 tasks loading at all, while the phasertng pair returned `None`. States the missing-dependency policy the rest of Part 2 should follow. |
| 2026-08-24 | C2 split into C2a/C2b after measuring that the proposed default would flip 268 of 320 calls, and confirming `appendErrorReport` severity does not gate job status. Added Phase 0 (baseline capture) and [The Phase 1 transition](#the-phase-1-transition); added `--predict-red-list` to the scan script. |
