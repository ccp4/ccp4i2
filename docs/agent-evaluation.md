# Evaluating agents that drive CCP4i2

> **Status.** Design note, 2026-10-06. Nothing here is built yet. It
> follows [agentic-knowledge.md](agentic-knowledge.md), which describes what
> an agent is given: per-task judgements, routes and the MCP facade. This
> note covers how to tell whether a change to those makes agents better
> across the board, not just on the structure that prompted it.

## 1. The problem: one dataset at a time

Until now the material has been tuned one trial at a time. An agent solves,
or fails to solve, one structure; we read what it did; we fix what misled
it. The week of CDK4/cyclin D1 (2026-10-05) shows both what this finds and
what it costs.

- Haiku placed CDK4 against a cyclin chain cut from another crystal's entry
  (6p8e), given to Phaser as a structure *already placed*. R-free stalled at
  0.42.
- Reading its transcript found three causes, none of them the one first
  guessed:
  - `list_tasks("phaser molrep molecular replacement")` returned nothing,
    because the query was matched as one phrase;
  - mrparse's judgement offered `phaser_simple_phil` first to every job.
    "One kind of molecule" was prose in the advice, not a condition;
  - nothing refused a fixed structure from another crystal.
- Each fix was checked on CDK4's files, and on nothing else.

Every fix to a judgement, a validity check or a facade tool may change what
an agent does on structures we are not looking at. A fix for dataset B can
quietly break dataset A, and the only check today is another expensive live
trial, on one dataset.

What is needed is a fixed set of scenarios that every change is measured
against, cheap enough to run on every change, and a routine for turning
failures into fixes that stay fixed.

## 2. What is being tested

Three things, which fail differently and are best tested separately.

| Layer | Example of a failure | Deterministic? |
|---|---|---|
| **The material**: judgements, validity checks | mrparse offers the single-model task for a complex | Yes: a pure function of a job's files |
| **The facade**: MCP tools | `list_tasks` finds nothing for a natural query | Yes, given the server's answers |
| **The agent's behaviour** given both | Haiku fills the fixed-structure slot with the second component | No: a model choosing, differently each run |

The first two can be tested like any other code. Only the third needs agent
runs, and those runs need not execute any crystallographic program, as §4
explains.

## 3. Tier 1: decision fixtures (no agent, no CCP4 programs)

Judgements read only what a job wrote: `program.xml`, `params.xml`,
`diagnostic.xml` and the job's KPIs, never `report_xml.xml`
(agentic-knowledge.md). Validity checks read a job's parameters and, at
submission, the headers of its input files. So a job's conclusions can be
recomputed in milliseconds from a snapshot of its directory, without
running anything.

**A fixture** is one scenario's jobs, recorded once from real runs. For
each job it keeps the files the judgement and validity checks read. The
large binary outputs are dropped. Input files are kept only as far as the
checks need them: an MTZ cut down to its header and a few reflections still
answers a cell check.

```
agent-fixtures/<scenario>/
  scenario.yaml            # what it is, its source, which entries to hide (§6)
  jobs/03_mrparse/         # params.xml, program.xml, diagnostic.xml, kpis.json
  jobs/05_phaser_simple_phil/
  inputs/                  # trimmed input files the validity checks open
  expect.yaml              # what each job should conclude
```

**The assertions**, in `expect.yaml`, state what a crystallographer would
check:

- each finished job's verdict outcome;
- the tasks its next steps offer, in order, and any that must not appear.
  For CDK4's mrparse: `phaser_pipeline_phil` first, `phaser_simple_phil`
  absent;
- for prepared-but-unrun jobs, the validity codes they must raise or must
  not. For CDK4: a Phaser job with 6p8e_B as the fixed structure raises
  210, a blocking error.

These are the "decision tables" that settle the regression worry directly.
A change to a judgement either keeps every expected decision across every
scenario or names the ones it changed.

**Where they come from.**
- The scripted help scenarios (`docs/user/tools/scenario_*.py`) already
  build real projects through i2run: MR, MR steps, experimental phasing,
  SHELX, predicted models, data reduction, refinement, ligands and more.
  Their job directories are the first fixtures.
- The held-out trials (rnase, cdk1 ternary, HypF Hg-SAD) add more.
- So does the CDK4/cyclin D1 project.

**Changing expectations.** When a judgement is changed on purpose, the run
lists every decision that moved. The reviewer accepts them by updating
`expect.yaml` in the same PR. A decision that moves without anyone meaning
it to is the regression this tier exists to catch.

**Where they run.** The judgement engine and the validity checks are
CCP4-free, so this tier joins the CI unit suites. It needs seconds, not
hours.

## 4. Tier 2: agent runs on replayed execution

Tier 1 cannot say whether an agent will *find* the right step, as the empty
task search showed. That takes an agent driving the real facade. What it
does not take is real crystallographic computing: a job whose inputs have
been run before can return what it returned then.

### 4.1 Why not record the MCP responses

The obvious harness records each tool call and its response, and replays
the response when the same call comes again. This breaks at the first fix.
A fix exists to make the agent do something different, and what it does
next was never recorded. After this week's changes Haiku would call
`phaser_pipeline_phil`, and no recording covers it. Replay keyed on request
bytes tests only that nothing changed.

### 4.2 Replay below the facade, keyed on meaning

So the harness replays at the one point where time is spent: running a job.
Everything above it is real:
- the Django server and its database (in a scratch home, never a live one);
- the REST API and the MCP facade;
- validation, filling in inputs from earlier jobs, and gleaning outputs;
- judgements.

**A replay job target.** The run-target registry (`ccp4i2/lib/dispatch/`)
already lets a deployment add a job target by naming its class in settings
(`CCP4I2_RUN_TARGETS`, `CCP4I2_JOB_TARGET`). The harness registers a
`replay` target in its own settings overlay, so no code in ccp4i2 knows it
exists. Its `run_job`:

1. computes the job's **key** from:
   - the task name;
   - a content hash of each input file, computed by the target. CCP4i2
     records a checksum only for imported files (`FileImport.checksum`),
     not for job outputs;
   - the parameters that were set, normalised (sorted, defaults dropped).
2. on a **hit**, copies the recorded job's outputs into the new job's
   directory, then lets the normal finishing path run: gleaning registers
   the output files, the KPIs are stored, and the status becomes the
   recorded status. The agent sees a real job that finished instantly.
3. on a **miss**, follows the run's policy:
   - *regression* runs stop the scenario and record where, and with what
     key, the agent left the recordings. That location is itself a finding:
     either a new route the agent took, or a fix that sent it somewhere
     unrecorded.
   - *exploration* runs execute the job for real (the `local` target),
     within a budget of CCP4 hours, and add it to the recordings.

Chained jobs keep their keys stable. A replayed job's outputs are the
recorded files, byte for byte, so the next job's key matches the recording
made after that same job.

A recording is one run of each program, so randomness in a program's own
search (Phaser, ModelCraft) is frozen. That is acceptable for testing the
agent. The occasional live run (§4.4) keeps the recordings honest.

### 4.3 Running the agent

- **The driver.** The agent is driven headless (the Claude Agent SDK, or
  `claude -p` with an MCP configuration pointing at the harness server).
  Its brief is the scenario's prompt and nothing else, the way the held-out
  trials were run. The transcript is saved; this week's diagnosis came from
  one.
- **The model.** By default the weakest model we intend to support (Haiku
  today). Material that steers Haiku right steers stronger models too, but
  not the other way round.
- **Repetition.** Each scenario runs 3–5 times. An agent is random: one run
  that passes, or one that fails, says little. Results are pass *rates*.

### 4.4 Tier 3: live runs

A small rotating subset (two or three scenarios a night) runs with real
execution throughout. This catches recordings that no longer match what the
programs do: a new CCP4 build, a wrapper change.

## 5. Grading

**The outcome** is checked against the deposited structure, automatically:

- the space group, or at least the point group, matches;
- every component the deposited model has was placed, matched by sequence
  to its chains;
- the final R-free is within a margin of the deposited value, or below a
  fixed bound when the deposited value is not comparable;
- the run finished, and did not give up.

**The process** is checked by rules over the transcript and the project's
jobs:

- validity errors hit, and whether the agent then corrected the job;
- tool errors, and tool calls returning nothing (an empty `list_tasks`);
- known misuses (§7): a homologue as a fixed structure, two components in
  one ensemble, a stage-by-stage route where a pipeline exists;
- the number of jobs, and repeated identical jobs.

A model grading the transcript against a rubric can come later, for things
rules cannot see: whether the agent said why it chose what it did, or
reported honestly. Rules come first, because they are cheap, repeatable and
cannot be argued with.

**The report** gives each scenario's pass rate and the misuses counted. It
also compares the totals with the last accepted run. A change is accepted
when tier 1 keeps every decision it did not mean to change, and the tier-2
totals do not fall. One scenario improving is not enough.

## 6. The scenario corpus

### 6.1 Coverage, not volume

About **50** scenarios to start, chosen to cover the kinds of problem, not
picked at random. A thousand would mostly be copies of easy molecular
replacement. Proposed coverage:

| Kind | About | Notes |
|---|---|---|
| MR, close homologue (> 60%) | 6 | including several copies in the AU |
| MR, distant homologue (25–40%) | 6 | where trimming (sculptor) matters |
| MR, predicted model only | 5 | AlphaFold; domain splitting (slicendice) |
| MR, complexes | 6 | two to four components, as CDK4/cyclin D1 |
| Nucleic acid, or protein with nucleic acid | 3 | |
| SAD / MAD | 6 | Se-Met, S, metals, halides; the HypF Hg case |
| Twinned or pseudo-symmetric | 4 | |
| Low resolution (> 3.2 Å) | 4 | |
| Ligands | 4 | where the ligand changes the route |
| Data reduction from unmerged data | 4 | §6.3 |
| Held out | ~10 | §6.4 |

The corpus grows when a real run fails in a way the corpus does not yet
contain, as CDK4 did this week.

### 6.2 Keeping the answer out of the search

A scenario made from a PDB entry has its answer in the PDB. mrparse would
find the target entry itself, at 100% identity, from the same crystal. Its
command-line options (`ccp4-20260702`) include none that excludes entries or
sets a date cutoff, so the exclusion has to be done by CCP4i2 or the harness:

- each `scenario.yaml` names the entries to hide:
  - the target entry;
  - other entries of the same crystal form;
  - for a fair test, entries deposited after it;
- a recorded mrparse job has those hits removed before it is stored, so
  replay never offers them;
- a live mrparse run needs the same filter in the wrapper. That means a
  parameter listing entries to drop from `homologs.json` before models are
  registered. It is a small change to `mrparse_wrapper.py`, and useful
  beyond evaluation (re-solving a known structure as a teaching exercise).

The agent's own memory cannot be filtered. A brief never names the entry,
and the facade's rule against fetching structures from outside CCP4i2
stands. Famous structures may still be recognised; scenarios where the
model is suspiciously quick deserve a look.

### 6.3 Data reduction needs unmerged data

Most PDB entries deposit merged data. Data-reduction scenarios therefore
come from:
- the demo data that ships unmerged (gamma, mdm2, the help scenarios'
  reduction projects);
- the few entries that deposit unmerged intensities;
- data of our own, such as CDK4/cyclin D1, which stay local. Nothing from
  such a scenario goes to an outside service beyond what the tasks
  themselves send (mrparse sends the sequence to EBI).

### 6.4 A held-out set

About ten scenarios are never looked at while the material is being edited.
They are run only to report. Otherwise the judgements come to fit the
suite, which is whack-a-mole at a larger scale. The held-out trials of
2026-10-03 did this informally; here it becomes a rule.

## 7. The loop

1. Run the suite: tier 1 always, tier 2 on the scenarios a change could
   touch, or all of them before a release.
2. **Read failures in bulk, then group them** by what went wrong, not by
   scenario: "empty task search", "homologue given as fixed", "stalled
   refinement accepted". Count each group.
3. Fix the largest group, at the right layer: a facade tool, a validity
   check, a judgement condition, or advice.
4. **Lock it in tier 1.** Every fix adds the decision it corrects as an
   expectation on a fixture, so the same failure cannot come back
   unnoticed.
5. Rerun. Accept on the totals (§5).

Step 2 is where this week went wrong. One transcript, read early, would
have shown the empty task search; the first explanation was guessed from
the judgement instead. Across fifty transcripts, groups make the commonest
failure obvious.

## 8. Cost

- **Recording:** one real run of each scenario's route. Exploration runs
  add more as agents find new routes, within a budget of CCP4 hours.
- **Regression:** agent tokens only. At Haiku's prices, 50 scenarios × 3
  runs is affordable per change. Job time is negligible on replay; what
  remains is the agent's own thinking.
- **Storage:** recordings are job directories without large binaries. Like
  the i2run tier's downloads (`~/.cache/ccp4i2-tests/downloads`), they are
  kept outside the repository, fetched by a script, and versioned by the
  CCP4 build they were recorded on. Tier-1 fixtures are small enough for
  the repository.

## 9. Plan

1. **Tier 1 on what we have.** Fixtures from the help scenarios, the
   held-out trials and CDK4/cyclin D1, about 10 scenarios. Expectations for
   every judgement decision they contain, plus the CDK4 validity cases. Run
   in CI.
2. **The replay target and the driver.** The settings overlay, the key, the
   hit and miss paths, and the headless driver with saved transcripts, on
   the same 10 scenarios.
3. **Grading and the report.** The outcome checks against deposited models;
   the process rules; comparison with the last accepted run.
4. **The corpus to ~50,** with the mrparse exclusion parameter and the
   held-out set.
5. **Nightly live runs** of a rotating subset.

## 10. Open questions

- **R-free margin.** How close to the deposited R-free counts as solved?
  It should vary with resolution.
- **Which models to measure.** Haiku only, or Sonnet as well, so we know
  whether a change helps the weakest model at the strongest one's expense?
- **Model-graded rubric.** Worth it, and for which judgements?
- **Recordings and CCP4 builds.** Recordings made on one build may not match
  another's programs. Re-record a build's whole corpus when we validate it,
  as the i2run baselines are re-run now?
