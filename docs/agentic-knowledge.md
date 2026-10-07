# Agentic knowledge: task judgement, routes and the MCP facade

> **Status.** Design note, 2026-10-02; pilot in progress. It takes up three
> rows of the table in [driving-ccp4i2-as-an-agent.md](driving-ccp4i2-as-an-agent.md)
> §7: *expert knowledge* (a task and recipe catalogue), *evaluate and
> critique* (helpers that turn a job's results into pass or fail), and *tools*
> (a curated MCP facade). The stance is that of
> [NLP_JOB_CONSTRUCTION_DISCUSSION.md](NLP_JOB_CONSTRUCTION_DISCUSSION.md) §0:
> CCP4i2 is the substrate a user's agent drives, whichever agent that is. So
> the knowledge is written as plain data any model can read, and the code that
> serves it is a layer over the REST API, never a second way in.

## 1. What an agent lacks today

An agent can already create, set up, validate and run any job through the
REST API. What it cannot get from CCP4i2 is the judgement a crystallographer
brings:

- **which task** to run, and when a task is the wrong one;
- **which inputs matter** among a task's hundreds of parameters, and how to
  choose them;
- **whether it worked**: which numbers, in which file, mean success, and at
  what values;
- **what comes next**, given the result.

The user help (`docs/user/`) now states much of this, in prose, checked
against runs (the scenarios, `docs/user/tools/scenario_*.py`). This design
turns that into data, kept beside the code it describes, and serves it.

## 2. The four layers

| Layer | What | Written by | Where |
|---|---|---|---|
| Task card | Inputs, outputs, parameter labels, the i2run tests: facts | Generated from the source | `docs/user/tools/taskcard.py --json`, served by the facade |
| Judgement | When to use, what to set, how to read the result, traps, next steps | Drafted by a model, checked against runs, reviewed by an expert | `<task dir>/script/<task>.agent.yaml`, beside the def.xml |
| Route | A path through several tasks, with the decision points between them | From the scenarios | `server/ccp4i2/agent/routes/*.yaml` |
| Facade | MCP tools and resources over the REST API | Code | `server/ccp4i2/agent/` (optional extra) |

Generated facts are never copied into a hand-written file: a judgement file
names parameters and results; the card says what they are.

The judgement file sits beside the def.xml so that it ships in the wheel (the
facade serves it from an installed app), and so that a change to a task's
code and to its judgement are seen together in review. The help page and the
judgement file state the same knowledge for two readers; the page is prose,
the file is what a program can evaluate.

## 3. The judgement file

One per task. YAML, because experts will review it by reading it.

```yaml
task: phaser_simple_phil
status: draft                     # draft | reviewed
reviewed_by: []                   # names, once an expert has read it
purpose: >
  One sentence: what the task does, in crystallographic terms.

use_when:
  - A search model of at least ~30% sequence identity (or a predicted model)
    covers most of one copy of the molecule.
not_when:
  - text: The data are anomalous and no model exists.
    instead: phaser_ep_auto_phil

inputs:                           # only the ones that need judgement
  - param: inputData.ASUFILE      # container-relative path, as the card names it
    advice: >
      How to choose it, and what goes wrong if it is wrong.

results:                          # the numbers that decide, and where they are
  TFZ:
    file: program.xml             # program.xml | a job file named here | kpi
    xpath: .//Solution/TFZ        # ElementTree subset, first match; or a list
                                  # of xpaths, the first that reads wins
    type: float                   # float | int | str ("4.01%" reads as 4.01)
    meaning: Translation-function Z-score of the top solution.
  TNCS:
    file: program.xml
    xpath: .//Analysis/TNCS
    attribute: tNCS               # an attribute of the element, not its text
    type: str
    optional: true                # absent is normal: not listed as missing
  LLG:
    file: program.xml
    xpath: .//Solution/LLG
    type: float

verdict:                          # evaluated in order; the first that holds wins
  - when: TFZ >= 8 and LLG > 60
    outcome: solved
    basis: "<the source: program documentation, a paper, or a run and the value seen>"
  - when: TFZ >= 6
    outcome: ambiguous
    basis: "<...>"
  - when: true
    outcome: failed

traps:
  - Text of a mistake that produces a plausible but wrong result, and how to
    see it.

next:
  - when: outcome == "solved"
    task: servalcat_pipe
    inputs:
      XYZIN: "[-1].XYZOUT[0]"     # the shared argument and file-use syntax
  - when: outcome == "ambiguous"
    advice: Refine anyway and judge by R-free; ...

sources:                          # where each claim was checked
  - docs/user/source/tasks/phaser_simple_phil/index.rst
  - "run: MDM2 job 25 (scenario_mr_steps.py)"
```

Rules:

- **`when` is a tiny language**, not Python: result names, numbers, strings,
  `+ - * /`, `== != < <= > >=`, `and`, `or`, `not`, parentheses, `true`
  (so `RFREE_START - RFREE >= 0.02`). A result that could not be read makes
  any condition on it unknown, and an unknown condition never holds; to ask
  about absence itself, `NAME == null` (absent) and `NAME != null` (read). The facade
  evaluates it with its own parser, never `eval`.
- **Every threshold has a `basis`**: the program's documentation, the
  literature, or our own runs (project and job, with the value seen). A
  threshold with no basis is the one an expert must look at first, and the
  review is aimed there.
- **`results` must point at something that exists** in a real job of the
  task. A check (`agent/check.py`) reads each one from the scenario jobs and
  fails on a path that finds nothing; a number that exists only in a log is a
  finding, fixed by adding it to the task's program.xml.
- **Results come from what the job wrote**: program.xml, params.xml, the
  KPIs and the job's own data files. Never from the rendered report
  (report_xml.xml), which is a presentation of those, made only when someone
  opens it. When a program records its numbers only in text (MrBUMP's quick
  mode, results.txt), the wrapper adds them to program.xml.
- **A finished job is judged once.** The verdict is kept in the job
  directory as `judgement.json` with the version of the judgement (a hash of
  its file and of the evaluating code) and reused until that changes. It
  records which version of a judgement said what about the job; a job still
  running is judged afresh each time.
- **`status: draft` until an expert has read it.** The facade says so with
  every verdict it gives from a draft.

## 4. Routes

A route is a scenario lifted out of Python: steps (task plus inputs, in the
file-use syntax both i2run and `set_parameter` accept), each followed by the
verdicts that continue, branch or stop. A scenario run is the route's test.
The first is MR from merged data: ASU contents, MR, refinement, model
building.

## 5. The facade (MCP)

A thin adapter over the REST API, so the server's validation still decides
everything. Two ways in:

- **HTTP, served by the app itself at `/mcp/ccp4i2`**
  (`server/ccp4i2/agent/http.py`, mounted in `config/asgi.py`).
  - **Path:** scoped like the REST API (`/api/ccp4i2`), so an application
    serving CCP4i2 can serve MCP servers of its own beside it. This is an
    ASGI prefix, which a host cannot re-route as it can a URLconf, so
    `CCP4I2_MCP_PATH` moves it; an empty or `/` value keeps the default.
  - **Where it is served:** on the desktop; in a deployment only with
    `CCP4I2_MCP=1`, so a deployment opts in rather than gaining an endpoint
    by updating CCP4i2. `CCP4I2_MCP=0` turns it off anywhere.
  - **Stateless,** so either of the desktop's two uvicorn workers can
    answer.
  - **Authentication:** these requests skip Django's middleware, but each
    tool call goes back to the same server's REST API with the caller's own
    `Authorization`, so the deployment's authentication and permissions
    decide everything. The caller's address goes with it as
    `X-Forwarded-For`. `/mcp/ccp4i2` itself refuses a request with no
    `Authorization`; on the desktop it checks the session token.
  - **Connecting:** Help > About shows the address and a setup line. Both
    change at every launch, deliberately: an agent's access lasts as long as
    the session the person started.

        claude mcp add --transport http ccp4i2 http://127.0.0.1:<port>/mcp/ccp4i2 \
            --header "Authorization: Bearer <token>"

- **stdio**, `i2-mcp` (or `python -m ccp4i2.agent.mcp_server`), for clients
  that only speak stdio; `CCP4I2_URL` and `CCP4I2_TOKEN` say where the
  server is.

`mcp` is a runtime dependency. Its SDK derives a key with a call CCP4's
`cryptography` (2.8, not ours to replace) cannot make; `agent/request_state.py`
supplies the same derivation from the standard library through the SDK's own
codec hook. If the MCP server cannot be built at all, the app serves the REST
API alone and says so in its log: the agent route never stops the app.

Tools, first set:

| Tool | Over |
|---|---|
| `list_tasks`, `describe_task` | the registry, the card, the judgement file |
| `create_job`, `set_parameter`, `validate`, `run` | the job endpoints |
| `job_status`, `judge_job` | job tree; the judgement file's `results` read from the job's files, its `verdict` evaluated |
| `what_next` | the judgement file's `next`, then the server's `what_next` |
| `list_files` | the project's files and the jobs they came from |

Resources: the help pages, the judgement files, the routes. Deleting,
changing job status and anything outside the project are not offered.

## 6. Pilot and the test of it

1. This note and the judgement-file format.
2. Judgement files for two tasks drafted by two models (Fable 5.1, Opus 5.5)
   from the same brief, compared against the job files; the user chooses the
   model for the rest.
3. Judgement files for the MR route; the result check; the facade.
4. The test: a fresh agent on a less capable model, given only the facade,
   solves an MR case on a scratch server. If it can, the knowledge is in the
   material, not in the model. Where it stalls is the next thing to write.
