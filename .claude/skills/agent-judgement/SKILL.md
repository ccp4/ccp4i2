---
name: agent-judgement
description: Write or revise a CCP4i2 task's judgement file (<task>.agent.yaml: when to use the task, which inputs need thought, where its deciding numbers are and what they mean, traps, next steps), or change the agent facade that serves them (server/ccp4i2/agent/, the MCP server i2-mcp, the jobs/<id>/judgement and parameters endpoints). Use when making CCP4i2 drivable by an AI agent.
---

# Task judgement for agents

The design is [docs/agentic-knowledge.md](../../../docs/agentic-knowledge.md);
the file format is its section 3. This skill holds how the files are made
and the judgement the format cannot express.

## Where the truth lives

| What | Where |
|---|---|
| A task's judgement | `<task dir>/script/<task>.agent.yaml`, beside its def.xml |
| Reading results, the `when` language, the verdict | `server/ccp4i2/agent/judgement.py`, `condition.py` |
| A job's parameters, one line each | `server/ccp4i2/agent/parameters.py`; `GET jobs/<id>/parameters/` |
| The MCP tools | `server/ccp4i2/agent/mcp_server.py` (`i2-mcp`, extra `agent`) |
| Is a file well formed | `tests/unit/lib/test_agent_judgement.py` (checks every shipped file) |
| Real jobs to read numbers from | the help scenarios' projects (`docs/user/tools/scenario_*.py`); `manage.py list_jobs` |
| What the task is, from source | `python3 docs/user/tools/taskcard.py <task>`; its help page |

## How a file is made: one model drafts, another reviews

Settled on 2026-10-02 by drafting the same two files (Phaser, MolRep) with
two models. Each made one serious error the other did not: one judged
Phaser on the *refined* TFZ (lenient; Phaser's 8 applies to the search TFZ)
and misread a database; the other gave a coarser verdict scale and missed
code defects. So:

1. **Draft** with Opus: [draft-brief.md](draft-brief.md).
2. **Review** with Fable: [review-brief.md](review-brief.md). The reviewer
   checks every result against the job files and every claim against its
   basis, and writes findings, not a rewrite.
3. **Resolve** in the main session: check each disagreement in the job files
   yourself, apply, run the test, judge the evidence jobs. Then commit as
   `status: draft`.

An expert's reading is what turns `draft` into `reviewed`; the facade says
"draft" on every verdict until then.

## Judgement the format cannot hold

- **Judge on the number the threshold was made for.** Phaser's TFZ 8 is for
  the search (`Placement/TFZ`), not the refined solution (`Solution/TFZ`).
  MOLREP's documented contrast bands (>3 definite) do not match the contrast
  column it prints (20-55 for right and wrong alike). A threshold is only as
  good as the number it is applied to.
- **A finished job is not a successful one,** and a program's own "solution
  found" is not proof: the wrong enantiomorph passes every MOLREP check;
  only refinement (R-free falling to a sensible value) separated the hands.
- **Distinguish placed from solved.** A correct placement of a distant or
  partial model ends at R-free 0.45-0.50; the next step is building, not
  another MR. Verdict scales are task-specific; say what each outcome means
  for the next step.
- **A number only in a log is a finding.** Fix the wrapper to put it in
  program.xml (as molrep_mr now records MOLREP's score and z-score), with a
  test, rather than pointing an agent at a log.
- **Placeholders are worse than nothing.** A field written with a constant
  (`mr_score 0.0000`) reads as a result. Leave a field out if the program did
  not produce it; the verdict then sees it as missing.
- **Thresholds from memory are marked as such** in their `basis`, and are
  the first thing an expert is asked to check.
- Declare every result's `type` (float, int, str): an untyped text value is
  read as a float, fails, and is silently missing.
- **Never name an interactive task as a `next` `task:`** (coot_rebuild,
  moorhen, coot1; anything with `interactive=True` or a GUI an agent cannot
  drive). Say in `advice:` what a person should do there, and give the
  non-interactive task that follows, if there is one.
- **One vocabulary per kind of task,** so an agent reasons the same way
  whichever program ran. Refinement (prosmart_refmac, servalcat_pipe):
  `refined` (done for now), `unconverged` (same again, more cycles),
  `not_improving` (R-free did not fall: something upstream is wrong),
  `needs_building` (fell, below 0.40: build next), `placed` (fell, at or
  above 0.40: a right but distant or partial model, or a wrong one; build
  decides), `failed`, `unjudged` (nothing to judge on). MR: `solved`,
  `placed`, `partial`, `ambiguous`, `failed`. Experimental phasing (crank2,
  shelx, phaser_ep_phil): `built` (phased, a model built and refined:
  refine or validate next), `phased` (substructure and hand decided, map
  interpretable, little or nothing built: build next), `ambiguous`
  (substructure found but the hand or the map unconvincing), `failed`.
  Density modification: `improved`, `no_better`, `failed`.
- **A pipeline that stops early is not a pipeline that finished.** Say
  which step's output to look at, and what "finished" means for the model
  (built and refined, or only phased). An agent trial (Haiku, 2026-10-03)
  phased Gamma with Crank2, stopped, and called the model ready to deposit.

