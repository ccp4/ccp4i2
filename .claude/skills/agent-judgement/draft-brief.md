# Brief: draft a CCP4i2 task judgement file

You are writing knowledge an AI agent will use to drive CCP4i2 without a
crystallographer present: whether to run the task, how to set its key
inputs, whether the job worked, what to do next. A threshold that is
plausible but wrong is the worst thing you can write: the agent will act on
it with confidence. Every claim is checked in job files or attributed to a
named source; where you are unsure, the file says so.

## Read first
- The format: `docs/agentic-knowledge.md`, section 3. Follow it exactly.
- The judgement the format cannot hold: `.claude/skills/agent-judgement/SKILL.md`,
  last section. It records mistakes already made; do not repeat them.
- Two finished examples: `server/ccp4i2/pipelines/phaser_simple_phil/script/phaser_simple_phil.agent.yaml`
  and `server/ccp4i2/pipelines/molrep_pipe/script/molrep_pipe.agent.yaml`.

## The task
Given per task by the caller: its name, code directory, help page, and the
finished jobs to use as evidence (project directory and job numbers). Also:
- `python3 docs/user/tools/taskcard.py <task>`: parameters, interface, tests.
- The job tree: `sqlite3 <home>/db.sqlite3` (read-only): `ccp4i2_job`
  (number, task_name, status, parent_id), KPIs in `ccp4i2_jobfloatvalue`.
- The CCP4 documentation: `/Users/nmemn/Developer/ccp4-20260702/html/`.

## What to produce
- `<OUT>/<task>.agent.yaml`.
- Every `results` entry has an xpath (ElementTree subset; `[last()]` is
  allowed) and a `type`, and you have READ it from each evidence job's
  program.xml (or the named file) and seen the value you expect:
  `python3 -c "import xml.etree.ElementTree as ET; print(ET.parse('<job>/program.xml').find('<xpath>').text)"`.
- A number that exists only in a log: a `# FINDING:` comment naming it and
  where it is; do not point at the log.
- `basis` on every threshold: documentation (file and line), literature, or
  our runs (project, job, value). Never a value you have not read.
- The verdict scale fits the task: say what each outcome means for the next
  step. Next steps use the file-use syntax (`"[-1].XYZOUT[0]"`) and real
  parameter names from the next task's def.xml.
- 3-6 traps that really bite an agent, not a list of every option.
- `<OUT>/<task>.notes.md`: (a) every number used and where you read it;
  (b) anything in the help page or code you think is wrong, with file:line;
  (c) the claims you are least sure of, for an expert.

Validate before you finish: from `server/`, load your file with
`ccp4i2.agent.judgement.load(path=...)`, check `problems()` is empty, and
`judge(task, <job dir>, judgement=...)` each evidence job; report the
outcome each gets and whether that is right.

Do not run CCP4 programs or jobs, do not edit the repository, do not contact
outside services.
