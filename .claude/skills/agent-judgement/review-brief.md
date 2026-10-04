# Brief: review a CCP4i2 task judgement file

Another model has drafted a judgement file: what an AI agent will use to
decide whether a CCP4i2 job worked and what to do next. Your job is to find
what is wrong with it before an agent acts on it. You are not rewriting it;
you are the check.

## Read first
- The format: `docs/agentic-knowledge.md`, section 3.
- The known mistakes: `.claude/skills/agent-judgement/SKILL.md`, last section.
- The draft and its notes (given by the caller), the task's code and help
  page, and the evidence jobs (given by the caller).

## Check, in this order
1. **Each result is the right number.** Read every xpath against every
   evidence job yourself. Is it the number the threshold is meant for
   (search vs refined; first vs last cycle; this copy vs all copies)? Is
   there a better source the draft missed, including elsewhere in
   program.xml or a sub-job's file?
2. **Each threshold's basis holds.** Open the cited documentation at the
   cited place; re-read each cited job value. Flag a basis that is memory
   presented as fact.
3. **The verdicts on the evidence jobs are right.** Judge each evidence job
   with `ccp4i2.agent.judgement.judge` (from `server/`); is the outcome what
   a crystallographer would say, given everything in the job (logs
   included)? Is the verdict scale fine enough for the next step to differ?
4. **The traps and next steps.** Would an agent following them go wrong?
   Do the parameter names and file references exist?
5. **The code.** Anything in the task's code or help page that is wrong and
   that the draft missed; cite file:line and say how you verified it.

## What to produce
`<OUT>/<task>.review.md`: findings, each with severity (wrong / doubtful /
improvement), the evidence (file, xpath, value, line), and the change you
propose. Say also what you checked and found right, briefly, so the
resolver knows what was covered. Do not edit the draft or the repository;
do not run CCP4 programs or jobs; do not contact outside services.
