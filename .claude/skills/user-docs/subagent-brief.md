# Brief for a page-drafting subagent

The main session's prompt is short: "Read
`.claude/skills/user-docs/subagent-brief.md` and follow it", then the
specifics below. Everything else is here, so every subagent works the same
way and the prompts cost little. One task (or one shared page) per
subagent. The main session has already run the scenario, fixed what the
outlines showed, and checked the facts it gives.

**Specifics the prompt gives**

- the task; the worktree and branch; the scratch home (`CCP4I2_HOME`);
- the project, the run's job number(s), the unrun clone's job number;
- the page directory, and whether it is a Qt page to rewrite or a new page
  (and then its category in `tasks/index.rst`);
- the scenario (`../../../tools/scenario_<route>.py`);
- the checked facts, each with where it came from;
- known gaps not to document as features; neighbours to compare with.

---

You are drafting a CCP4i2 user-help page. Read
`.claude/skills/user-docs/SKILL.md` and `traps.md` first; follow them. Read
`docs/user/README.md` for `capture.mjs` and `shots.json`, and use
`docs/user/source/tasks/import_files/` and `.../pairef/` as models.

**What you have**

- The task card: `python3 docs/user/tools/taskcard.py <task>`.
- The app, running against the scratch home (Django 3421, Next 3420). Do
  not start, stop or restart it. If every shot says "Nothing to click
  labelled View", it is down: stop and report.
- The jobs, in `$CCP4I2_HOME/projects/<project>/CCP4_JOBS/job_<n>/`.
  Outline them first (`"outline": true` shots): cheap, and it gives the
  sections, labels and values the page and the shots need.
- The checked facts. Check them too: a fact in a brief can be wrong (a His
  tag from another project's sequence reached a brief in the first trial;
  the subagent caught it in the job's files). State nothing else as a
  result of the run without checking it in the job's files.

**What to do**

1. Write `shots.json` (`"draft": true`, the scenario, the project) and
   capture from a shell that has not sourced CCP4:
   `cd docs/user/source/tasks/<dir> && node ../../../tools/capture.mjs shots.json`.
   Captures queue on a lock, one at a time across all subagents, so a wait
   ("Waiting for another capture") is normal. Typically one input shot
   (the clone) and one report shot (the run), with numbered callouts.
2. Look at every picture (Read the PNG). Fix the shots, not the app.
3. Write the page: what the task is for and when to use it rather than its
   neighbours; the input, with callouts; the results, with what the numbers
   mean and what to do next. Keep the old page's prose where still true.
   Delete old pictures the page no longer uses, and any stray `index.html`.
4. `python3 docs/user/tools/compress_images.py docs/user/source/tasks/<dir>`,
   then `python3 docs/user/tools/stamp.py stamp <dir>`.

**Do not** change any file outside the page directory (and the index entry
for a new page); do not run i2run, scenarios, the test suites or the Sphinx
build; do not commit or touch git; leave no process running when you
report.

**Report back**, briefly:
- the files written;
- every defect you saw in the app (raw parameter names, unannotated or
  mislabelled outputs, layout faults, report errors), with the picture or
  outline line that shows it;
- every statement on the page that is not in the checked facts, with its
  source, saying which you could not check.
