# Brief for a page-drafting subagent

Fill in the brackets and send as the subagent's prompt. One task (or one
shared page) per subagent. The main session has already run the scenario,
fixed what the outlines showed, and checked the numbers below.

---

You are drafting the CCP4i2 user-help page for **[task]** in the worktree
`[path]`, branch `[branch]`. Read `.claude/skills/user-docs/SKILL.md` and
`traps.md` first; follow them.

**What you have**

- Task card: run `python3 docs/user/tools/taskcard.py [task]`.
- The app is running against the scratch home `[CCP4I2_HOME]` (Django 3421,
  Next 3420). Do not start, stop or restart it.
- Project **[project]**: the run is job **[n]**, its unrun clone job **[m]**.
  Outlines of both are in `[outline files]`.
- Page directory: `docs/user/source/tasks/[dir]/` ([Qt page to rewrite |
  new page; add it to `tasks/index.rst` under [category]]).
- Checked facts you may state: [the numbers and conclusions, with where
  each came from, and the job directory's path]. Check them too: a fact in
  a brief can be wrong (a His tag from another project's sequence reached
  a brief in the first trial; the subagent caught it in the job's files).
  State nothing else as a result of this run without checking it in the
  job's files, and say in your report what you checked.

**What to do**

1. Write `shots.json` (draft: true; scenario `[scenario]`; project
   `[project]`) and capture from a shell that has not sourced CCP4:
   `cd docs/user/source/tasks/[dir] && node ../../../tools/capture.mjs shots.json`.
2. Look at every picture. Fix the shots, not the app.
3. Write the page: what the task is for and when to use it rather than its
   neighbours; the input, with callouts; the results, with what the numbers
   mean and what to do next. Keep the old page's prose where still true.
4. `python3 docs/user/tools/compress_images.py docs/user/source/tasks/[dir]`,
   then `python3 docs/user/tools/stamp.py stamp [dir]`.

Captures are slow when several run at once (ten minutes a shot with five
in parallel, two alone): outline first, capture each shot once it is
right, and leave no capture or other process running when you report.

**Do not** change any file outside `docs/user/source/tasks/[dir]/` and the
index entry; do not run i2run, scenarios or the test suites; do not commit.

**Report back**, briefly:
- the files written;
- every defect you saw in the app (raw parameter names, unannotated or
  mislabelled outputs, layout faults, report errors), with the picture or
  outline line that shows it;
- every statement on the page that is not in the checked facts above, with
  its source.
