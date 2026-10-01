# CCP4i2 user documentation

The task help, as Sphinx/RST. Imported from
[gitlab.com/ccp4i2/rstdocs](https://gitlab.com/ccp4i2/rstdocs), which stays the
help for the Qt interface; this copy is the help for this one.

## Building

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r docs/user/requirements.txt
.venv-docs/bin/sphinx-build -b html docs/user/source /tmp/userdocs
```

## Pictures from the running app

A page whose figures come from the app has a `shots.json` next to its
`index.rst`, naming what to capture, and a scenario script in `tools/` that
builds the project with real jobs the pictures are taken from. Nothing is
added to `demo_data`; the scenario fetches its data from the PDB and PDB-REDO.

```bash
# 1. The scenario, into a scratch home (never a live one: the script refuses).
cd server
env CCP4I2_HOME=/tmp/docs-home DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
    ccp4-python manage.py migrate
env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_parrot.py

# 2. The app against that home, on ports of its own (Django 3421, Next 3420).
CCP4I2_HOME=/tmp/docs-home docs/user/tools/devserver.sh start
#    ... restart django   after a server change (it runs --noreload)
#    ... clear-reports <project> <job>...   after a report class changes
#    ... stop

# 3. The pictures (Node 22+, Google Chrome; CHROME=path if not the macOS default),
#    from a shell that has not sourced CCP4 (its setup puts Node 20 first).
cd docs/user/source/tasks/parrot && node ../../../tools/capture.mjs shots.json
```

`capture.mjs` opens each job page, turns developer mode off, crops to one
section and numbers the fields the text refers to. It finds fields by the label
the user sees, so when a label changes the capture fails and names it: that is
the moment to reread the page's text too. A shot with `"outline": true` prints
the page as text instead (tabs, sections, fields and their values, report
tables): the cheap way to learn what a new page's sections and labels are.
Captures run one at a time: each waits on a lock (`ccp4i2-capture.lock` in
the temporary directory) while another runs, since several at once against
one dev server slowed every shot five-fold.

To learn about a task before writing its page:

```bash
python3 docs/user/tools/taskcard.py <task>      # registry, def.xml, interface, tests, page state
cd server && env CCP4I2_HOME=... DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
    ccp4-python manage.py list_jobs <project>          # the scenario's jobs
    ccp4-python manage.py list_files --project <project> --job <n> --format json
```

Which pages are converted, and which still show the Qt interface, is on the
generated status page (`tools/status.py`); `tools/status.py todo` lists what
is left, by the chooser's categories.

## Keeping a page true: stamps

A page describes its task as it was when it was checked. Each `shots.json`
records, under `"sources"`, the git blob id of every file the page was checked
against: the task's def.xml, script and report (from `core/tasks.py`), its
interface and what that imports from `task-interfaces/`, and the scenario.
`tools/stamp.py` finds these files itself; nothing is listed by hand.

```bash
python3 docs/user/tools/stamp.py check          # which pages are stale, and why
python3 docs/user/tools/stamp.py diff parrot    # what changed since it was stamped
python3 docs/user/tools/stamp.py stamp parrot   # it still holds (or was re-captured)
```

`check` also reports any label a shot is keyed on (section, callout, from/until)
that no longer appears in those files, which would break the next capture.

A pull request that touches a stamped file runs the user-docs CI job, which
warns on the PR, naming the page and the files. It never fails: many changes
leave the page right. Read the diff, then re-capture or restamp in the same
pull request. The status page marks stale pages too.
